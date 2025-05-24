import os
import sys
import time
from typing import Dict, Optional
import shutil

import numpy as np
import pandas as pd
import pyranges as pr
import scipy.sparse as ssp
import hicstraw


def df_to_pyranges(
    df,
    start_col="start",
    end_col="end",
    chr_col="chr",
    start_slop=0,
    end_slop=0,
    chrom_sizes_map: Optional[Dict[str, int]] = None,
):
    df["Chromosome"] = df[chr_col]
    df["Start"] = df[start_col] - start_slop
    df["End"] = df[end_col] + end_slop
    if start_slop or end_slop:
        assert chrom_sizes_map, "Must pass in chrom_sizes_map if using slop"
        df["chr_sizes"] = df[chr_col].apply(lambda x: chrom_sizes_map[x])
        df["Start"] = df["Start"].apply(lambda x: max(x, 0))
        df["End"] = df[["End", "chr_sizes"]].min(axis=1)
        df.drop("chr_sizes", axis=1, inplace=True)
    return pr.PyRanges(df)


def get_hic_file(chromosome, hic_dir, allow_vc=True, hic_type="juicebox"):
    if hic_type == "juicebox":
        is_vc = False
        filetypes = ["KR", "INTERSCALE"]
        for filetype in filetypes:
            hic_file = os.path.join(
                hic_dir, chromosome, chromosome + f".{filetype}observed.gz"
            )
            hic_norm = os.path.join(
                hic_dir, chromosome, chromosome + f".{filetype}norm.gz"
            )
            if hic_exists(hic_file):
                print("Using: " + hic_file)
                return hic_file, hic_norm, False

        if allow_vc:
            hic_file = os.path.join(hic_dir, chromosome, chromosome + ".VCobserved.gz")
            hic_norm = os.path.join(hic_dir, chromosome, chromosome + ".VCnorm.gz")
            if hic_exists(hic_file):
                print(
                    f"Could not find KR normalized hic file. Using VC normalized hic file: {hic_file}"
                )
                return hic_file, hic_norm, True

        raise RuntimeError(
            f"Could not find {', '.join(filetypes)} or VC normalized hic files"
        )

    elif hic_type == "bedpe":
        hic_file = os.path.join(hic_dir, chromosome, chromosome + ".bedpe.gz")
        return hic_file, None, None
    elif hic_type == "avg":
        hic_file = os.path.join(hic_dir, chromosome, chromosome + ".bed.gz")
        return hic_file, None, None


def hic_exists(file):
    if not os.path.exists(file):
        return False
    elif file.endswith("gz"):
        # gzip file still have some size. This is a hack
        return os.path.getsize(file) > 100
    else:
        return os.path.getsize(file) > 0


def load_hic_juicebox(
    hic_file,
    hic_norm_file,
    hic_is_vc,
    hic_resolution,
    tss_hic_contribution,
    window,
    min_window,
    gamma,
    scale=None,
    interpolate_nan=True,
    apply_diagonal_bin_correction=True,
):
    print("Loading HiC Juicebox")
    HiC_sparse_mat = hic_to_sparse(hic_file, hic_norm_file, hic_resolution)
    HiC = process_hic(
        hic_mat=HiC_sparse_mat,
        hic_norm_file=hic_norm_file,
        hic_is_vc=hic_is_vc,
        resolution=hic_resolution,
        tss_hic_contribution=tss_hic_contribution,
        window=window,
        min_window=min_window,
        gamma=gamma,
        interpolate_nan=interpolate_nan,
        apply_diagonal_bin_correction=apply_diagonal_bin_correction,
        scale=scale,
    )
    return HiC


def load_hic_bedpe(hic_file):
    print("Loading HiC bedpe")
    return pd.read_csv(
        hic_file,
        sep="\t",
        names=["chr1", "x1", "x2", "chr2", "y1", "y2", "name", "hic_contact"],
    )


def load_hic_avg(hic_file, hic_resolution):
    print("Loading HiC avg")
    cols = {"x1": np.int64, "x2": np.int64, "hic_contact": np.float64}
    HiC = pd.read_csv(
        hic_file, sep="\t", names=cols.keys(), usecols=cols.keys(), dtype=cols
    )
    HiC["x1"] = np.floor(HiC["x1"] / hic_resolution).astype(int)
    HiC["x2"] = np.floor(HiC["x2"] / hic_resolution).astype(int)
    HiC.rename(columns={"x1": "bin1", "x2": "bin2"}, inplace=True)
    return HiC


# def juicebox_to_bedpe(hic, chromosome, resolution):
#     hic['chr'] = chromosome
#     hic['x1'] = hic['bin1'] * resolution
#     hic['x2'] = (hic['bin1'] + 1) * resolution
#     hic['y1'] = hic['bin2'] * resolution
#     hic['y2'] = (hic['bin2'] + 1) * resolution

#     return(hic)


def process_hic(
    hic_mat,
    hic_norm_file,
    hic_is_vc,
    resolution,
    tss_hic_contribution,
    window,
    min_window=0,
    hic_is_doubly_stochastic=False,
    apply_diagonal_bin_correction=True,
    interpolate_nan=True,
    gamma=None,
    kr_cutoff=0.25,
    scale=None,
):
    # Make doubly stochastic.
    # Juicer produces a matrix with constant row/column sums. But sum is not 1 and is variable across chromosomes
    t = time.time()

    if not hic_is_doubly_stochastic and not hic_is_vc:
        # Any row with Nan in it will sum to nan
        # So need to calculate sum excluding nan
        temp = hic_mat
        temp.data = np.nan_to_num(temp.data, copy=False)
        sums = temp.sum(axis=0)
        sums = sums[~np.isnan(sums)]
        # assert(np.max(sums[sums > 0])/np.min(sums[sums > 0]) < 1.001)
        mean_sum = np.mean(sums[sums > 0])

        if abs(mean_sum - 1) < 0.001:
            print(
                "HiC Matrix has row sums of {}, continuing without making doubly stochastic".format(
                    mean_sum
                )
            )
        else:
            print(
                "HiC Matrix has row sums of {}, making doubly stochastic...".format(
                    mean_sum
                )
            )
            hic_mat = hic_mat.multiply(1 / mean_sum)

    # Adjust diagonal of matrix based on neighboring bins
    # First and last rows need to be treated differently
    if apply_diagonal_bin_correction:
        last_idx = hic_mat.shape[0] - 1
        nonzero_diag = hic_mat.nonzero()[0][
            hic_mat.nonzero()[0] == hic_mat.nonzero()[1]
        ]
        nonzero_diag = list(
            set(nonzero_diag) - set(np.array([last_idx])) - set(np.array([0]))
        )

        for ii in nonzero_diag:
            hic_mat[ii, ii] = (
                max(hic_mat[ii, ii - 1], hic_mat[ii, ii + 1])
                * tss_hic_contribution
                / 100
            )

        if hic_mat[0, 0] != 0:
            hic_mat[0, 0] = hic_mat[0, 1] * tss_hic_contribution / 100

        if hic_mat[last_idx, last_idx] != 0:
            hic_mat[last_idx, last_idx] = (
                hic_mat[last_idx, last_idx - 1] * tss_hic_contribution / 100
            )

    # Any entries with low KR norm entries get set to NaN. These will be interpolated below
    hic_mat = apply_kr_threshold(hic_mat, hic_norm_file, kr_cutoff)

    # Remove lower triangle
    if not hic_is_vc:
        hic_mat = ssp.triu(hic_mat)
    else:
        hic_mat = process_vc(hic_mat)

    # Turn into dataframe
    hic_mat = hic_mat.tocoo(copy=False)
    hic_df = pd.DataFrame(
        {"bin1": hic_mat.row, "bin2": hic_mat.col, "hic_contact": hic_mat.data}
    )

    # Prune to window
    hic_df = hic_df.loc[
        np.logical_and(
            abs(hic_df["bin1"] - hic_df["bin2"]) <= window / resolution,
            abs(hic_df["bin1"] - hic_df["bin2"]) >= min_window / resolution,
        )
    ]
    print(
        "HiC has {} rows after windowing between {} and {}".format(
            hic_df.shape[0], min_window, window
        )
    )

    hic_df["juicebox_contact_values"] = hic_df["hic_contact"]
    # Fill NaN
    if interpolate_nan:
        interpolate_nan_in_hic(hic_df, gamma, scale, resolution)

    print("process.hic: Elapsed time: {}".format(time.time() - t))

    return hic_df


def interpolate_nan_in_hic(hic_df, gamma, scale, resolution) -> None:
    # NaN in the KR normalized matrix are not zeros. They are entries where the KR algorithm did not converge (or low KR norm)
    # So need to fill these. Use powerlaw.
    # Not ideal obviously but the scipy interpolation algos are either very slow or don't work since the nan structure implies that not all nans are interpolated
    nan_loc = np.isnan(hic_df["hic_contact"])
    hic_df.loc[nan_loc, "hic_contact"] = get_powerlaw_at_distance(
        abs(hic_df.loc[nan_loc, "bin1"] - hic_df.loc[nan_loc, "bin2"]) * resolution,
        gamma,
        scale,
        min_distance=resolution,
    )


def apply_kr_threshold(hic_mat, hic_norm_file, kr_cutoff):
    # Convert all entries in the hic matrix corresponding to low kr norm entries to NaN
    # Note that in scipy sparse matrix multiplication 0*nan = 0
    # So this doesn't convert 0's to nan only nonzero to nan
    norms = np.loadtxt(hic_norm_file)
    norms[norms < kr_cutoff] = np.nan
    norms[norms >= kr_cutoff] = 1
    norm_mat = ssp.dia_matrix((1.0 / norms, [0]), (len(norms), len(norms)))

    return norm_mat * hic_mat * norm_mat


def hic_to_sparse(filename, norm_file, resolution, hic_is_doubly_stochastic=False):
    t = time.time()
    HiC = pd.read_table(
        filename,
        names=["bin1", "bin2", "hic_contact"],
        header=None,
        engine="c",
        memory_map=True,
    )

    # verify our assumptions
    assert np.all(HiC.bin1 <= HiC.bin2)

    # Need load norms here to know the dimensions of the hic matrix
    norms = pd.read_csv(norm_file, header=None)
    hic_size = norms.shape[0]

    # convert to sparse matrix in CSR (compressed sparse row) format, chopping
    # down to HiC bin size.  note that conversion to scipy sparse matrices
    # accumulates repeated indices, so this will do the right thing.
    row = np.floor(HiC.bin1.values / resolution).astype(int)
    col = np.floor(HiC.bin2.values / resolution).astype(int)
    dat = HiC.hic_contact.values

    # JN: Need both triangles in order to compute row/column sums to make double stochastic.
    # If juicebox is upgraded to return DS matrices, then can remove one triangle
    # TO DO: Remove one triangle when juicebox is updated.
    # we want a symmetric matrix.  Easiest to do that during creation, but have to be careful of diagonal
    if not hic_is_doubly_stochastic:
        mask = row != col  # off-diagonal
        row2 = col[mask]  # note the row/col swap
        col2 = row[mask]
        dat2 = dat[mask]

        # concat and create
        row = np.hstack((row, row2))
        col = np.hstack((col, col2))
        dat = np.hstack((dat, dat2))

    print("hic.to.sparse: Elapsed time: {}".format(time.time() - t))

    return ssp.csr_matrix((dat, (row, col)), (hic_size, hic_size))


def get_powerlaw_at_distance(distances, gamma, scale, min_distance=5000):
    assert gamma > 0
    assert scale > 0

    # The powerlaw is computed for distances > 5kb. We don't know what the contact freq looks like at < 5kb.
    # So just assume that everything at < 5kb is equal to 5kb.
    # TO DO: get more accurate powerlaw at < 5kb
    distances = np.clip(distances, min_distance, np.Inf)
    log_dists = np.log(distances + 1)

    powerlaw_contact = np.exp(scale + -1 * gamma * log_dists)
    return powerlaw_contact


def process_vc(hic):
    # For a vc normalized matrix, need to make rows sum to 1.
    # Assume rows correspond to genes and cols to enhancers

    row_sums = hic.sum(axis=0)
    row_sums[row_sums == 0] = 1
    norm_mat = ssp.dia_matrix(
        (1.0 / row_sums, [0]), (row_sums.shape[1], row_sums.shape[1])
    )

    # left multiply to operate on rows
    hic = norm_mat * hic

    return hic


def make_predictions(
    chromosome, enhancers, genes, window, tss_slop, hic_file, hic_type, hic_resolution, 
    hic_gamma, hic_scale, hic_gamma_reference, tss_hic_contribution, scale_hic_using_powerlaw,
    chrom_sizes_map
):
    pred = make_pred_table(chromosome, enhancers, genes, window, chrom_sizes_map)
    pred = annotate_predictions(pred, tss_slop)
    pred = add_powerlaw_to_predictions(pred, hic_gamma_reference, hic_gamma, hic_scale)
    # if Hi-C file is not provided, only powerlaw model will be computed
    if hic_file:
        if hic_type == "hic":
            pred = add_hic_from_hic_file(
                pred, hic_file, chromosome, hic_resolution, window
            )
        else:
            pred = add_hic_from_directory(
                chromosome,
                enhancers,
                genes,
                pred,
                hic_file,
                hic_type,
                hic_resolution,
                tss_hic_contribution,
                window,
                hic_gamma,
                hic_scale,
                chrom_sizes_map,
            )

        # Remove all NaN values so we have valid scores
        pred.fillna(value={"hic_contact": 0}, inplace=True)
        # Add powerlaw scaling
        pred = scale_hic_with_powerlaw(pred, scale_hic_using_powerlaw)
        # Add pseudocount
        pred = add_hic_pseudocount(pred)
        print("HiC Complete")

        pred = compute_score(
            pred,
            [pred["activity_base"], pred["hic_contact_pl_scaled_adj"]],
            "ABC",
            adjust_self_promoters=True,
        )
    pred = compute_score(
        pred,
        [pred["activity_base"], pred["powerlaw_contact"]],
        "powerlaw",
        adjust_self_promoters=True,
    )
    return pred


def make_pred_table(chromosome, enh, genes, window, chrom_sizes_map: Dict[str, int]):
    print("Making putative predictions table...")
    t = time.time()
    enh["enh_midpoint"] = (enh["start"] + enh["end"]) / 2
    enh["enh_idx"] = enh.index
    genes["gene_idx"] = genes.index
    enh_pr = df_to_pyranges(enh)
    genes_pr = df_to_pyranges(
        genes,
        start_col="TargetGeneTSS",
        end_col="TargetGeneTSS",
        start_slop=window,
        end_slop=window,
        chrom_sizes_map=chrom_sizes_map,
    )

    pred = enh_pr.join(genes_pr).df.drop(
        ["Start_b", "End_b", "chr_b", "Chromosome", "Start", "End"], axis=1
    )
    pred["distance"] = abs(pred["enh_midpoint"] - pred["TargetGeneTSS"])
    pred = pred.loc[pred["distance"] < window, :]  # for backwards compatability

    print(
        "Done. There are {} putative enhancers for chromosome {}".format(
            pred.shape[0], chromosome
        )
    )
    print("Elapsed time: {}".format(time.time() - t))

    return pred


def create_df_from_records(records, hic_resolution):
    bin_data = [[r.binX, r.binY, r.counts] for r in records]
    df = pd.DataFrame(bin_data, columns=["binX", "binY", "counts"])
    df["binX"] = np.floor(df["binX"] / hic_resolution).astype(int)
    df["binY"] = np.floor(df["binY"] / hic_resolution).astype(int)
    return df


def get_chrom_format(hic: hicstraw.HiCFile, chromosome):
    """
    hic files can have 'chr1' or just '1' as the chromosome name
    we need to make sure we're using the format consistent with
    the hic file
    """
    hic_chrom_names = [chrom.name for chrom in hic.getChromosomes()]
    if hic_chrom_names[1].startswith("chr"):  # assume index 1 should be chr1
        return chromosome
    else:
        return chromosome[3:]


def add_hic_from_hic_file(pred, hic_file, chromosome, hic_resolution, window):
    pred["enh_bin"] = np.floor(pred["enh_midpoint"] / hic_resolution).astype(int)
    pred["tss_bin"] = np.floor(pred["TargetGeneTSS"] / hic_resolution).astype(int)
    pred["binX"] = np.min(pred[["enh_bin", "tss_bin"]], axis=1)
    pred["binY"] = np.max(pred[["enh_bin", "tss_bin"]], axis=1)
    hic = hicstraw.HiCFile(hic_file)
    chromosome = get_chrom_format(hic, chromosome)
    matrix_object = hic.getMatrixZoomData(
        chromosome, chromosome, "observed", "SCALE", "BP", hic_resolution
    )
    start_loci = pred["end"].min()
    end_loci = pred["end"].max()
    step_size = 10000 * hic_resolution  # ~10k bins at a time
    for i in range(start_loci, end_loci, step_size):
        start = i
        end = start + step_size
        records = matrix_object.getRecords(start, end, start - window, end + window)
        df = create_df_from_records(records, hic_resolution)
        pred = pred.merge(df, how="left", on=["binX", "binY"], suffixes=(None, "_"))
        if "counts_" in pred:
            pred["counts"] = np.max(pred[["counts", "counts_"]], axis=1)
            pred.drop("counts_", inplace=True, axis=1)

    pred.drop(
        [
            "binX",
            "binY",
            "enh_idx",
            "gene_idx",
            "enh_midpoint",
            "tss_bin",
            "enh_bin",
        ],
        inplace=True,
        axis=1,
        errors="ignore",
    )

    return pred.rename(columns={"counts": "hic_contact"})


def add_hic_from_directory(
    chromosome,
    enh,
    genes,
    pred,
    hic_dir,
    hic_type,
    hic_resolution,
    tss_hic_contribution,
    window,
    hic_gamma,
    hic_scale,
    chrom_sizes_map,
):
    hic_file, hic_norm_file, hic_is_vc = get_hic_file(
        chromosome, hic_dir, hic_type=hic_type
    )
    print("Begin HiC")
    # Add hic to pred table
    # At this point we have a table where each row is an enhancer/gene pair.
    # We need to add the corresponding HiC matrix entry.
    # If the HiC is provided in juicebox format (ie constant resolution), then we can just merge using the indices
    # But more generally we do not want to assume constant resolution. In this case hic should be provided in bedpe format

    t = time.time()
    if hic_type == "bedpe":
        HiC = load_hic_bedpe(hic_file)
        # Use pyranges to compute overlaps between enhancers/genes and hic bedpe table
        # Consider each range of the hic matrix separately - and merge each range into both enhancers and genes.
        # Then remerge on hic index

        HiC["hic_idx"] = HiC.index
        hic1 = df_to_pyranges(HiC, start_col="x1", end_col="x2", chr_col="chr1")
        hic2 = df_to_pyranges(HiC, start_col="y1", end_col="y2", chr_col="chr2")

        # Overlap in one direction
        enh_hic1 = (
            df_to_pyranges(
                enh,
                start_col="enh_midpoint",
                end_col="enh_midpoint",
                end_slop=1,
                chrom_sizes_map=chrom_sizes_map,
            )
            .join(hic1)
            .df
        )
        genes_hic2 = (
            df_to_pyranges(
                genes,
                start_col="TargetGeneTSS",
                end_col="TargetGeneTSS",
                end_slop=1,
                chrom_sizes_map=chrom_sizes_map,
            )
            .join(hic2)
            .df
        )
        ovl12 = enh_hic1[["enh_idx", "hic_idx", "hic_contact"]].merge(
            genes_hic2[["gene_idx", "hic_idx"]], on="hic_idx"
        )

        # Overlap in the other direction
        enh_hic2 = (
            df_to_pyranges(
                enh,
                start_col="enh_midpoint",
                end_col="enh_midpoint",
                end_slop=1,
                chrom_sizes_map=chrom_sizes_map,
            )
            .join(hic2)
            .df
        )
        genes_hic1 = (
            df_to_pyranges(
                genes,
                start_col="TargetGeneTSS",
                end_col="TargetGeneTSS",
                end_slop=1,
                chrom_sizes_map=chrom_sizes_map,
            )
            .join(hic1)
            .df
        )
        ovl21 = enh_hic2[["enh_idx", "hic_idx", "hic_contact"]].merge(
            genes_hic1[["gene_idx", "hic_idx"]], on=["hic_idx"]
        )

        # Concatenate both directions and merge into preditions
        ovl = pd.concat([ovl12, ovl21]).drop_duplicates()
        pred = pred.merge(ovl, on=["enh_idx", "gene_idx"], how="left")
    elif hic_type == "juicebox" or hic_type == "avg":
        if hic_type == "juicebox":
            HiC = load_hic_juicebox(
                hic_file=hic_file,
                hic_norm_file=hic_norm_file,
                hic_is_vc=hic_is_vc,
                hic_resolution=hic_resolution,
                tss_hic_contribution=tss_hic_contribution,
                window=window,
                min_window=0,
                gamma=hic_gamma,
                scale=hic_scale,
            )
        else:
            HiC = load_hic_avg(hic_file, hic_resolution)

        # Merge directly using indices
        # Could also do this by indexing into the sparse matrix (instead of merge) but this seems to be slower
        # Index into sparse matrix
        # pred['hic_contact'] = [HiC[i,j] for (i,j) in pred[['enh_bin','tss_bin']].values.tolist()]

        pred["enh_bin"] = np.floor(pred["enh_midpoint"] / hic_resolution).astype(
            int
        )
        pred["tss_bin"] = np.floor(pred["TargetGeneTSS"] / hic_resolution).astype(
            int
        )
        if not hic_is_vc:
            # in this case the matrix is upper triangular.
            #
            pred["bin1"] = np.amin(pred[["enh_bin", "tss_bin"]], axis=1)
            pred["bin2"] = np.amax(pred[["enh_bin", "tss_bin"]], axis=1)
            pred = pred.merge(HiC, how="left", on=["bin1", "bin2"])
        else:
            # The matrix is not triangular, its full
            # For VC assume genes correspond to rows and columns to enhancers
            pred = pred.merge(
                HiC,
                how="left",
                left_on=["tss_bin", "enh_bin"],
                right_on=["bin1", "bin2"],
            )
        # QC juicebox HiC
        pred = qc_hic(pred)

    pred.drop(
        [
            "x1",
            "x2",
            "y1",
            "y2",
            "bin1",
            "bin2",
            "enh_idx",
            "gene_idx",
            "hic_idx",
            "enh_midpoint",
            "tss_bin",
            "enh_bin",
        ],
        inplace=True,
        axis=1,
        errors="ignore",
    )

    print("HiC added to predictions table. Elapsed time: {}".format(time.time() - t))
    return pred


def scale_hic_with_powerlaw(pred, scale_hic_using_powerlaw):
    # Scale hic values to reference powerlaw

    if not scale_hic_using_powerlaw:
        #        values = pred.loc[pred['hic_contact']==0].index.astype('int')
        #        pred.loc[values, 'hic_contact'] = pred.loc[values, 'powerlaw_contact']
        pred["hic_contact_pl_scaled"] = pred["hic_contact"]
    else:
        pred["hic_contact_pl_scaled"] = pred["hic_contact"] * (
            pred["powerlaw_contact_reference"] / pred["powerlaw_contact"]
        )

    return pred


def add_powerlaw_to_predictions(pred, hic_gamma_reference, hic_gamma, hic_scale):
    pred["powerlaw_contact"] = get_powerlaw_at_distance(
        pred["distance"].values, hic_gamma, hic_scale
    )

    # 4.80 and 11.63 come from a linear regression of scale on gamma across 20
    # hic cell types at 5kb resolution. Do the params change across resolutions?
    hic_scale_reference = -4.80 + 11.63 * hic_gamma_reference
    pred["powerlaw_contact_reference"] = get_powerlaw_at_distance(
        pred["distance"].values, hic_gamma_reference, hic_scale_reference
    )

    return pred


def add_hic_pseudocount(pred):
    # Add a pseudocount based on the powerlaw expected count at a given distance

    pseudocount = pred[["powerlaw_contact", "powerlaw_contact_reference"]].min(axis=1)
    pred["hic_pseudocount"] = pseudocount
    pred["hic_contact_pl_scaled_adj"] = pred["hic_contact_pl_scaled"] + pseudocount

    return pred


def qc_hic(pred, threshold=0.01):
    # Genes with insufficient hic coverage should get nan'd

    summ = (
        pred.loc[pred["isSelfPromoter"], :]
        .groupby(["TargetGene"])
        .agg({"hic_contact": "sum"})
    )
    bad_genes = summ.loc[summ["hic_contact"] < threshold, :].index

    pred.loc[pred["TargetGene"].isin(bad_genes), "hic_contact"] = np.nan

    return pred


def compute_score(enhancers, product_terms, prefix, adjust_self_promoters=True):
    scores = np.column_stack(product_terms).prod(axis=1)

    enhancers[prefix + ".Score.Numerator"] = scores
    enhancers[prefix + ".Score"] = enhancers[
        prefix + ".Score.Numerator"
    ] / enhancers.groupby(["TargetGene", "TargetGeneTSS"])[
        prefix + ".Score.Numerator"
    ].transform(
        "sum"
    )

    # Self promoters by definition regulate the gene, so we
    # want to make sure they have a high score
    if adjust_self_promoters:
        self_promoters = enhancers[enhancers["isSelfPromoter"]]
        enhancers.loc[self_promoters.index, prefix + ".Score"] = 1

    return enhancers


def annotate_predictions(pred, tss_slop=500):
    # TO DO: Add is self genic
    pred["isSelfPromoter"] = np.logical_and.reduce(
        (
            pred["class"] == "promoter",
            pred.start - tss_slop < pred.TargetGeneTSS,
            pred.end + tss_slop > pred.TargetGeneTSS,
        )
    )

    return pred


def make_gene_prediction_stats(pred, score_column, threshold, output_file):
    summ1 = pred.groupby(["chr", "TargetGene", "TargetGeneTSS"]).agg(
        {
            "TargetGeneIsExpressed": lambda x: set(x).pop(),
            score_column: lambda x: all(np.isnan(x)),
            "name": "count",
        }
    )
    summ1.columns = ["geneIsExpressed", "geneFailed", "nEnhancersConsidered"]

    summ2 = (
        pred.loc[pred["class"] != "promoter", :]
        .groupby(["chr", "TargetGene", "TargetGeneTSS"])
        .agg({score_column: lambda x: sum(x > threshold)})
    )
    summ2.columns = ["nDistalEnhancersPredicted"]
    summ1 = summ1.merge(summ2, left_index=True, right_index=True)

    summ1.to_csv(output_file, sep="\t", index=True)



def test_variant_overlap(outdir, score_column, all_putative, chrom_sizes):
    variant_overlap_file = os.path.join(
        outdir,
        "EnhancerPredictionsAllPutative.ForVariantOverlap.shrunk150bp.tsv.gz",
    )
    # generate predictions for variant overlap
    score_t = all_putative[score_column] > 0.015
    not_promoter = all_putative["class"] != "promoter"
    is_promoter = all_putative["class"] == "promoter"
    score_one = all_putative[score_column] > 0.1
    all_putative[(score_t & not_promoter) | (is_promoter & score_one)]
    variant_overlap = all_putative[(score_t & not_promoter) | (is_promoter & score_one)]
    # remove nan predictions
    variant_overlap_pred = variant_overlap.dropna(subset=[score_column])
    variant_overlap = variant_overlap_pred.loc[
        variant_overlap_pred["distance"] <= 2000000
    ]
    variant_overlap.to_csv(
        variant_overlap_file + ".tmp",
        sep="\t",
        index=False,
        header=True,
        compression="gzip",
        float_format="%.6f",
    )
    # shrink regions
    os.system(
        "zcat {}.tmp 2>/dev/null | head -1 | gzip > {}".format(
            variant_overlap_file, variant_overlap_file
        )
    )
    os.system(
        "zcat {}.tmp | sed 1d | bedtools slop -b -150 -g {} | gzip >> {}".format(
            variant_overlap_file, chrom_sizes, variant_overlap_file
        )
    )


def determine_expressed_genes(genes, expression_cutoff, activity_quantile_cutoff):
    # Evaluate whether a gene should be considered 'expressed' so that it runs through the model

    # A gene is runnable if:
    # It is expressed OR (there is no expression AND its promoter has high activity)

    genes["isExpressed"] = np.logical_or(
        genes.Expression >= expression_cutoff,
        np.logical_and(
            np.isnan(genes.Expression),
            genes.PromoterActivityQuantile >= activity_quantile_cutoff,
        ),
    )

    return genes


def main(enhancers, genes, chrom_sizes, hic_file, outdir):
    shutil.rmtree(outdir, ignore_errors=True)
    os.makedirs(outdir)

    score_column = "ABC.Score"

    accessibility_feature = "ATAC"

    cellType = None

    hic_resolution = 5000
    hic_type = "hic"
    scale_hic_using_powerlaw = True

    hic_gamma = 1.024238616787792
    hic_scale = 5.9594510043736655
    hic_gamma_reference = 0.87

    expression_cutoff = 1
    promoter_activity_quantile_cutoff = 0.3

    window = 5000000
    tss_hic_contribution = 100

    tss_slop = 500

    chromosomes = "all"
    include_chrY = True

    genes = pd.read_csv(genes, sep="\t")
    genes = determine_expressed_genes(
        genes, expression_cutoff, promoter_activity_quantile_cutoff
    )

    genes_columns_to_subset = [
        "chr",
        "symbol",
        "tss",
        "Expression",
        "PromoterActivityQuantile",
        "isExpressed",
        "Ensembl_ID",
    ]
    genes_column_names = [
        "chr",
        "TargetGene",
        "TargetGeneTSS",
        "TargetGeneExpression",
        "TargetGenePromoterActivityQuantile",
        "TargetGeneIsExpressed",
        "TargetGeneEnsembl_ID",
    ]
    print("reading enhancers")
    enhancers_full = pd.read_csv(enhancers, sep="\t")
    enhancers_column_names = ["chr", "start", "end", "name", "class", "activity_base"]
    if accessibility_feature not in {"ATAC", "DHS"}:
        raise ValueError("The feature has to be either ATAC or DHS!")
    normalized_activity_col = f"normalized_{accessibility_feature.lower()}"

    genes_subset_columns = genes_columns_to_subset + [
        f"{accessibility_feature}.RPKM.quantile.TSS1Kb"
    ]
    genes = genes.loc[:, genes_subset_columns]
    genes.columns = genes_column_names + [normalized_activity_col]

    enh_subset_columns = enhancers_column_names + [normalized_activity_col]
    enhancers = enhancers_full.loc[:, enh_subset_columns]

    enhancers["activity_base_squared"] = enhancers["activity_base"] ** 2
    # Initialize Prediction files
    all_pred_file_expressed = os.path.join(
        outdir, "EnhancerPredictionsAllPutative.tsv.gz"
    )
    all_pred_file_nonexpressed = os.path.join(
        outdir, "EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz"
    )
    all_putative_list = []

    # Make predictions
    if chromosomes == "all":
        chromosomes = set(genes["chr"]).intersection(set(enhancers["chr"]))
        if not include_chrY:
            chromosomes.discard("chrY")
        chromosomes = sorted(chromosomes)
    else:
        chromosomes = chromosomes.split(",")

    chrom_sizes_map = pd.read_csv(
        chrom_sizes, sep="\t", header=None, index_col=0
    ).to_dict()[1]

    for chromosome in chromosomes:
        print("Making predictions for chromosome: {}".format(chromosome))
        t = time.time()
        this_enh = enhancers.loc[enhancers["chr"] == chromosome, :].copy()
        this_genes = genes.loc[genes["chr"] == chromosome, :].copy()

        this_chr = make_predictions(
            chromosome,
            this_enh,
            this_genes,
            window, 
            tss_slop, 
            hic_file, 
            hic_type, 
            hic_resolution,
            hic_gamma,
            hic_scale,
            hic_gamma_reference, 
            tss_hic_contribution, 
            scale_hic_using_powerlaw,
            chrom_sizes_map,
        )
        all_putative_list.append(this_chr)

        print(
            "Completed chromosome: {}. Elapsed time: {} \n".format(
                chromosome, time.time() - t
            )
        )

    # Subset predictions
    all_putative = pd.concat(all_putative_list)
    all_putative["CellType"] = cellType
    if hic_file:
        all_putative["hic_contact_squared"] = all_putative["hic_contact"] ** 2

    all_putative.loc[all_putative.TargetGeneIsExpressed, :].to_csv(
        all_pred_file_expressed,
        sep="\t",
        index=False,
        header=True,
        compression="gzip",
        float_format="%.6f",
        na_rep="NaN",
    )
    all_putative.loc[~all_putative.TargetGeneIsExpressed, :].to_csv(
        all_pred_file_nonexpressed,
        sep="\t",
        index=False,
        header=True,
        compression="gzip",
        float_format="%.6f",
        na_rep="NaN",
    )

    test_variant_overlap(outdir, score_column, all_putative, chrom_sizes)


enhancers = snakemake.input.enhancers
genes = snakemake.input.genes
chrom_sizes = snakemake.input.chrom_sizes
hic_file = snakemake.input.hic

outdir = os.path.dirname(snakemake.output.predictions)

main(enhancers, genes, chrom_sizes, hic_file, outdir)

