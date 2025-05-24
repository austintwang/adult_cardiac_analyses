
# read count and offset files (text)
ytxt=snakemake@input$ytxt
ktxt=snakemake@input$ktxt
xtxt=NA
# read tables
y=read.table(ytxt,as.is=T)
k=read.table(ktxt,as.is=T)
# if(!is.na(xtxt)){x=read.table(xtxt,as.is=T)}
# output binary file names
ybin=snakemake@output$ybin
kbin=snakemake@output$kbin
if(!is.na(xtxt)){xbin=gsub("txt", "bin", xtxt)}
# open files
fybin=file(ybin,"wb")
fkbin=file(kbin,"wb")
if(!is.na(xtxt)){fxbin=file(xbin,"wb")}
# write tables as binary
writeBin(as.double(c(t(y[,-1]))), fybin)
writeBin(as.double(c(t(k[,-1]))), fkbin)
if(!is.na(xtxt)){writeBin(as.double(c(as.matrix(x))), fxbin)}
# close files
close(fybin)
close(fkbin)
if(!is.na(xtxt)){close(fxbin)}