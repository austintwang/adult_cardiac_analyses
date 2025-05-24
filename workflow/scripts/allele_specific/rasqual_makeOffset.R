
Quantile <-
function(x,k=20){x=rank(x,ties="random"); z=rep(0,length(x));for(i in 1:k){z=z+as.numeric(x<=quantile(x,i/k,na.rm=T))};k-z}
solveDuplication <- 
function(x){
	x2=x
	for(i in as.numeric(names(table(x)[table(x)>1]))){
		ids=seq(length(x))[x==i];
		if(ids[1]>1){
			a = (x[ids[1]-1]+x[ids[1]])/2
		}else{
			a = x[1]
		}
		if(ids[length(ids)]<length(x)){
			b = (x[ids[length(ids)]]+x[ids[length(ids)]+1])/2
		}else{
			b = x[length(x)]
		}
		w = (b-a)/(length(ids)+1)
		x2[ids] = a+w*(1:length(ids))
	}
	x2
}
gcCor <-
function(Y,gcvec,PLOT=F){
	bin=Quantile(gcvec,200);
	x=sort(unlist(lapply(split(gcvec,bin),mean)))
	x=solveDuplication(x)
	S=apply(Y,2,function(y){unlist(lapply(split(y,bin),sum))[as.character(0:199)]});
	Fs=log(t(t(S)/apply(S,2,sum))/apply(S,1,sum)*sum(S));
	Gs=apply(Fs,2,function(y){smooth.spline(x,y,spar=1)$y}); 
	if(PLOT){
		par(mfcol=c(5,5),mar=c(2,2,2,2)); 
		for(i in 1:ncol(Y)){
			plot(Fs[,i])
			lines(Gs[,i],col=2)
		}
		matplot(x,Gs,type="l",col=2,lty=1)
	}
	exp(Gs[bin+1,])
}


# read count and gc content vector (text)
ytxt=snakemake@input$ytxt
gcc=NA

# read tables
Y=read.table(ytxt,as.is=T)
fid=Y[[1]]
Y=as.matrix(Y[,-1])
K=array(1, dim(Y))

# GC correction
if(!is.na(gcc)){
	gcc=scan(gcc)
	K=gcCor(Y,gcc)
}

# library size
sf=apply(Y,2,sum)
sf=sf/mean(sf)

# write K as binary
output=snakemake@output$ktxt
write.table(data.frame(fid, t(t(K)*sf)+1e-4), col=F, row=F, sep="\t", file=output, quote=F)