########################
## Benjamin Haibe-Kains
## May 09, 2011
########################

rm(list=ls())

library(Hmisc)

saveres <- "saveres"
if(!file.exists(saveres)) { system(sprintf("mkdir %s", saveres)) }

rdir <- "all_ranking"
rfile <- "all2.13091"

rfn <- dir(rdir)
rfn <- rfn[grep(pattern=rfile, x=rfn)]
rfnn <- sapply(strsplit(rfn, "[.]"), function(x) { return(paste("mFS",x[[3]], x[[4]], sep="")) })

## load our affy datasets
load("../data.RData")
annots <- datas.m$annots[[1]]

load(sprintf("%s/%s", rdir, rfn[1]))
## mFS is a list with one item containing the indices of all the genes that passed the filtering
ddn <- sort(annots[[1]][mFS$msel])

## estimation of the concordance index (actually, just the sign) for all genes
mydd <- c("UPP", "STK", "VDX", "UNT", "MAINZ", "TRANSBIG")
mystat <- matrix(NA, ncol=length(datas.m$data), nrow=length(ddn), dimnames=list(ddn, names(datas.m$data)))
for(i in 1:length(mydd)) {
	mystat[ , i] <- apply(datas.m$data[[mydd[i]]][ , ddn, drop=FALSE], 2, rcorr.cens, S=datas.m$demo[[mydd[i]]][ , "surv.bin"], outx=TRUE)["Dxy", ]
}
mystat <- sign(apply(mystat, 1, mean, na.rm=TRUE))
rm(list=c("datas.m"))
gc()

## load the "causal" rankings and compare them with the traditional ranking
pdf(sprintf("%s/plot_causal_vs_traditional_ranking5.pdf", saveres))

## traditional ranking
## affy probesets
mymFS <- rownames(annots)[mFS$usel]
write.table(cbind(mymFS, (length(mymFS):1) * mystat[mymFS]), file=sprintf("%s/mFS00_mimo.rnk", saveres), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
## Entrez gene ids
mymFS2 <- annots[mFS$usel, "EntrezGene.ID"]
write.table(cbind(mymFS2, (length(mymFS2):1) * mystat[mymFS]), file=sprintf("%s/mFS00_entrez_mimo.rnk", saveres), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
rm(mFS)

for(i in 1:length(rfn)) {
	load(sprintf("%s/%s", rdir, rfn[i]))
	## affy probesets
	mymFS <- rownames(annots)[mFS$msel]
	write.table(cbind(mymFS, (length(mymFS):1) * mystat[mymFS]), file=sprintf("%s/%s_mimo.rnk", saveres, rfnn[i]), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
	## Entrez gene ids
	mymFS2 <- annots[mFS$msel, "EntrezGene.ID"]
	write.table(cbind(mymFS2, (length(mymFS2):1) * mystat[mymFS]), file=sprintf("%s/%s_entrez_mimo.rnk", saveres, rfnn[i]), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
	tt <- table(mFS$usel != mFS$msel)
	plot(mFS$usel, mFS$msel, xlab="Traditional ranking (usel)", ylab="Causal ranking (msel)", main=sprintf("LAMBDA = %s.%s", substr(rfnn[i], nchar(rfnn[i])-1, nchar(rfnn[i])-1), substr(rfnn[i], nchar(rfnn[i]), nchar(rfnn[i]))), sub=sprintf("%.3g%% of the ranking is different", (tt[2] / sum(tt) * 100)))
	rm(mFS)
}

dev.off()

