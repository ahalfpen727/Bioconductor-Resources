########################
## Benjamin Haibe-Kains
## May 09, 2011
########################

## directories
saveres1 <- "saveres"
saveres <- "saveres2"
if(!file.exists(saveres)) { system(sprintf("mkdir %s", saveres)) }

## GSEA parameters
exe.path <- "/Applications/gsea2-2.07.jar"
#gmt.path <- "~/MyWork/GSEA/c5.bp.v3.0.symbols.gmt"
gmt.path <- "c5.bp.v3.0.entrez.gmt"
#gsea.collapse <- "true -mode Max_probe -norm meandiv"
gsea.collapse <- "false"
nperm <- 1000
gsea.seed <- 54321
gsea.out <- saveres

#resn.all <- sprintf("mFS%s_mimo", c("00", "02", "05", "07", "10"))
lambdan <- c("00", "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "20")
resn.all <- sprintf("mFS%s_entrez_mimo", lambdan)

## run preranked GSEA analyses for each lambda
gsea.res <- NULL
for(i in 1:length(resn.all)) {
	resn <- resn.all[i]
	## build GSEA command
	rnk.path <- sprintf("%s/%s.rnk", saveres1, resn)
	gsea.report <- resn
	rest <- dir(saveres)
	rest <- rest[grep(pattern=resn, x=rest)[1]]
	if(length(rest) == 0 || is.na(rest)) {
		#gsea.cmd <- sprintf("java -Xmx512m -cp %s xtools.gsea.GseaPreranked -gmx %s -collapse %s -nperm %i -rnk %s -scoring_scheme weighted -rpt_label %s -chip %s -include_only_symbols true -make_sets false -plot_top_x 20 -rnd_seed %i -set_max 500 -set_min 15 -zip_report false -out %s -gui false", exe.path, gmt.path, gsea.collapse, nperm, rnk.path, gsea.report, gsea.chip, gsea.seed, gsea.out)
		gsea.cmd <- sprintf("java -Xmx512m -cp %s xtools.gsea.GseaPreranked -gmx %s -collapse %s -nperm %i -rnk %s -scoring_scheme weighted -rpt_label %s -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed %i -set_max 500 -set_min 15 -zip_report false -out %s -gui false", exe.path, gmt.path, gsea.collapse, nperm, rnk.path, gsea.report, gsea.seed, gsea.out)
		system(gsea.cmd)
		cat("\n-------------------------------------\n\n")
		## read results
		rest <- dir(saveres)
		rest <- rest[grep(pattern=resn, x=rest)[1]]
	}
	restn <- sapply(strsplit(rest, "[.]"), function(x) { return(x[length(x)]) })
	tt <- rbind(read.table(sprintf("%s/%s/gsea_report_for_na_pos_%s.xls", saveres, rest, restn), stringsAsFactors=FALSE, sep="\t", header=TRUE), read.table(sprintf("%s/%s/gsea_report_for_na_neg_%s.xls", saveres, rest, restn), stringsAsFactors=FALSE, sep="\t", header=TRUE))
	gsea.res <- c(gsea.res, list(tt))
}
names(gsea.res) <- resn.all

## extract the most relevant enriched gene sets
topn2 <- 3
topn <- 50
if(topn > nrow(gsea.res[[1]])) { topn <- nrow(gsea.res[[1]]) }
maxp <- 2
nn <- NULL
for(i in 1:length(gsea.res)) {
	tt <- gsea.res[[i]][ , "NAME"][gsea.res[[i]][ , "FWER.p.val"] < maxp]
	tt <- tt[order(abs(gsea.res[[i]][match(tt, gsea.res[[i]][ , "NAME"]), "NES"]), decreasing=TRUE)][1:topn]
	nn <- c(nn, tt)
}
nn <- unique(nn[!is.na(nn)])

matres <- matpres <- matrix(NA, ncol=length(gsea.res), nrow=length(nn), dimnames=list(nn, names(gsea.res)))
for(i in 1:length(gsea.res)) {
	matres[ , i] <- gsea.res[[i]][match(nn, gsea.res[[i]][ ,"NAME"]), "NES"]
	matpres[ , i] <- gsea.res[[i]][match(nn, gsea.res[[i]][ ,"NAME"]), "FWER.p.val"]
}

## resort based on the max NES
tt <- tt[order(apply(matres[tt, , drop=FALSE], 1, max, na.rm=TRUE), decreasing=TRUE)]
write.csv(matres[tt, , drop=FALSE], file=sprintf("%s/gsea_res_all.csv", saveres))

rrlm <- apply(matres, 1, function(x, y) { return(lm(abs(x) ~ y)) }, y=as.numeric(lambdan)/10)
rrcoef <- sapply(rrlm, function(x) { return(x$coefficients[2]) })

## top more causal
tt <- order(rrcoef, decreasing=TRUE)[1:topn2]
tt <- tt[rrcoef[tt] > 0]
tt <- tt[order(matres[tt, 1], decreasing=TRUE)]

matres2 <- rep(NA, length(tt) * ncol(matres))
for(i in 1:ncol(matres)) {
	matres2[seq(i, length(tt) * ncol(matres), by=ncol(matres))] <- matres[tt, i, drop=FALSE]
}
mysp <- rep(0.1, length(matres2))
mysp[seq(ncol(matres), length(mysp) - 1, by=ncol(matres))] <- 1
mycol <- rep(rev(rainbow(ncol(matres), v=0.8, alpha=0.5)), length(tt))
mylab <- rep(NA, length(matres2))
mylab[seq(1, length(matres2), by=ncol(matres))] <- dimnames(matres[tt, ,drop=FALSE])[[1]]
names(matres2) <- mylab

## with labels
pdf(sprintf("%s/GSEA_res_morecausal_plot_labels.pdf", saveres), width=10, height=5)
par(las=1, mar=c(5, 26, 4, 2) + 0.1)
barplot(height=rev(matres2), space=rev(mysp), horiz=TRUE, col=mycol, xlab="Normalized Enrichment Score", border=NA, xlim=c(-2.5, 2.5))
dev.off()

## without lablels
pdf(sprintf("%s/GSEA_res_morecausal_plot.pdf", saveres), width=5, height=5)
barplot(height=rev(matres2), space=rev(mysp), horiz=TRUE, col=mycol, xlab="Normalized Enrichment Score", border=NA, names.arg="", xlim=c(-2.5, 2.5))
dev.off()

## top less causal
tt <- order(rrcoef, decreasing=FALSE)[1:topn2]
tt <- tt[rrcoef[tt] < 0]
tt <- tt[order(matres[tt, 1], decreasing=TRUE)]

matres2 <- rep(NA, length(tt) * ncol(matres))
for(i in 1:ncol(matres)) {
	matres2[seq(i, length(tt) * ncol(matres), by=ncol(matres))] <- matres[tt, i, drop=FALSE]
}
mysp <- rep(0.1, length(matres2))
mysp[seq(ncol(matres), length(mysp) - 1, by=ncol(matres))] <- 1
mycol <- rep(rev(rainbow(ncol(matres), v=0.8, alpha=0.5)), length(tt))
mylab <- rep(NA, length(matres2))
mylab[seq(1, length(matres2), by=ncol(matres))] <- dimnames(matres[tt, ,drop=FALSE])[[1]]
names(matres2) <- mylab

## with labels
pdf(sprintf("%s/GSEA_res_lesscausal_plot_labels.pdf", saveres), width=10, height=5)
par(las=1, mar=c(5, 26, 4, 2) + 0.1)
barplot(height=rev(matres2), space=rev(mysp), horiz=TRUE, col=mycol, xlab="Normalized Enrichment Score", border=NA, xlim=c(-2.5, 2.5))
dev.off()

## without lablels
pdf(sprintf("%s/GSEA_res_lesscausal_plot.pdf", saveres), width=5, height=5)
barplot(height=rev(matres2), space=rev(mysp), horiz=TRUE, col=mycol, xlab="Normalized Enrichment Score", border=NA, names.arg="", xlim=c(-2.5, 2.5))
dev.off()

## top same causal
tt <- order(abs(rrcoef), decreasing=FALSE)[1:topn2]
tt <- tt[order(matres[tt, 1], decreasing=TRUE)]

matres2 <- rep(NA, length(tt) * ncol(matres))
for(i in 1:ncol(matres)) {
	matres2[seq(i, length(tt) * ncol(matres), by=ncol(matres))] <- matres[tt, i, drop=FALSE]
}
mysp <- rep(0.1, length(matres2))
mysp[seq(ncol(matres), length(mysp) - 1, by=ncol(matres))] <- 1
mycol <- rep(rev(rainbow(ncol(matres), v=0.8, alpha=0.5)), length(tt))
mylab <- rep(NA, length(matres2))
mylab[seq(1, length(matres2), by=ncol(matres))] <- dimnames(matres[tt, ,drop=FALSE])[[1]]
names(matres2) <- mylab

## with labels
pdf(sprintf("%s/GSEA_res_samecausal_plot_labels.pdf", saveres), width=10, height=5)
par(las=1, mar=c(5, 26, 4, 2) + 0.1)
barplot(height=rev(matres2), space=rev(mysp), horiz=TRUE, col=mycol, xlab="Normalized Enrichment Score", border=NA, xlim=c(-2.5, 2.5))
dev.off()

## without lablels
pdf(sprintf("%s/GSEA_res_samecausal_plot.pdf", saveres), width=5, height=5)
barplot(height=rev(matres2), space=rev(mysp), horiz=TRUE, col=mycol, xlab="Normalized Enrichment Score", border=NA, names.arg="", xlim=c(-2.5, 2.5))
dev.off()

