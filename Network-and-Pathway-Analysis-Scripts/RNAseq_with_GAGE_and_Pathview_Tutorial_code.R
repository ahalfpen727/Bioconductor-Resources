###################################################
###################################################
##step 0: Setup (you don't need to map reads outside R)
###################################################

## set working directory wherever you like
setwd("~/Dropbox/Feina/2015MasterVHIR/RNAseq_Pathways")

## install Bioconductor and packages in R
## load the script from the internet that is used in install bioconductor
source("http://bioconductor.org/biocLite.R")

## commands that tell Bioconductor to download and install each package
biocLite(c("pathview", "gage", "gageData"))
biocLite(c("edgeR", "limma"))
##biocLite(c("GenomicAlignments","TxDb.Hsapiens.UCSC.hg19.knownGene"))

###################################################
### WARNING MESSAGE: to avoid waiting too long time
### answer no if asked for updating other packages
###################################################



###################################################
###################################################
##step 1: Count the reads mapped to each gene
##        (you don't need to run this code)
###################################################

#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#exByGn <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")
#
#library(GenomicAlignments)
#
#fls <- list.files("tophat_all/", pattern="bam$", full.names=T)
#bamfls <- BamFileList(fls)
#
#flag <- scanBamFlag(isNotPrimaryRead=FALSE, isProperPair=TRUE)
#param <- ScanBamParam(flag=flag)
#
#gnCnt <- summarizeOverlaps(exByGn, bamfls, mode="Union",
#                           ignore.strand=TRUE, single.end=FALSE, param=param)
#hnrnp.cnts <- assay(gnCnt)

###################################################
## All commands from previous chunk (step 1), as well as tophat alignment... 
## can be replaced, in this tutorial, by the command:
# data(hnrnp.cnts)
## in order to load the counts table and avoid performing all the hard work ;)
###################################################


###################################################
###################################################
##step 2: Normalize and process read counts
###################################################

library(gageData)
data(hnrnp.cnts)
cnts <- hnrnp.cnts
dim(cnts)

## check elements from counts table with rowSums value is not zero
sel.rn <- rowSums(cnts) != 0

## subset counts table taking only those cases with rowSums value not zero
cnts <- cnts[sel.rn,]
dim(cnts)

############################################################
##joint workflow with DEseq/edgeR/limma/Cufflinks forks here
############################################################

## compute library size for each sample (each column)
libsizes <- colSums(cnts)
show(libsizes)

## plot lib sizes
pdf("hnrnp.lib.sizes.pdf")
barplot(libsizes, ylim=c(0, 4e+07), las=3, cex.names=.7, col="green",
        main="Library size per sample")
dev.off()

## compute size factor in order to normalize by lib size
size.factor <- libsizes/exp(mean(log(libsizes)))
show(size.factor)

## normalize by using size factors and log2 transform
cnts.norm <- t(t(cnts)/size.factor)
range(cnts.norm)
head(cnts.norm)

## add a small yet appropriate positive count number (+8) to all genes before
## doing log2 transformation, in order to avoid -Inf and stablize the variance 
## at low expression end
cnts.norm <- log2(cnts.norm+8)

dim(cnts.norm)
range(cnts.norm)
head(cnts.norm)
tail(cnts.norm)

## generate optional MA plot
pdf("hnrnp.cnts.maplot.pdf")
plot((cnts.norm[,1]+cnts.norm[,5])/2, (cnts.norm[,1]-cnts.norm[,5]),
     main="MA plot Knockdown vs Control", xlab="mean", ylab="change",
     ylim=c(-5,5), xlim=c(0,20), lwd=1)
abline(h=0, lwd=2, col="red", lty="dashed")
dev.off()


###################################################
###################################################
##step 3: Differential expression with edgeR and Limma
###################################################
## Here we assume that we have gone through the 
## earlier steps till the fork point right after 
## the code line "cnts <- cnts[sel.rn,]" that you 
## can find in the previous chunk (step 2)

library(edgeR)

## define sample conditions for sample 1 to 8
grp.idx <- rep(c("knockdown", "control"), each=4)
grp.idx

## obtain differential expression list
dge <- DGEList(counts=cnts, group=factor(grp.idx))
dge <- calcNormFactors(dge)

library(limma)

## define the model and fit it
design <- model.matrix(~grp.idx)
log2.cpm <- voom(dge, design)
fit <- lmFit(log2.cpm, design)
fit <- eBayes(fit)

## obtain results as a top table (with logFC, pvalues...)
limma.res <- topTable(fit, coef=2, n=Inf, sort="p")
dim(limma.res)
head(limma.res) # notice that row names will be Gene IDs
tail(limma.res) # notice that row names will be Gene IDs

## write out the top table of DE genes
write.csv(limma.res, file="topTable.csv")

#limma.fc <- limma.res$logFC
#names(limma.fc) <- limma.res$ID
#exp.fc <- limma.fc
#exp.fc
#out.suffix <- "limma"


###################################################
###################################################
##step 4: Gene set test with GAGE
###################################################

## joint workflow with DEseq/edgeR/limma/Cufflinks merges around here
library(gage)
ref.idx <- 5:8
samp.idx <- 1:4
data(kegg.gs)
## knockdown and control samples are unpaired
cnts.kegg.p <- gage(cnts.norm, gsets=kegg.gs, ref=ref.idx,
                    samp=samp.idx, compare="unpaired")


###################################################
###################################################
##step 5: Pathway visualization with Pathview
###################################################

## differential expression: log2 ratio or fold change, unpaired samples
cnts.d <-  cnts.norm[, samp.idx] - rowMeans(cnts.norm[, ref.idx])
dim(cnts.d)
head(cnts.d)
tail(cnts.d)

## select up-regulated pathways (q-value < 0.1)
sel.up <- cnts.kegg.p$greater[, "q.val"] < 0.1 & !is.na(cnts.kegg.p$greater[,"q.val"])

## obtain Ids of top up-regulated pathways
path.ids.up <- rownames(cnts.kegg.p$greater)[sel.up]
path.ids.up
path.ids2.up <- substr(path.ids.up, 1, 8)

## up-regulated pathways (top 3) generated by pathview
library(pathview)
pv.out.list.up <- sapply(path.ids2.up[1:3], function(pid) pathview(
  gene.data = cnts.d, pathway.id = pid,
  species = "hsa"))

## select down-regulated pathways (q-value < 0.1)
sel.dwn <- cnts.kegg.p$less[, "q.val"] < 0.1 & !is.na(cnts.kegg.p$less[,"q.val"])

## obtain Ids of top down-regulated pathways
path.ids.dwn <- rownames(cnts.kegg.p$less)[sel.dwn]
path.ids.dwn
path.ids2.dwn <- substr(path.ids.dwn, 1, 8)

## down-regulated pathways (top 3) visualized by pathview
pv.out.list.dwn <- sapply(path.ids2.dwn[1:3], function(pid) pathview(
  gene.data = cnts.d, pathway.id = pid,
  species = "hsa"))


###################################################
###################################################
##step 6: Analyses of GO terms
###################################################
library(gageData)
data(go.sets.hs)
data(go.subs.hs)
lapply(go.subs.hs, head)

## Molecular Function (MF) analysis is quicker, hence run as demo
cnts.mf.p <- gage(cnts.norm, gsets = go.sets.hs[go.subs.hs$MF],
                  ref = ref.idx, samp = samp.idx, compare ="unpaired")

## Biological Process (BP) takes a few minutes if you try it
#cnts.bp.p <- gage(cnts.norm, gsets = go.sets.hs[go.subs.hs$BP], 
#                  ref = ref.idx, samp = samp.idx, compare ="unpaired")

## Loop to analyze (and plot, as pdf) expression perturbations in the 
## most significant GO terms
for (gs in rownames(cnts.mf.p$less)[1:3]) {
    outname = gsub(" |:|/", "_", substr(gs, 12, 100))
    geneData(genes = go.sets.hs[[gs]], exprs = cnts.norm, ref = ref.idx,
             samp = samp.idx, outname = outname, txt = T, heatmap = T,
             limit = 3, scatterplot = T)
}




###################################################
###################################################
###################################################
