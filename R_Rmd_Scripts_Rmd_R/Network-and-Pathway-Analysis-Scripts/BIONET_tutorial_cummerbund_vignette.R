### R code from vignette source 'Tutorial.Rnw'
#source("http://bioconductor.org/biocLite.R")
#biocLite("GOstats")
biocLite("GEOquery")
biocLite("affy")
biocLite("gcrma")
biocLite("hugene10stv1cdf")
biocLite("hugene10stv1probe")
biocLite("hugene10stprobeset.db")
biocLite("hugene10sttranscriptcluster.db")
###################################################
### code chunk number 2: Tutorial.Rnw:55-59
###################################################
library(xtable)
library(cummeRbund)
library(BioNet)
library(DLBCL)
data(dataLym)
data(interactome)
str(dataLym)
colnames(dataLym)
rownames(dataLym)


#Load the necessary libraries
library(GEOquery)
library(affy)
library(gcrma)
library(hugene10stv1cdf)
library(hugene10stv1probe)
library(hugene10stprobeset.db)
library(hugene10sttranscriptcluster.db)

#Set working directory for download
setwd("/Users/ogriffit/Dropbox/BioStars")


###################################################
##################################################
 geneListString='CXCL12 TGFB MYC RAS CXCR4 IL8 IL6'
 cuff<-readCufflinks(gtfFile="cuffcmp.combined.gtf",genome="genome.fa",rebuild=F)
 cuff
 refgtf="cuffcmp.combined.gtf"
 genome_path="hg19"
 over="LUTS"
 under="CTRL"
###########################################################################
### Chunk_2: Significantly differentially expressed features: method 1
############################################################################
mySigGenes<-getSig(cuff,x=over,y=under,alpha=.05,level='genes')
head(mySigGenes)
length(mySigGenes)
mySigGeneTable<-getSigTable(cuff,alpha=0.05,level='genes')
head(mySigGeneTable,20)
length(mySigGeneTable)
sigGenes<-getGenes(cuff, mySigGenes)
length(sigGenes)
sig_genes_exp.diff<-diffData(sigGenes)
sig_genes_annot<-sigGenes@annotation
head(sig_genes_exp.diff)
head(sig_genes_annot)
sig_genes_exp.diff["gene_id"]<-sig_genes_annot["gene_short_name"]
head(sig_genes_exp.diff)

mySigCDS<-getSig(cuff,x=over,y=under,alpha=.05,level='CDS')
head(mySigCDS)
length(mySigCDS)
mySigCDSTable<-getSigTable(cuff,alpha=0.05,level='CDS')
head(mySigCDSTable,20)
length(mySigCDSTable)
sigCDS<-getGenes(cuff, mySigCDS)
sigCDS
length(sigCDS)
sig_cds_exp.diff<-diffData(sigCDS)
sig_cds_annot<-sigCDS@annotation
head(sig_cds_exp.diff)
head(sig_cds_annot)
sig_cds_exp.diff["gene_id"]<-sig_cds_annot["gene_short_name"]
head(sig_cds_exp.diff)

mySigIsos<-getSig(cuff,x=over,y=under,alpha=.05,level='isoforms')
head(mySigIsos)
length(mySigIsos)
mySigIsoTable<-getSigTable(cuff,alpha=0.05,level='isoforms')
head(mySigIsoTable,20)
length(mySigIsoTable)
sigIsos<-getGenes(cuff, mySigIsos)
sigIsos
length(sigIsos)
sig_isos_exp.diff<-diffData(sigIsos)
sig_isos_annot<-sigIsos@annotation
head(sig_isos_exp.diff)
head(sig_isos_annot)
sig_isos_exp.diff["gene_id"]<-sig_isos_annot["gene_short_name"]
head(sig_isos_exp.diff)

mySigTSS<-getSig(cuff,x=over,y=under,alpha=.05,level='TSS')
head(mySigTSS)
length(mySigTSS)
mySigTSS<-getSig(cuff,x=over,y=under,alpha=.05,level='TSS')
head(mySigTSS)
length(mySigTSS)
mySigTSSTable<-getSigTable(cuff,alpha=0.05,level='TSS')
head(mySigTSSTable,20)
length(mySigTSSTable)
sigTSS<-getGenes(cuff, mySigTSS)
sigTSS
length(sigTSS)
sig_tss_exp.diff<-diffData(sigTSS)
sig_tss_annot<-sigTSS@annotation
head(sig_tss_exp.diff)
head(sig_tss_annot)
#################################################################
### code chunk number 8: repFPKMMatrix
##################################################################

genediff<-diffTable(genes(cuff))
head(genediff)
tssdiff<-diffTable(TSS(cuff))
head(tssdiff)
isodiff<-diffTable(isoforms(cuff), )
head(isodiff)
cdsdiff<-diffTable(CDS(cuff))
head(cdsdiff)
###################################################
runInfo(cuff)
idpatient<-replicates(cuff)
conditions(cuff)
colnames(idpatient)
sample_digit<-idpatient[,"file"]
g.rep.matrix<-repFpkmMatrix(genes(cuff))
head(g.rep.matrix)
colnames(g.rep.matrix)<-idpatient[,"file"]
head(g.rep.matrix)

gene_scale_table<-rbind2(c(idpatient[,"total_mass"],idpatient[,"norm_mass"],
						  idpatient[,"internal_scale"], idpatient[,"rep_name"], idpatient[,"file"]), g.rep.matrix)
head(gene_scale_table)

gene.xloc.matrix<-featureNames(genes(cuff))
gene.list<-getGenes(cuff,geneId = gene.xloc.matrix)
gene.list
gene_annotation_data<-featureNames(gene.list)
head(gene_annotation_data)
gene.rep.matrix<-cbind2(gene_annotation_data, g.rep.matrix)
gene_table = file.path("./RepFpkmMatrix_leaveout19.genes.txt")
#gene_table = file.path(Sys.getenv("CUMMERBUND_FILE"))
write.table(gene.rep.matrix,gene_table, sep="  ", row.names = F , col.names = T,quote = F)
head(gene.rep.matrix)
length(colnames(gene.rep.matrix))
gene_rep.df<-gene.rep.matrix[,2:20]
gene_rep.df<-as.data.frame(gene_rep.df, colnames=T, rownames=F, sep="\t")
dim(gene_rep.df)
head(gene_rep.df)
siggenes_table<-as.data.frame(sig_tss_exp.diff, colnames=T,
							  rownames=F, sep="\t")
head(siggenes_table)
or<-sort(siggenes_table[,"q_value"], decreasing=F)
siggenes_table[,"gene_id"]<-sig_tss_annot[,"gene_short_name"]
head(siggenes_table)
###################################################
### code chunk number 3: Tutorial.Rnw:63-66
###################################################
pvals <- cbind(t=dataLym$t.pval, s=dataLym$s.pval)
rownames(pvals) <- dataLym$label
pval <- aggrPvals(pvals, order=2, plot=FALSE)
head(dataLym)
head(pvals)
class(pval)
pvals <- cbind(t=siggenes_table$p.pval, s=siggenes_table$q.pval)
#rownames(pvals) <- siggenes_table$gene_id
#pval <- aggrPvals(pvals, order=2, plot=FALSE)

pval<-siggenes_table[,"q_value"]
pval<-cbind2(siggenes_table["gene_id"], pval)
length(pval)
colnames(pval)<-c("Gene_Name", "P_Value")
head(pval)
###################################################
### code chunk number 4: Tutorial.Rnw:70-73
###################################################
subnet <- subNetwork(dataLym$label, interactome)
subnet <- rmSelfLoops(subnet)
subnet

lutsubnet<-subNetwork(genediff$gene_short_name, interactome)
luts_subnet<-rmSelfLoops(lutsubnet)
luts_subnet

###################################################
### code chunk number 5: Tutorial.Rnw:79-81
###################################################
fb <- fitBumModel(pval, plot=FALSE)
scores <- scoreNodes(subnet, fb, fdr=0.001)


###################################################
### code chunk number 6: Tutorial.Rnw:88-91
###################################################
module <- runFastHeinz(subnet, scores)
logFC <- dataLym$diff
names(logFC) <- dataLym$label


###################################################
### code chunk number 7: Tutorial.Rnw:98-99
###################################################
plotModule(module, scores=scores, diff.expr=logFC)
exprLym


###################################################
### code chunk number 10: Tutorial.Rnw:173-174
###################################################
interactome


###################################################
### code chunk number 11: Tutorial.Rnw:181-183
###################################################
network <- subNetwork(featureNames(exprLym), interactome)
network


###################################################
### code chunk number 12: Tutorial.Rnw:188-190
###################################################
network <- largestComp(network)
network


###################################################
### code chunk number 13: Tutorial.Rnw:200-204
###################################################
###################################################
### code chunk number 8: Tutorial.Rnw:154-158
###################################################
library(BioNet)
library(DLBCL)
data(exprLym)
data(interactome)
library(genefilter)
library(impute)
expressions <- impute.knn(exprs(exprLym))$data
t.test <- rowttests(expressions, fac=exprLym$Subgroup)


###################################################
### code chunk number 14: Tutorial.Rnw:207-208
###################################################
t.test[1:10, ]


###################################################
### code chunk number 15: Tutorial.Rnw:215-218
###################################################
top.table <- xtable(t.test[1:10,], display=c("s", "f", "f", "f"))
print(top.table, floating=FALSE)


###################################################
### code chunk number 16: Tutorial.Rnw:228-233
###################################################
data(dataLym)
ttest.pval <- t.test[, "p.value"]
surv.pval <- dataLym$s.pval
names(surv.pval) <- dataLym$label
pvals <- cbind(ttest.pval, surv.pval)


###################################################
### code chunk number 17: Tutorial.Rnw:242-243
###################################################
pval <- aggrPvals(pvals, order=2, plot=FALSE)


###################################################
### code chunk number 18: Tutorial.Rnw:250-252
###################################################
fb <- fitBumModel(pval, plot=FALSE)
fb


###################################################
### code chunk number 19: Tutorial.Rnw:255-260
###################################################
dev.new(width=13, height=7)
par(mfrow=c(1,2))
hist(fb)
plot(fb)
dev.off()


###################################################
### code chunk number 20: Tutorial.Rnw:275-276
###################################################
plotLLSurface(pval, fb)


###################################################
### code chunk number 21: Tutorial.Rnw:285-286
###################################################
scores <- scoreNodes(network=network, fb=fb, fdr=0.001)


###################################################
### code chunk number 22: Tutorial.Rnw:293-296
###################################################
network <- rmSelfLoops(network)
writeHeinzEdges(network=network, file="lymphoma_edges_001", use.score=FALSE)
writeHeinzNodes(network=network, file="lymphoma_nodes_001", node.scores = scores)


###################################################
### code chunk number 23: Tutorial.Rnw:313-315
###################################################
datadir <- file.path(path.package("BioNet"), "extdata")
dir(datadir)


###################################################
### code chunk number 24: Tutorial.Rnw:321-324
###################################################
module <- readHeinzGraph(node.file=file.path(datadir, "lymphoma_nodes_001.txt.0.hnz"), network=network)
diff <- t.test[, "dm"]
names(diff) <- rownames(t.test)


###################################################
### code chunk number 25: Tutorial.Rnw:327-328
###################################################
plotModule(module, diff.expr=diff, scores=scores)


###################################################
### code chunk number 26: Tutorial.Rnw:350-353
###################################################
sum(scores[nodes(module)])
sum(scores[nodes(module)]>0)
sum(scores[nodes(module)]<0)

# install the GEO libraries
biocLite("GEOquery")
library(GEOquery)
getGEOSuppFiles("GSE20986")
untar("GSE20986/GSE20986_RAW.tar", exdir="data")
cels <- list.files("data/", pattern = "[gz]")
sapply(paste("data", cels, sep="/"), gunzip)
cels
library(simpleaffy)
celfiles<-read.affy(covdesc="phenodata.txt", path="data")
celfiles
celfiles.gcrma<-gcrma(celfiles)
celfiles.gcrma
# load colour libraries
library(RColorBrewer)
# set colour palette
cols <- brewer.pal(8, "Set1")
# plot a boxplot of unnormalised intensity values
boxplot(celfiles, col=cols)
# plot a boxplot of normalised intensity values, affyPLM required to interrogate celfiles.gcrma
library(affyPLM)
boxplot(celfiles.gcrma, col=cols)
# the boxplots are somewhat skewed by the normalisation algorithm
# and it is often more informative to look at density plots
# Plot a density vs log intensity histogram for the unnormalised data
hist(celfiles, col=cols)
# Plot a density vs log intensity histogram for the normalised data
hist(celfiles.gcrma, col=cols)
# Perform probe-level metric calculations on the CEL files:
celfiles.qc <- fitPLM(celfiles)
# Create an image of GSM24662.CEL:
image(celfiles.qc, which=1, add.legend=TRUE)
# Create an image of GSM524665.CEL
# There is a spatial artifact present
image(celfiles.qc, which=4, add.legend=TRUE)
# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.  GSM524665.CEL is an outlier
RLE(celfiles.qc, main="RLE")
# We can also use NUSE (Normalised Unscaled Standard Errors).
# The median standard error should be 1 for most genes.
# GSM524665.CEL appears to be an outlier on this plot too
NUSE(celfiles.qc, main="NUSE")
eset <- exprs(celfiles.gcrma)
distance <- dist(t(eset),method="maximum")
clusters <- hclust(distance)
plot(clusters)
celfiles.filtered <- nsFilter(celfiles.gcrma, require.entrez=FALSE, remove.dupEntrez=FALSE)
# What got removed and why?
celfiles.filtered$filter.log
samples <- celfiles.gcrma$Target
# check the results of this
samples
# convert into factors
samples<- as.factor(samples)
# check factors have been assigned
samples
# set up the experimental design
design = model.matrix(~0 + samples)
colnames(design) <- c("choroid", "huvec", "iris", "retina")
# inspect the experiment design
design
library(limma)
# fit the linear model to the filtered expression set
fit <- lmFit(exprs(celfiles.filtered$eset), design)
# set up a contrast matrix to compare tissues v cell line
contrast.matrix <- makeContrasts(huvec_choroid = huvec - choroid, huvec_retina = huvec - retina, huvec_iris = huvec - iris, levels=design)
# check the contrast matrix
contrast.matrix
# Now the contrast matrix is combined with the per-probeset linear model fit.
huvec_fits <- contrasts.fit(fit, contrast.matrix)
huvec_ebFit <- eBayes(huvec_fits)
# return the top 10 results for any given contrast
# coef=1 is huvec_choroid, coef=2 is huvec_retina
topTable(huvec_ebFit, number=10, coef=1)
nrow(topTable(huvec_ebFit, coef=1, number=10000, lfc=5))
nrow(topTable(huvec_ebFit, coef=1, number=10000, lfc=4))
nrow(topTable(huvec_ebFit, coef=1, number=10000, lfc=3))
nrow(topTable(huvec_ebFit, coef=1, number=10000, lfc=2))
# Get a list for probesets with a four fold change or more
probeset.list <- topTable(huvec_ebFit, coef=1, number=10000, lfc=4)
biocLite("hgu133plus2.db")
library(hgu133plus2.db)
library(annotate)
gene.symbols <- getSYMBOL(probeset.list$ID, "hgu133plus2")
results <- cbind(probeset.list, gene.symbols)
write.table(results, "results.txt", sep="\t", quote=FALSE)
# install additional bioconductor libraries, if not already installed

#Download the CEL file package for this dataset (by GSE - Geo series id)
getGEOSuppFiles("GSE27447")

#Unpack the CEL files
setwd("/Users/ogriffit/Dropbox/BioStars/GSE27447")
untar("GSE27447_RAW.tar", exdir="data")
cels = list.files("data/", pattern = "CEL")
sapply(paste("data", cels, sep="/"), gunzip)
cels = list.files("data/", pattern = "CEL")

setwd("/Users/ogriffit/Dropbox/BioStars/GSE27447/data")
raw.data=ReadAffy(verbose=TRUE, filenames=cels, cdfname="hugene10stv1") #From bioconductor

#perform RMA normalization (I would normally use GCRMA but it did not work with this chip)
data.rma.norm=rma(raw.data)

#Get the important stuff out of the data - the expression estimates for each array
rma=exprs(data.rma.norm)

#Format values to 5 decimal places
rma=format(rma, digits=5)

#Map probe sets to gene symbols or other annotations
#To see all available mappings for this platform
ls("package:hugene10stprobeset.db") #Annotations at the exon probeset level
ls("package:hugene10sttranscriptcluster.db") #Annotations at the transcript-cluster level (more gene-centric view)

#Extract probe ids, entrez symbols, and entrez ids
probes=row.names(rma)
Symbols = unlist(mget(probes, hugene10sttranscriptclusterSYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hugene10sttranscriptclusterENTREZID, ifnotfound=NA))

#Combine gene annotations with raw data
rma=cbind(probes,Symbols,Entrez_IDs,rma)

#Write RMA-normalized, mapped data to file
write.table(rma, file = "rma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
df <- read.table("/path/to/file/rpkm.txt")
    dim(df) #32291   352
    df <- df[,-c(1,2)] # first 2 columns have accessory data

    library(rrcov)
    pcaHub <- PcaHubert(t(df))
    outliers <- which(pcaHub@flag=='FALSE')
###################################################
### code chunk number 29: Tutorial.Rnw:396-397
###################################################
library(BioNet)
library(DLBCL)
library(ALL)
data(ALL)
data(interactome)
interactome
ALL

###################################################
### code chunk number 30: Tutorial.Rnw:407-409
###################################################
mapped.eset <- mapByVar(ALL, network=interactome, attr="geneID")
mapped.eset[1:5,1:5]
length(intersect(rownames(mapped.eset), nodes(interactome)))


###################################################
### code chunk number 32: Tutorial.Rnw:424-429
###################################################
network <- subNetwork(rownames(mapped.eset), interactome)
network
network <- largestComp(network)
network <- rmSelfLoops(network)
network


###################################################
### code chunk number 33: Tutorial.Rnw:441-449
###################################################
library(limma)
design <- model.matrix(~ -1+ factor(c(substr(unlist(ALL$BT), 0, 1))))
colnames(design)<- c("B", "T")
contrast.matrix <- makeContrasts(B-T, levels=design)
contrast.matrix
fit <- lmFit(mapped.eset, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


###################################################
### code chunk number 34: Tutorial.Rnw:454-455
###################################################
pval <- fit2$p.value[,1]


###################################################
### code chunk number 35: Tutorial.Rnw:466-468
###################################################
fb <- fitBumModel(pval, plot=FALSE)
fb


###################################################
### code chunk number 36: Tutorial.Rnw:471-475
###################################################
dev.new(width=13, height=7)
par(mfrow=c(1,2))
hist(fb)
plot(fb)


###################################################
### code chunk number 37: Tutorial.Rnw:487-488
###################################################
scores <- scoreNodes(network=network, fb=fb, fdr=1e-14)


###################################################
### code chunk number 38: Tutorial.Rnw:494-496
###################################################
writeHeinzEdges(network=network, file="ALL_edges_001", use.score=FALSE)
writeHeinzNodes(network=network, file="ALL_nodes_001", node.scores = scores)


###################################################
### code chunk number 39: Tutorial.Rnw:517-519
###################################################
datadir <- file.path(path.package("BioNet"), "extdata")
module <- readHeinzGraph(node.file=file.path(datadir, "ALL_nodes_001.txt.0.hnz"), network=network)


###################################################
### code chunk number 40: Tutorial.Rnw:524-529
###################################################
nodeDataDefaults(module, attr="diff") <- ""
nodeData(module, n=nodes(module), attr="diff") <- fit2$coefficients[nodes(module),1]
nodeDataDefaults(module, attr="score") <- ""
nodeData(module, n=nodes(module), attr="score") <- scores[nodes(module)]
nodeData(module)[1]


###################################################
### code chunk number 41: Tutorial.Rnw:535-536
###################################################
saveNetwork(module, file="ALL_module", type="XGMML")


###################################################
### code chunk number 42: Tutorial.Rnw:568-576 (eval = FALSE)
###################################################
## j.repl <- 100
## resampling.pvals <- list()
## for(i in 1:j.repl)
## {
##   resampling.result <- resamplingPvalues(exprMat=mapped.eset, groups=factor(c(substr(unlist(ALL$BT), 0, 1))), resampleMat=FALSE, alternative="two.sided")
##   resampling.pvals[[i]] <- resampling.result$p.values
##   print(i)
## }


###################################################
### code chunk number 43: Tutorial.Rnw:585-591 (eval = FALSE)
###################################################
## fb <- lapply(resampling.pvals, fitBumModel, plot=FALSE, starts=1)
## resampling.scores <- c()
## for(i in 1:j.repl)
## {
##   resampling.scores[[i]] <- scoreNodes(network=network, fb=fb[[i]], fdr=1e-14)
## }


###################################################
### code chunk number 44: Tutorial.Rnw:596-598 (eval = FALSE)
###################################################
## score.mat <- as.data.frame(resampling.scores)
## colnames(score.mat) <- paste("resample", (1:j.repl), sep="")


###################################################
### code chunk number 45: Tutorial.Rnw:603-605 (eval = FALSE)
###################################################
## writeHeinzEdges(network=network, file="ALL_e_resample", use.score=FALSE)
## writeHeinzNodes(network=network, file="ALL_n_resample", node.scores = score.mat)


###################################################
### code chunk number 46: Tutorial.Rnw:618-620
###################################################
datadir <- file.path(path.package("BioNet"), "extdata")
modules <- readHeinzGraph(node.file=file.path(datadir, "ALL_n_resample.txt.0.hnz"), network=network)


###################################################
### code chunk number 47: Tutorial.Rnw:629-631
###################################################
cons.scores <- consensusScores(modules, network)
writeHeinz(network=network, file="ALL_cons", node.scores=cons.scores$N.scores, edge.scores=cons.scores$E.scores)


###################################################
### code chunk number 48: Tutorial.Rnw:638-643
###################################################
datadir <- file.path(path.package("BioNet"), "extdata")
cons.module <- readHeinzGraph(node.file=file.path(datadir, "ALL_cons_n.txt.0.hnz"), network=network)
cons.edges <- sortedEdgeList(cons.module)
E.width <- 1+cons.scores$E.freq[cons.edges]*10
N.size <- 1+cons.scores$N.freq[nodes(cons.module)]*10


###################################################
### code chunk number 49: Tutorial.Rnw:647-648
###################################################
plotModule(cons.module, edge.width=E.width, vertex.size=N.size, edge.label=cons.scores$E.freq[cons.edges]*100, edge.label.cex=0.6)



