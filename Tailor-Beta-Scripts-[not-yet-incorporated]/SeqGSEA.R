### R code from vignette source 'SeqGSEA.Rnw'
source("https://bioconductor.org/biocLite.R")
biocLite("GSVAdata")
library(GSAR)


###################################################
### code chunk number 1: SeqGSEA
###################################################
library(SeqGSEA)


###################################################
### code chunk number 2: help (eval = FALSE)
###################################################
## ? SeqGSEA


###################################################
### code chunk number 3: ReadCountSet
###################################################
rcounts <- cbind(t(sapply(1:10, function(x) {rnbinom(5, size=10, prob=runif(1))} )), 
                 t(sapply(1:10, function(x) {rnbinom(5, size=10, prob=runif(1))} )))
colnames(rcounts) <- c(paste("S", 1:5, sep=""), paste("C", 1:5, sep="")) 
geneIDs <- c(rep("G1", 4), rep("G2", 6))
exonIDs <- c(paste("E", 1:4, sep=""), paste("E", 1:6, sep=""))
RCS <- newReadCountSet(rcounts, exonIDs, geneIDs)
RCS 


###################################################
### code chunk number 4: RCS_example
###################################################
data(RCS_example, package="SeqGSEA")
RCS_example


###################################################
### code chunk number 5: geneID
###################################################
length(unique(geneID(RCS_example)))


###################################################
### code chunk number 6: exonTestability
###################################################
RCS_example <- exonTestability(RCS_example, cutoff = 5)


RCS_example <- estiExonNBstat(RCS_example)
RCS_example <- estiGeneNBstat(RCS_example)
head(fData(RCS_example)[, c("exonIDs", "geneIDs", "testable", "NBstat")])


permuteMat <- genpermuteMat(RCS_example, times=20)
RCS_example <- DSpermute4GSEA(RCS_example, permuteMat)
head(RCS_example@permute_NBstat_gene)
DSscore.normFac <- normFactor(RCS_example@permute_NBstat_gene)
DSscore <- scoreNormalization(RCS_example@featureData_gene$NBstat, 
                              DSscore.normFac)
DSscore.perm <- scoreNormalization(RCS_example@permute_NBstat_gene, 
                                   DSscore.normFac)
DSscore[1:5]
DSscore.perm[1:5,1:10]


###################################################
### code chunk number 10: DSpval
###################################################
RCS_example <- DSpermutePval(RCS_example, permuteMat)
head(DSresultGeneTable(RCS_example))


###################################################
### code chunk number 11: genecount
###################################################
geneCounts <- getGeneCount(RCS_example)
dim(geneCounts) # 182  20
head(geneCounts)


###################################################
### code chunk number 12: DEseq
###################################################
label <- label(RCS_example)
DEG <- runDESeq(geneCounts, label)


###################################################
### code chunk number 13: DSNB
###################################################
DEGres <- DENBStat4GSEA(DEG)
DEGres[1:5, "NBstat"]


###################################################
### code chunk number 14: DSNBpermut
###################################################
DEpermNBstat <- DENBStatPermut4GSEA(DEG, permuteMat)
DEpermNBstat[1:5, 1:10]


###################################################
### code chunk number 15: DEscores
###################################################
DEscore.normFac <- normFactor(DEpermNBstat)
DEscore <- scoreNormalization(DEGres$NBstat, DEscore.normFac)
DEscore.perm <- scoreNormalization(DEpermNBstat, DEscore.normFac)
DEscore[1:5]
DEscore.perm[1:5, 1:10]


###################################################
### code chunk number 16: DEpval
###################################################
DEGres <- DEpermutePval(DEGres, DEpermNBstat) 
DEGres[1:6, c("NBstat", "perm.pval", "perm.padj")]


###################################################
### code chunk number 17: DEpval2
###################################################
DEGres <- DENBTest(DEG)
DEGres <- DEpermutePval(DEGres, DEpermNBstat) 
DEGres[1:6, c("NBstat", "pval", "padj", "perm.pval", "perm.padj")]


###################################################
### code chunk number 18: genescore_l
###################################################
gene.score <- geneScore(DEscore, DSscore, method="linear", DEweight = 0.3)
gene.score.perm <- genePermuteScore(DEscore.perm, DSscore.perm, 
                                    method="linear",  DEweight=0.3)
plotGeneScore(gene.score, gene.score.perm)


###################################################
### code chunk number 19: genescore_r
###################################################
gene.score <- geneScore(DEscore, DSscore, method="rank", DEweight = 0.3)
gene.score.perm <- genePermuteScore(DEscore.perm, DSscore.perm, 
                                    method="rank",  DEweight=0.3)
plotGeneScore(gene.score, gene.score.perm)


###################################################
### code chunk number 20: genescore_gr
###################################################
combine <- rankCombine(DEscore, DSscore, DEscore.perm, DSscore.perm, DEweight=0.3) 
gene.score <- combine$geneScore
gene.score.perm <- combine$genePermuteScore
plotGeneScore(gene.score, gene.score.perm)


###################################################
### code chunk number 21: seqgeneset
###################################################
data(GS_example, package="SeqGSEA") 
GS_example


###################################################
### code chunk number 22: GSEAmain
###################################################
GS_example <- GSEnrichAnalyze(GS_example, gene.score, gene.score.perm)
topGeneSets(GS_example, 5)


###################################################
### code chunk number 23: plotES
###################################################
plotES(GS_example)


###################################################
### code chunk number 24: plotSig
###################################################
plotSig(GS_example)


###################################################
### code chunk number 25: plotSigGS
###################################################
plotSigGeneSet(GS_example, 9, gene.score) # 9th gene set is the most significant one.


###################################################
### code chunk number 26: wrightSigGS
###################################################
writeSigGeneSet(GS_example, 9, gene.score) # 9th gene set is the most significant one.


###################################################
### code chunk number 27: doParallel1
###################################################
library(doParallel)
a <- matrix(1:16, 4, 4)
b <- t(a)
foreach(b=iter(b, by='col'), .combine=cbind) %dopar%
  (a %*% b)


###################################################
### code chunk number 28: doParallel2
###################################################
library(doParallel)
cl <- makeCluster(2) # specify 2 cores to be used in this computing 
registerDoParallel(cl)
getDoParWorkers() # 2
a <- matrix(1:16, 4, 4)
b <- t(a)
foreach(b=iter(b, by='col'), .combine=cbind) %dopar%
  (a %*% b)


###################################################
### code chunk number 29: sysfile
###################################################
system.file("extscripts", package="SeqGSEA", mustWork=TRUE)


###################################################
### code chunk number 30: Initialization
###################################################
rm(list=ls())
# input count data files
data.dir <- system.file("extdata", package="SeqGSEA", mustWork=TRUE)
case.pattern <- "^SC"  # file name starting with "SC"
ctrl.pattern <- "^SN"  # file name starting with "SN"
case.files <- dir(data.dir, pattern=case.pattern, full.names = TRUE)
control.files <- dir(data.dir, pattern=ctrl.pattern, full.names = TRUE)
# gene set file
geneset.file <- system.file("extdata", "gs_symb.txt",  
                            package="SeqGSEA", mustWork=TRUE)
# output file prefix 
output.prefix <- "SeqGSEA.test"
# setup parallel backend 
library(doParallel)
cl <- makeCluster(2) # specify 2 cores to be used in computing 
registerDoParallel(cl) # parallel backend registration
# setup permutation times
perm.times <- 20 # change the number to >= 1000 in your analysis 


###################################################
### code chunk number 31: DS_analysis
###################################################
# load exon read count data
RCS <- loadExonCountData(case.files, control.files)
# remove genes with low expression 
RCS <- exonTestability(RCS, cutoff=5)
geneTestable <- geneTestability(RCS)
RCS <- subsetByGenes(RCS, unique(geneID(RCS))[ geneTestable ])
# get gene IDs, which will be used in initialization of gene set
geneIDs <- unique(geneID(RCS))
# calculate DS NB statistics
RCS <- estiExonNBstat(RCS)
RCS <- estiGeneNBstat(RCS)
# calculate DS NB statistics on the permutation data sets
permuteMat <- genpermuteMat(RCS, times=perm.times)
RCS <- DSpermute4GSEA(RCS, permuteMat)


###################################################
### code chunk number 32: DE_analysis
###################################################
# get gene read counts
geneCounts <- getGeneCount(RCS)
# calculate DE NB statistics 
label <- label(RCS)
DEG <-runDESeq(geneCounts, label)
DEGres <- DENBStat4GSEA(DEG)
# calculate DE NB statistics on the permutation data sets
DEpermNBstat <- DENBStatPermut4GSEA(DEG, permuteMat) # permutation


###################################################
### code chunk number 33: score_int
###################################################
# DE score normalization
DEscore.normFac <- normFactor(DEpermNBstat)
DEscore <- scoreNormalization(DEGres$NBstat, DEscore.normFac)
DEscore.perm <- scoreNormalization(DEpermNBstat, DEscore.normFac)
# DS score normalization
DSscore.normFac <- normFactor(RCS@permute_NBstat_gene)
DSscore <- scoreNormalization(RCS@featureData_gene$NBstat, DSscore.normFac)
DSscore.perm <- scoreNormalization(RCS@permute_NBstat_gene, DSscore.normFac)
# score integration
gene.score <- geneScore(DEscore, DSscore, DEweight=0.5)
gene.score.perm <- genePermuteScore(DEscore.perm, DSscore.perm, DEweight=0.5)
# visilization of scores 
# NOT run in the example; users to uncomment the following 6 lines to run 
#plotGeneScore(DEscore, DEscore.perm, pdf=paste(output.prefix,".DEScore.pdf",sep=""),
#              main="Expression")
#plotGeneScore(DSscore, DSscore.perm, pdf=paste(output.prefix,".DSScore.pdf",sep=""), 
#              main="Splicing")
#plotGeneScore(gene.score, gene.score.perm, 
#              pdf=paste(output.prefix,".GeneScore.pdf",sep=""))


###################################################
### code chunk number 34: main_gsea
###################################################
# load gene set data
gene.set <- loadGenesets(geneset.file, geneIDs, geneID.type="ensembl",
                         genesetsize.min = 5, genesetsize.max = 1000)
# enrichment analysis
gene.set <- GSEnrichAnalyze(gene.set, gene.score, gene.score.perm, weighted.type=1)
# format enrichment analysis results
GSEAres <- GSEAresultTable(gene.set, TRUE)
# output results 
# NOT run in the example; users to uncomment the following 4 lines to run 
#write.table(GSEAres, paste(output.prefix,".GSEA.result.txt",sep=""), 
#            quote=FALSE, sep="\t", row.names=FALSE)
#plotES(gene.set, pdf=paste(output.prefix,".GSEA.ES.pdf",sep=""))
#plotSig(gene.set, pdf=paste(output.prefix,".GSEA.FDR.pdf",sep=""))


###################################################
### code chunk number 35: Initialization_DE
###################################################
rm(list=ls())
# input count data files
data.dir <- system.file("extdata", package="SeqGSEA", mustWork=TRUE)
count.file <- paste(data.dir,"geneCounts.txt",sep="/")
# gene set file
geneset.file <- system.file("extdata", "gs_symb.txt",  
                            package="SeqGSEA", mustWork=TRUE)
# output file prefix 
output.prefix <- "SeqGSEA.test"
# setup parallel backend 
library(doParallel)
cl <- makeCluster(2) # specify 2 cores to be used in computing 
registerDoParallel(cl) # parallel backend registration
# setup permutation times
perm.times <- 20 # change the number to >= 1000 in your analysis 


###################################################
### code chunk number 36: DE_analysis_DE
###################################################
# load gene read count data
geneCounts <- read.table(count.file)
# speficify the labels of each sample
label <- as.factor(c(rep(1,10), rep(0,10)))
# calculate DE NB statistics 
DEG <-runDESeq(geneCounts, label)
DEGres <- DENBStat4GSEA(DEG)
# calculate DE NB statistics on the permutation data sets
permuteMat <- genpermuteMat(label, times=perm.times)
DEpermNBstat <- DENBStatPermut4GSEA(DEG, permuteMat) # permutation


###################################################
### code chunk number 37: score_int_DE
###################################################
# DE score normalization
DEscore.normFac <- normFactor(DEpermNBstat)
DEscore <- scoreNormalization(DEGres$NBstat, DEscore.normFac)
DEscore.perm <- scoreNormalization(DEpermNBstat, DEscore.normFac)
# score integration - DSscore can be null 
gene.score <- geneScore(DEscore, DEweight=1)
gene.score.perm <- genePermuteScore(DEscore.perm, DEweight=1)  # visilization of scores 
# NOT run in the example; users to uncomment the following 6 lines to run 
#plotGeneScore(DEscore, DEscore.perm, pdf=paste(output.prefix,".DEScore.pdf",sep=""),
#              main="Expression")
#plotGeneScore(gene.score, gene.score.perm, 
#              pdf=paste(output.prefix,".GeneScore.pdf",sep=""))


###################################################
### code chunk number 38: main_gsea_DE
###################################################
# load gene set data
geneIDs <- rownames(geneCounts)
gene.set <- loadGenesets(geneset.file, geneIDs, geneID.type="ensembl",
                         genesetsize.min = 5, genesetsize.max = 1000)
# enrichment analysis
gene.set <- GSEnrichAnalyze(gene.set, gene.score, gene.score.perm, weighted.type=1)
# format enrichment analysis results
GSEAres <- GSEAresultTable(gene.set, TRUE)
# output results 
# NOT run in the example; users to uncomment the following 4 lines to run 
#write.table(GSEAres, paste(output.prefix,".GSEA.result.txt",sep=""), 
#            quote=FALSE, sep="\t", row.names=FALSE)
#plotES(gene.set, pdf=paste(output.prefix,".GSEA.ES.pdf",sep=""))
#plotSig(gene.set, pdf=paste(output.prefix,".GSEA.FDR.pdf",sep=""))


###################################################
### code chunk number 39: runSeqGSEA
###################################################
### Initialization ###
# input file location and pattern
data.dir <- system.file("extdata", package="SeqGSEA", mustWork=TRUE)
case.pattern <- "^SC" # file name starting with "SC"
ctrl.pattern <- "^SN" # file name starting with "SN"
# gene set file and type
geneset.file <- system.file("extdata", "gs_symb.txt",
                            package="SeqGSEA", mustWork=TRUE)
geneID.type <- "ensembl"
# output file prefix
output.prefix <- "SeqGSEA.example"
# analysis parameters
nCores <- 8
perm.times <- 1000 # >= 1000 recommended
DEonly <- FALSE
DEweight <- c(0.2, 0.5, 0.8)  # a vector for different weights
integrationMethod <- "linear"

### one step SeqGSEA running ###
# NOT run in the example; uncomment the following 4 lines to run 
# CAUTION: running the following lines will generate lots of files in your working dir
#runSeqGSEA(data.dir=data.dir, case.pattern=case.pattern, ctrl.pattern=ctrl.pattern, 
#           geneset.file=geneset.file, geneID.type=geneID.type, output.prefix=output.prefix,
#           nCores=nCores, perm.times=perm.times, integrationMethod=integrationMethod,
#           DEonly=DEonly, DEweight=DEweight)


###################################################
### code chunk number 40: sessionInfo
###################################################
sessionInfo()


###################################################
### code chunk number 41: <closeConnetions
###################################################
allCon <- showConnections()
socketCon <- as.integer(rownames(allCon)[allCon[, "class"] == "sockconn"])
sapply(socketCon, function(ii) close.connection(getConnection(ii)) )


