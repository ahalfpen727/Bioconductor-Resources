### R code from vignette source 'Rsubread.Rnw'
source("http://bioconductor.org/biocLite.R")
biocLite(pkgs=c("Rsubread","limma","edgeR"))
### R code from vignette source 'seqc.Rnw'
source("https://bioconductor.org/biocLite.R")
biocLite("seqc")
###################################################
### code chunk number 1: seqc.Rnw:31-34
###################################################
library(seqc)
#options(width=110, digits=2)
#ls(2)
library("Rsubread")
library("limma")
library("edgeR")

###################################################
### code chunk number 3: seqc script
###################################################
dim(ROC_aceview_gene_MGP)
ROC_aceview_gene_MGP[1:15,]
colnames(ILM_aceview_gene_BGI)
ILM_aceview_gene_BGI[1:15,1:7]
ILM_junction_AGR_B[1:15,]
colnames(taqman)
taqman[1:15, 1:9]
entrez_to_symbol<-taqman[,1:2]
###################################################
### code chunk number 1: Rsubread.Rnw:67-70
###################################################

ref <- system.file("extdata","reference.fa",package="Rsubread")
buildindex(basename="reference_index",reference=ref)
featureCounts(files=targets$OutputFile,annot.inbuilt="hg19")
buildindex(basename="chr1",reference="hg19_chr1.fa")
###################################################
### code chunk number 2: Rsubread.Rnw:87-89
###################################################
reads <- system.file("extdata","reads.txt.gz",package="Rsubread")
align(index="reference_index",readfile1=reads,output_file="alignResults.BAM",phredOffset=64)

targets <- readTargets()
align(index="chr1",readfile1=targets$InputFile,input_format="gzFASTQ",output_format="BAM",
output_file=targets$OutputFile,unique=TRUE,indels=5)
###################################################
### code chunk number 3: Rsubread.Rnw:96-100
###################################################
reads1 <- system.file("extdata","reads1.txt.gz",package="Rsubread")
reads2 <- system.file("extdata","reads2.txt.gz",package="Rsubread")
align(index="reference_index",readfile1=reads1,readfile2=reads2,
output_file="alignResultsPE.BAM",phredOffset=64)

fc <- featureCounts(files=targets$OutputFile,annot.inbuilt="hg19")
fc$counts[1:5,]
fc$annotation[1:5,c("GeneID","Length")]
x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])
x_rpkm <- rpkm(x,x$genes$Length,prior.count=0)
isexpr <- rowSums(cpm(x) > 10) >= 2
x <- x[isexpr,]

# creating a Design Matrix
celltype <- factor(targets$CellType)
design <- model.matrix(~0+celltype)
colnames(design) <- levels(celltype)

#Normalization
y <- voom(x,design,plot=TRUE)
# sample clustering
plotMDS(y,xlim=c(-2.5,2.5))
# linear modeling
fit <- lmFit(y,design)
contr <- makeContrasts(BvsA=B-A,levels=design)
fit.contr <- eBayes(contrasts.fit(fit,contr))
dt <- decideTests(fit.contr)
summary(dt)
# top ten diff genes
options(digits=2)
topTable(fit.contr)

###################################################
### code chunk number 4: Rsubread.Rnw:121-131
###################################################
ann <- data.frame(
GeneID=c("gene1","gene1","gene2","gene2"),
Chr="chr_dummy",
Start=c(100,1000,3000,5000),
End=c(500,1800,4000,5500),
Strand=c("+","+","-","-"),
stringsAsFactors=FALSE)
ann
fc_SE <- featureCounts("alignResults.BAM",annot.ext=ann)
fc_SE


###################################################
### code chunk number 5: Rsubread.Rnw:138-140
###################################################
fc_PE <- featureCounts("alignResultsPE.BAM",annot.ext=ann,isPairedEnd=TRUE)
fc_PE


###################################################
### code chunk number 6: Rsubread.Rnw:160-162
###################################################
x <- qualityScores(filename=reads,offset=64,nreads=1000)
x[1:10,1:10]


###################################################
### code chunk number 7: Rsubread.Rnw:187-188
###################################################
propmapped("alignResults.BAM")


