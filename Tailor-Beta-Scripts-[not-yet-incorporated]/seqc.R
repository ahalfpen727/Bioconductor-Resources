### R code from vignette source 'seqc.Rnw'
# To run this case study, you should have R version of 3.0.2 or later
library(cummeRbund)

# load libraries
library(Rsubread)
library(limma)
library(edgeR)
library(seqc)
buildindex(basename="grch38_indexed_genome",reference="genome.fa")

cuff<-readCufflinks(gtfFile="cuffcmp.combined.gtf",genome="genome.fa",rebuild=F)
 refgtf="cuffcmp.combined.gtf"
 genome_path="genome.fa"
 over="LUTS"
 under="CTRL"
# create a design matrix
# # read in target file
options(digits=2)
targets <- readTargets()
targets
celltype <- factor(targets$CellType)
design <- model.matrix(~celltype)

# build an index for reference sequence (Chr1 in hg19)
buildindex(basename="chr1",reference="chr1.fa")
buildindex(basename="chr2",reference="chr2.fa")
buildindex(basename="chr3",reference="chr3.fa")
buildindex(basename="chr4",reference="chr4.fa")
buildindex(basename="chr5",reference="chr5.fa")
buildindex(basename="chr6",reference="chr6.fa")
buildindex(basename="chr7",reference="chr7.fa")
buildindex(basename="chr8",reference="chr8.fa")
buildindex(basename="chr9",reference="chr9.fa")
buildindex(basename="chr10",reference="chr10.fa")
buildindex(basename="chr11",reference="chr11.fa")
buildindex(basename="chr12",reference="chr12.fa")
buildindex(basename="chr13",reference="chr13.fa")
buildindex(basename="chr14",reference="chr14.fa")
buildindex(basename="chr15",reference="chr15.fa")
buildindex(basename="chr16",reference="chr16.fa")
buildindex(basename="chr17",reference="chr17.fa")
buildindex(basename="chr18",reference="chr18.fa")
buildindex(basename="chr19",reference="chr19.fa")
buildindex(basename="chr20",reference="chr20.fa")
buildindex(basename="chr21",reference="chr21.fa")
buildindex(basename="chr22",reference="chr22.fa")
buildindex(basename="chr23",reference="chr23.fa")
buildindex(basename="chrM",reference="chrM.fa")
buildindex(basename="chrX",reference="chrX.fa")
buildindex(basename="chrY",reference="chrYfa")

# align reads
align(index="chr1",readfile1=targets$InputFile,input_format="gzFASTQ",output_format="BAM",output_file=targets$OutputFile,unique=TRUE,indels=5)

# count numbers of reads mapped to NCBI Refseq genes
# cuffmerge produ
fc <- featureCounts(files="GRCh38_unguidedMERGE.sam", 
					annot.ext="cuffcmp.combined.gtf",
					isGTFAnnotationFile=T,GTF.attrType="exons")

x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])

# generate RPKM values if you need them
x_rpkm <- rpkm(x,x$genes$Length)

# filter out low-count genes
isexpr <- rowSums(cpm(x) > 10) >= 2
x <- x[isexpr,]

# perform voom normalization
y <- voom(x,design,plot=TRUE)

# cluster libraries
plotMDS(y,xlim=c(-2.5,2.5))

# fit linear model and assess differential expression
fit <- eBayes(lmFit(y,design))
topTable(fit,coef=2)

###################################################
### code chunk number 1: seqc.Rnw:31-34
###################################################
library(seqc)
options(width=110, digits=2)
ls(2)
ROC_aceview_gene_MGP[1:15,]
colnames(ILM_aceview_gene_BGI)
ILM_aceview_gene_BGI[1:15,1:7]
ILM_junction_AGR_B[1:15,]
colnames(taqman)
taqman[1:15, 1:9]
###################################################
# NOTE: Only given for reference
#       Not supposed to be executed during the lab
#samtools view accepted_hits_137_1.bam | \
# sort > accepted_hits_prehtseq_137_1.sam
#htseq-count -s no -q accepted_hits_prehtseq_137_1.sam \
# Homo_sapiens.GRCh37.71.gtf > 137_1.counts

#The code to run in R:

# read in the data
counts <- read.delim("count_table.txt")
head(counts)

# define colors:
col.def<-c("red","blue","green","magenta")
sample.def<-c("ctrl", "t2h", "t6h", "t24h")
colors <- rep(col.def, each=3)

#PCA

myPca <- prcomp(t(log2(counts+1)))

#This creates a list that contains:

  # the samples mapping to each PC in myPca$x
  #  PC contribution to variance in myPca$sdev
  #  PC loadings for each gene in myPca$rotation
#  plot of the first two principal components (PC1 vs PC2):
plot(myPca$x[,1],myPca$x[,2],col=colors,pch=1)
legend("topright",sample.def,pch=1,col=col.def)
dev.off()

pdf('pca*plot*5pc.pdf')
tmpPcaData <- as.data.frame(myPca$x[,1:5])
plot(tmpPcaData, col=colors,pch=1)
dev.off()

#Another thing to look at is the pairwise correlation between all the samples
#grouping based on correlation.
#Let’s create one matrix with all pairwise Pearson correlations (again in log-space).

nSamples<-ncol(counts)
C<-cor(log2(counts+1),method="pearson")

# Now create a hierarchical clustering dendrogram using 1-correlation as the distance measure:

	d <- as.dist(1-C)
h <- hclust(d,method="ward.D2")
dendro <- as.dendrogram(h)

#Now you will plot a heatmap with the correlations:

pdf('correlation_heatmap.pdf')
heatmap(C,Colv=dendro,Rowv=dendro)
dev.off()


