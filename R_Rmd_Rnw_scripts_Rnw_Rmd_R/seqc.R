### R code from vignette source 'seqc.Rnw'
# To run this case study, you should have R version of 3.0.2 or later

# load libraries
library(Rsubread)
library(limma)
library(edgeR)

# read in target file
options(digits=2)
targets <- readTargets()
targets

# create a design matrix
celltype <- factor(targets$CellType)
design <- model.matrix(~celltype)

# build an index for reference sequence (Chr1 in hg19)
buildindex(basename="chr1",reference="hg19_chr1.fa")

# align reads
align(index="chr1",readfile1=targets$InputFile,input_format="gzFASTQ",output_format="BAM",output_file=targets$OutputFile,unique=TRUE,indels=5)

# count numbers of reads mapped to NCBI Refseq genes
fc <- featureCounts(files=targets$OutputFile,annot.inbuilt="hg19")
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


###################################################
### code chunk number 2: seqc.Rnw:88-89
###################################################
ROC_aceview_gene_MGP[1:15,]


###################################################
### code chunk number 3: seqc.Rnw:91-93
###################################################
colnames(ILM_aceview_gene_BGI)
ILM_aceview_gene_BGI[1:15,1:7]


###################################################
### code chunk number 4: seqc.Rnw:103-104
###################################################
ILM_junction_AGR_B[1:15,]


###################################################
### code chunk number 5: seqc.Rnw:114-116
###################################################
colnames(taqman)
taqman[1:15, 1:9]

# NOTE: Only given for reference
#       Not supposed to be executed during the lab
samtools view accepted_hits_137_1.bam | \
 sort > accepted_hits_prehtseq_137_1.sam
htseq-count -s no -q accepted_hits_prehtseq_137_1.sam \
 Homo_sapiens.GRCh37.71.gtf > 137_1.counts

This was run for each of the samples and the counts were combined into a single table. You can get the count table from the data directory. You can run R on UPPMAX, or download the file to your local computer and do the analysis locally if you prefer.

The code to run in R:

# read in the data
counts <- read.delim("count_table.txt")
head(counts)

As you can see, the samples are ordered with the 3 replicates from each group next to each other. So when we are to define colors for the samples we only have to repeat each color 3 times (this may not always be the case!)

# define colors:
col.def<-c("red","blue","green","magenta")
sample.def<-c("ctrl", "t2h", "t6h", "t24h")
colors <- rep(col.def, each=3)

Start with a PCA to se the general distribution. PCA of RNA-seq data is usually performed in log-scale. We also add a pseudo-count of +1 to avoid logging zero (gives infinity). You need to transpose - t() - the data matrix, otherwise you will run PCA on the genes instead of samples.

myPca <- prcomp(t(log2(counts+1)))

This creates a list that contains:

    the samples mapping to each PC in myPca$x
    PC contribution to variance in myPca$sdev
    PC loadings for each gene in myPca$rotation

Now some plotting. In R you can either plot into a default window or direct all your output to a “device”, that can be pdf, png, tiff etc. To open a new pdf device:

pdf('pca_plot.pdf')
# once you have plotted all you want to put into that file,
# close it with dev.off()

Let’s first make a simple plot of the first two principal components (PC1 vs PC2):

plot(myPca$x[,1],myPca$x[,2],col=colors,pch=1)
legend("topright",sample.def,pch=1,col=col.def)
dev.off()

Sometimes the first two PCs may not be the ones that will best separate the sample groups, so it is a good idea to look at more PCs. Here is one example that shows how to plot the top 5 PCs:

pdf('pca*plot*5pc.pdf')
tmpPcaData <- as.data.frame(myPca$x[,1:5])
plot(tmpPcaData, col=colors,pch=1)
dev.off()

Another thing to look at is the pairwise correlation between all the samples and see how they group based on correlation. Let’s create one matrix with all pairwise Pearson correlations (again in log-space).

nSamples<-ncol(counts)
C<-cor(log2(counts+1),method="pearson")

Now create a hierarchical clustering dendrogram using 1-correlation as the distance measure:

	d <- as.dist(1-C)
h <- hclust(d,method="ward.D2")
dendro <- as.dendrogram(h)

Now you will plot a heatmap with the correlations:

pdf('correlation_heatmap.pdf')
heatmap(C,Colv=dendro,Rowv=dendro)
dev.off()


