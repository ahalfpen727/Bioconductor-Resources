source("https://bioconductor.org/biocLite.R")
biocLite(c("Rsubread","limma","edgeR","cummeRbund","Org.Hs.eg.db"))
library(readxl)
library(Rsubread)
library(limma)
library(edgeR)
A_over_1021 <- read_excel("A-over-1021.xlsx")
A_over_1021$GeneName
#############################################################################
# 20 samples

plotMDSTopGenes <- function (x, top = 500, labels = NULL, pch = NULL, cex = 1, dim.plot = c(1, 2), ndim = max(dim.plot), gene.selection = "pairwise", xlab = NULL, ylab = NULL, ...)
{
    x <- as.matrix(x)
    nsamples <- ncol(x)
    if (nsamples < 3) 
        stop(paste("Only", nsamples, "columns of data: need at least 3"))
    cn <- colnames(x)
    bad <- rowSums(is.finite(x)) < nsamples
    if (any(bad)) 
        x <- x[!bad, , drop = FALSE]
    nprobes <- nrow(x)
    top <- min(top, nprobes)
    if (is.null(pch) & is.null(labels)) {
        labels <- colnames(x)
        if (is.null(labels)) 
            labels <- 1:nsamples
    }
    if (!is.null(labels)) 
        labels <- as.character(labels)
    dim.plot <- unique(as.integer(dim.plot))
    if (length(dim.plot) != 2L) 
        stop("dim.plot must specify two dimensions to plot")
    if (ndim < 2L) 
        stop("Need at least two dim.plot")
    if (nsamples < ndim) 
        stop("ndim is greater than number of samples")
    if (nprobes < ndim) 
        stop("ndim is greater than number of rows of data")
    gene.selection <- match.arg(gene.selection, c("pairwise", 
                                                  "common"))
    dd <- matrix(0, nrow = nsamples, ncol = nsamples, dimnames = list(cn, 
                                                          cn))
    if (gene.selection == "pairwise") {
        topindex <- nprobes - top + 1L
        for (i in 2L:(nsamples)) for (j in 1L:(i - 1L)) dd[i, 
                                                           j] = sqrt(mean(sort.int((x[, i] - x[, j])^2, partial = topindex)[topindex:nprobes]))
        axislabel <- "Leading logFC dim"
    }
    else {
        if (nprobes > top) {
            s <- rowMeans((x - rowMeans(x))^2)
            o <- order(s, decreasing = TRUE)
            x <- x[o[1L:top], , drop = FALSE]
        }
        for (i in 2L:(nsamples)) {
            dd[i, 1L:(i - 1L)] = sqrt(colMeans((x[, i] - x[, 1:(i - 1), drop = FALSE])^2))
        }
        axislabel <- "Principal Component"
    }
    a1 <- suppressWarnings(cmdscale(as.dist(dd), k = ndim))
    mds <- new("MDS", list(dim.plot = dim.plot,
                           distance.matrix = dd,
                           cmdscale.out = a1,
                           top = top,
                           gene.selection = gene.selection))
    if (dim.plot[1] > ncol(a1)) {
        mds$x <- rep.int(0, nsamples)
        warning(paste("dimension", dim.plot[1], "is degenerate or all zero"))
    }
    else mds$x <- a1[, dim.plot[1]]
    if (dim.plot[2] > ncol(a1)) {
        mds$y <- rep.int(0, nsamples)
        warning(paste("dimension", dim.plot[2], "is degenerate or all zero"))
    }
    else mds$y <- a1[, dim.plot[2]]
    mds$top <- top
    mds$topGenes <- x
    mds$distances <- s[o[1L:top]]
    mds$axislabel <- axislabel
    plotMDS(mds, labels = labels, pch = pch, cex = cex, xlab = xlab, 
            ylab = ylab, ...)
}

## plotMDS.default <-
##     function (x, top = 500, labels = colnames(x), col = NULL, cex = 1,
##               dim.plot = c(1, 2), ndim = max(dim.plot), gene.selection = "pairwise",
##               xlab = paste("Dimension", dim.plot[1]), ylab = paste("Dimension", dim.plot[2]))
## {
##     x <- as.matrix(x)
##     ok <- is.finite(x)
##     if (!all(ok))
##         x <- x[apply(ok, 1, all), ]
##     if (is.null(labels))
##         labels <- 1:dim(x)[2]
##     nprobes <- nrow(x)
##     nsamples <- ncol(x)
##     if (ndim < 2)
##         stop("Need at least two dim.plot")
##     if (nsamples < ndim)
##         stop("Two few samples")
##     gene.selection <- match.arg(gene.selection, c("pairwise",
##                                                   "common"))
##     cn <- colnames(x)
##     dd <- matrix(0, nrow = nsamples, ncol = nsamples, dimnames = list(cn,
##                                                           cn))
##     topindex <- nprobes - top + 1
##     if (gene.selection == "pairwise") {
##         for (i in 2:(nsamples))
##             for (j in 1:(i - 1))
##                 dd[i, j] = sqrt(meansort.int((x[,i] - x[, j])^2, partial = topindex)[topindex:nprobes]))
## }
## else {

##                                         # Same genes used for all comparisons ,"common"
##     s <- rowMeans((x - rowMeans(x))^2)
##     q <- quantile(s, p = (topindex - 1.5)/(nprobes - 1))

##     x <- x[s >= q, ]
##                                         # an extra line
##     ind.top.genes<-which(s >= q)

##     for (i in 2:(nsamples))
##         dd[i, 1:(i - 1)] = sqrt(colMeans((x[,i] - x[, 1:(i - 1), drop = FALSE])^2))
## }
## a1 <- cmdscale(as.dist(dd), k = ndim)
## mds <- new("MDS", list(dim.plot = dim.plot, distance.matrix = dd,
##                        cmdscale.out = a1, top = top, gene.selection = gene.selection))
## mds$x <- a1[, dim.plot[1]]
## mds$y <- a1[, dim.plot[2]]
## mdsPlot<-plotMDS(mds, labels = labels, col = col, cex = cex, xlab = xlab,
##                  ylab = ylab, ...)
## list (mds=mds, ind.top.genes=ind.top.genes)
## }


##                                         # The example from "plotMDS" documentation

## sd <- 0.3*sqrt(4/rchisq(1000,df=4))
## x <- matrix(rnorm(1000*6,sd=sd),1000,6)
## rownames(x) <- paste("Gene",1:1000)
## x[1:50,4:6] <- x[1:50,4:6] + 2
                                        # without labels, indexes of samples are plotted.


fc = featureCounts(

    files = "./createLinks_results_default/Sample_10211_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_10212_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_10213_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_10214_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_10215_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_1021B1_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_1021B2_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_1021B3_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_1021B4_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_1021B5_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_A1_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_A2_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_A3_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_A4_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_A5_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_AB1_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_AB2_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_AB3_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_AB4_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_AB5_R1.fastq.gz.subread.BAM",

    ##                                     # BAM/SAM files
    ## files = "./createLinks_results_default/Sample_10212_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_10214_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_10215_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_1021B2_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_1021B3_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_1021B4_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_1021B5_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_A2_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_A4_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_A5_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_AB1_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_AB2_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_AB3_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_AB5_R1.fastq.gz.subread.BAM",

    ## files = "./tophat_results_Smeliloti_1021_unguided/Sample_10212_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_10214_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_10215_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_1021B2_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_1021B3_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_1021B4_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_1021B5_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_A2_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_A4_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_A5_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_AB1_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_AB2_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_AB3_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_AB5_out/accepted_hits.bam",

    ## files = c(
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_10212_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_10214_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_10215_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_1021B2_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_1021B3_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_1021B4_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_1021B5_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_A2_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_A4_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_A5_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_AB1_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_AB2_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_AB3_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_AB5_out/accepted_hits.bam
    ## ),
                                        # annotation
    annot.inbuilt=NULL,
    annot.ext="Smeliloti_1021_refseq_bowtie_assembly/1021_genome.gtf",
    isGTFAnnotationFile=TRUE,
    GTF.featureType="CDS",
    GTF.attrType="Parent",
    chrAliases=NULL,
                                        # level of summarization
    useMetaFeatures=TRUE,
                                        # overlap between reads and features
    allowMultiOverlap=FALSE,
    minOverlap=1,
    ## fracOverlap=0,
    largestOverlap=FALSE,
    readExtension5=0,
    readExtension3=0,
    read2pos=NULL,
                                        # multi-mapping reads
    countMultiMappingReads=FALSE,
                                        # fractional counting
    fraction=FALSE,
                                        # read filtering
    minMQS=0,
    splitOnly=FALSE,
    nonSplitOnly=FALSE,
    primaryOnly=FALSE,
    ignoreDup=FALSE,
                                        # strandness
    strandSpecific=0,
                                        # exon-exon junctions
    juncCounts=FALSE,
    genome=NULL,
                                        # parameters specific to paired end reads
    isPairedEnd=TRUE,
    requireBothEndsMapped=TRUE,
    checkFragLength=FALSE,
    minFragLength=50,
    maxFragLength=1600,
    countChimericFragments=TRUE,
    autosort=TRUE,
                                        # number of CPU threads
    nthreads=1,
                                        # miscellaneous
    maxMOp=10,
    ## tmpDir=".",
    reportReads=FALSE

)

## saveRDS(fc, "featureCounts.list")
fc = readRDS("featureCounts.list")
rownames(fc$counts)
colnames(fc$counts)
rownames(fc$annotation)<-fc$annotation$GeneID
colnames(fc$annotation)
rownames(fc$annotation)


GeneSpecs<-fc$annotation[,c("GeneID","Length")]
GeneAnnotations<-A_over_1021[,1:2]
head(GeneAnnotations)

GeneIDandName<-rownames(fc$annotation) %in% GeneAnnotations$GeneID
fc$annotation<-fc$annotation[GeneIDandName,]

GeneIDandName<-rownames(fc$counts)%in% GeneAnnotations$GeneID
fc$counts<-fc$counts[GeneIDandName,]
fc$annotation<-merge(GeneAnnotations,fc$annotation, by="GeneID")

## write.table(fc$annotation, file="fc.annotation.txt", sep="\t", row.names=FALSE, quote=FALSE)
fc$annotation = read.table(file="fc.annotation.names.txt", sep="\t", header=TRUE)

## Create a DGEList object.
x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneName","GeneID","Length")])
x <- DGEList(counts=fc$counts, genes=fc$annotation)
x
#Calculate RPKM (reads per kilobases of exon per million reads mapped) values for genes:
x_rpkm <- rpkm(x,x$genes$Length,prior.count=0)
x_rpkm[1:5,]
## A_1.bam A_2.bam B_1.bam B_2.bam
## 653635        939   905.0     709     736
## 47
## 100422834      19     0.0       0       0
## 645520         11     8.1       0       0
## 79501           0     0.0       0       0
## 729737         62    64.9      19      16

## Filtering. Only keep in the analysis those genes which had
## 10 reads per million mapped reads in at least two libraries.
## isexpr <- rowSums(cpm(x) > 10) >= 2
isexpr <- rowSums(cpm(x) > 10) >= 3
x <- x[isexpr,]

## Create a design matrix:
fc$group = c('wt1021','wt1021','wt1021','wt1021','wt1021','wt1021B','wt1021B','wt1021B','wt1021B','wt1021B','A','A','A','A','A','AB','AB','AB','AB','AB')

groupsFactor <- factor(fc$group)
design <- model.matrix(~0+groupsFactor)
colnames(design) <- levels(groupsFactor)
## celltype <- factor(targets$CellType)
## design <- model.matrix(~0+celltype)
## colnames(design) <- levels(celltype)

## Normalization. Perform voom normalization:
## y <- voom(x,design,plot=FALSE)
y <- voom(x,design,plot=TRUE)

## Multi-dimensional scaling (MDS) plot
num=5000
num=1000
mds = plotMDSTopGenes(y,
        ## xlim=c(-2.5,2.5),
        ## xlim=c(-3.5,3.5),
        labels=c('1021','1021','1021','1021','1021','1021B','1021B','1021B','1021B','1021B','A','A','A','A','A','AB','AB','AB','AB','AB'),
        pch=21,
        gene.selection='common',
        top=num,
        main=paste0('MDS Plot using top ',num,' pairwise distances'),
        col=c('blue', 'blue', 'blue', 'blue', 'blue',
            'red', 'red', 'red', 'red', 'red',
            'green', 'green', 'green', 'green', 'green',
            'magenta', 'magenta', 'magenta', 'magenta', 'magenta')
        )

 
length(mds$topGenes)
dim(mds$topGenes)
rownames(mds$topGenes)
topGenesWithDistances = data.frame(Genes=rownames(mds$topGenes), Distances=mds$distances)
## write.table(topGenesWithDistances, file="top25Genes.txt", sep="\t", col.names=FALSE, row.names=TRUE, quote=FALSE)
write.table(topGenesWithDistances, file="top1000Genes.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

quartz()
plot(1:150,
     main="Elbow plot of mean Euclidean distances per gene for 1021 vs 1021B A vs AB",
     xlab="Gene",
     ylab="Mean Euclidean Distance",
     ## cex.lab=2,
     cex.main=2,
     topGenesWithDistances$Distances[1:150],
     type="o",
     lty=3,
     pch=15)

############################################################################################################

design.A.AB = design[11:20,1:2]
design.A.AB
x.A.AB = x[,11:20]
dim(x.A.AB)

y.A.AB <- voom(x.A.AB,design.A.AB,plot=FALSE)

## Multi-dimensional scaling (MDS) plot
num=1000
mds = plotMDSTopGenes(y.A.AB,
        ## xlim=c(-2.5,2.5),
        ## xlim=c(-3.5,3.5),
        labels=c('A','A','A','A','A','AB','AB','AB','AB','AB'),
#        labels=c('A','AB'),
        pch=21,
        gene.selection='common',
        top=num,
        main=paste0('MDS Plot using top ',num,' pairwise distances'),
        col=c('green', 'green', 'green', 'green', 'green',
            'magenta', 'magenta', 'magenta', 'magenta', 'magenta')
        )
length(mds$topGenes)
dim(mds$topGenes)
rownames(mds$topGenes)
topGenesWithDistances = data.frame(Genes=rownames(mds$topGenes), Distances=mds$distances)
## write.table(topGenesWithDistances, file="top25Genes.txt", sep="\t", col.names=FALSE, row.names=TRUE, quote=FALSE)
write.table(topGenesWithDistances, file="top1000Genes.A.AB.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

quartz()
plot(1:100,
     main="Elbow plot of mean Euclidean distances per gene for A vs AB",
     xlab="Gene",
     ylab="Mean Euclidean Distance",
     ## cex.lab=2,
     cex.main=2,
     topGenesWithDistances$Distances[1:100],
     type="o",
     lty=3,
     pch=15)

############################################################################################################
############################################################################################################
design
design.A.wt1021 = design[c(1:5,11:15),c(1,3)]
design.A.wt1021
x.A.wt1021 = x[,c(1:5,11:15)]
dim(x.A.wt1021)

y.A.wt1021 <- voom(x.A.wt1021,design.A.wt1021,plot=FALSE)

## Multi-dimensional scaling (MDS) plot
num=1000
mds = plotMDSTopGenes(y.A.wt1021,
                      ## xlim=c(-2.5,2.5),
                      ## xlim=c(-3.5,3.5),
                      labels=c('A','A','A','A','A','wt1021','wt1021','wt1021','wt1021','wt1021'),
                      pch=21,
                      gene.selection='common',
                      top=num,
                      main=paste0('MDS Plot using top ',num,' pairwise distances'),
                      col=c('green', 'green', 'green', 'green', 'green',
                            'magenta', 'magenta', 'magenta', 'magenta', 'magenta')
)
length(mds$topGenes)
dim(mds$topGenes)
rownames(mds$topGenes)
topGenesWithDistances = data.frame(Genes=rownames(mds$topGenes), Distances=mds$distances)
## write.table(topGenesWithDistances, file="top25Genes.txt", sep="\t", col.names=FALSE, row.names=TRUE, quote=FALSE)
write.table(topGenesWithDistances, file="top1000Genes.A.wt1021.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

quartz()
plot(1:100,
     main="Elbow plot of mean Euclidean distances per gene for A vs 1021",
     xlab="Gene",
     ylab="Mean Euclidean Distance",
     ## cex.lab=2,
     cex.main=2,
     topGenesWithDistances$Distances[1:100],
     type="o",
     lty=3,
     pch=15)

############################################################################################################

design.1021.1021B = design[1:10,3:4]
design.1021.1021B
x.1021.1021B = x[,1:10]
dim(x.1021.1021B)

y.1021.1021B <- voom(x.1021.1021B,design.1021.1021B,plot=FALSE)

## Multi-dimensional scaling (MDS) plot
num=1000
mds = plotMDSTopGenes(y.1021.1021B,
        ## xlim=c(-2.5,2.5),
        ## xlim=c(-3.5,3.5),
        labels=c('1021','1021','1021','1021','1021','1021B','1021B','1021B','1021B','1021B'),
        pch=21,
        gene.selection='common',
        top=num,
        main=paste0('MDS Plot using top ',num,' pairwise distances'),
        col=c('blue', 'blue', 'blue', 'blue', 'blue',
            'red', 'red', 'red', 'red', 'red')
        )

length(mds$topGenes)
dim(mds$topGenes)
rownames(mds$topGenes)
topGenesWithDistances = data.frame(Genes=rownames(mds$topGenes), Distances=mds$distances)
## write.table(topGenesWithDistances, file="top25Genes.txt", sep="\t", col.names=FALSE, row.names=TRUE, quote=FALSE)
write.table(topGenesWithDistances, file="top1000Genes.1021.1021B.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

quartz()
plot(1:100,
     main="Elbow plot of mean Euclidean distances per gene for 1021 vs 1021B",
     xlab="Gene",
     ylab="Mean Euclidean Distance",
     ## cex.lab=2,
     cex.main=2,
     topGenesWithDistances$Distances[1:100],
     type="o",
     lty=3,
     pch=15)


## Linear model fitting and differential expression analysis.
## Fit linear models to genes and assess differential expression using eBayes moderated t statistic.  Here we compare sample
## B vs sample A.
fit <- lmFit(y,design)
options(digits=3)








contr <- makeContrasts(wt1021Bvswt1021=wt1021B - wt1021,levels=design)
# write.fit(fit.contr, file="wt1021B-over-wt1021.txt")
fit.contr <- eBayes(contrasts.fit(fit,contr))
dt <- decideTests(fit.contr)
summary(dt)

## List top 10 differentially expressed genes:
topTable(fit.contr)

## topTable(fit.contr, sort="none",n=Inf)
topOutput = topTable(fit.contr, n=Inf)
topOutput$FC = 2^(topOutput$logFC)
write.table(topOutput, file="wt1021B-over-wt1021.txt", sep="\t", quote=FALSE)








contr <- makeContrasts(Avswt1021=A-wt1021,levels=design)
# write.fit(fit.contr, file="A-over-wt1021.txt")
fit.contr <- eBayes(contrasts.fit(fit,contr))
dt <- decideTests(fit.contr)
summary(dt)

## List top 10 differentially expressed genes:
topTable(fit.contr)

## topTable(fit.contr, sort="none",n=Inf)
topOutput = topTable(fit.contr, n=Inf)
topOutput$FC = 2^(topOutput$logFC)
write.table(topOutput, file="A-over-wt1021.txt", sep="\t", quote=FALSE)









contr <- makeContrasts(ABvswt1021=AB-wt1021,levels=design)
# write.fit(fit.contr, file="AB-over-wt1021.txt")
fit.contr <- eBayes(contrasts.fit(fit,contr))
dt <- decideTests(fit.contr)
summary(dt)

## List top 10 differentially expressed genes:
topTable(fit.contr)

## topTable(fit.contr, sort="none",n=Inf)
topOutput = topTable(fit.contr, n=Inf)
topOutput$FC = 2^(topOutput$logFC)
write.table(topOutput, file="AB-over-wt1021.txt", sep="\t", quote=FALSE)









contr <- makeContrasts(AvsAB=A-AB,levels=design)
# write.fit(fit.contr, file="A-over-AB.txt")
fit.contr <- eBayes(contrasts.fit(fit,contr))
dt <- decideTests(fit.contr)
summary(dt)

## List top 10 differentially expressed genes:
topTable(fit.contr)

## topTable(fit.contr, sort="none",n=Inf)
topOutput = topTable(fit.contr, n=Inf)
topOutput$FC = 2^(topOutput$logFC)
write.table(topOutput, file="A-over-AB.txt", sep="\t", quote=FALSE)








contr <- makeContrasts(ABvswt1021B=AB-wt1021B,levels=design)
# write.fit(fit.contr, file="AB-over-wt1021B.txt")
fit.contr <- eBayes(contrasts.fit(fit,contr))
dt <- decideTests(fit.contr)
summary(dt)

## List top 10 differentially expressed genes:
topTable(fit.contr)

## topTable(fit.contr, sort="none",n=Inf)
topOutput = topTable(fit.contr, n=Inf)
topOutput$FC = 2^(topOutput$logFC)
write.table(topOutput, file="AB-over-wt1021B.txt", sep="\t", quote=FALSE)









## GeneID Length logFC AveExpr   t P.Value adj.P.Val  B
## 100131754 100131754   1019   1.6      16 113 3.5e-28   6.3e-25 54
## 2023           2023   1812  -2.7      13 -91 2.2e-26   1.9e-23 51
## 2752           2752   4950   2.4      13  82 1.5e-25   9.1e-23 49
## 22883         22883   5192   2.3      12  64 1.8e-23   7.9e-21 44
## 6135           6135    609  -2.2      12 -62 3.1e-23   9.5e-21 44
## 6202           6202    705  -2.4      12 -62 3.2e-23   9.5e-21 44
