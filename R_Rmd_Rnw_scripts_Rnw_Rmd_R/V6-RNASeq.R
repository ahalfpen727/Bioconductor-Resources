## ----style, echo = FALSE, results = 'asis'--------------------------------------------------------
BiocStyle::markdown()
options(width=100, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")))

## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE--------------------------------------------
suppressPackageStartupMessages({
    library(DESeq2)
    library(limma)
    library(airway)
    library(gplots)
    library(RColorBrewer)
    library(ggplot2)
    library(genefilter)
    library(org.Hs.eg.db)
})

## ----configure-test-------------------------------------------------------------------------------
stopifnot(
    getRversion() >= '3.2' && getRversion() < '3.3',
    BiocInstaller::biocVersion() == "3.2"
)

## -------------------------------------------------------------------------------------------------
library(airway)
data("airway")
se <- airway

## -------------------------------------------------------------------------------------------------
head(assay(se))

## -------------------------------------------------------------------------------------------------
colSums(assay(se))

## -------------------------------------------------------------------------------------------------
colData(se)

## ----rowRanges`-----------------------------------------------------------------------------------
rowRanges(se)

## -------------------------------------------------------------------------------------------------
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ cell + dex)

## -------------------------------------------------------------------------------------------------
rld <- rlog(dds)
head(assay(rld))

## ----rldplot, fig.width=10, fig.height=5----------------------------------------------------------
opar <- par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)
plot( log2( 1 + counts(dds, normalized=TRUE)[ , 1:2] ),
     col=rgb(0,0,0,.2), pch=16, cex=0.3 )
plot( assay(rld)[ , 1:2],
     col=rgb(0,0,0,.2), pch=16, cex=0.3 )
par(opar)

## -------------------------------------------------------------------------------------------------
sampleDists <- dist( t( assay(rld) ) )
sampleDists

## -------------------------------------------------------------------------------------------------
library("gplots")
library("RColorBrewer")

## ----distheatmap, fig.width=8---------------------------------------------------------------------
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hc <- hclust(sampleDists)
heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none", col=colors,
          margins=c(2,10), labCol=FALSE )

## ----plotpca, fig.width=6, fig.height=4.5---------------------------------------------------------
plotPCA(rld, intgroup = c("dex", "cell"))

## -------------------------------------------------------------------------------------------------
(data <- plotPCA(rld, intgroup = c( "dex", "cell"), returnData=TRUE))
percentVar <- round(100 * attr(data, "percentVar"))

## -------------------------------------------------------------------------------------------------
library("ggplot2")

## ----ggplotpca, fig.width=6, fig.height=4.5-------------------------------------------------------
qplot(PC1, PC2, color=dex, shape=cell, data=data) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance"))

## -------------------------------------------------------------------------------------------------
dds$dex <- relevel(dds$dex, "untrt")

## -------------------------------------------------------------------------------------------------
dds <- DESeq(dds)

## -------------------------------------------------------------------------------------------------
(res <- results(dds))

## -------------------------------------------------------------------------------------------------
mcols(res, use.names=TRUE)

## -------------------------------------------------------------------------------------------------
summary(res)

## -------------------------------------------------------------------------------------------------
results(dds, contrast=c("cell", "N061011", "N61311"))

## -------------------------------------------------------------------------------------------------
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))

## -------------------------------------------------------------------------------------------------
sum(res$padj < 0.1, na.rm=TRUE)

## -------------------------------------------------------------------------------------------------
resSig <- subset(res, padj < 0.1)
head(resSig[ order( resSig$log2FoldChange ), ])

## -------------------------------------------------------------------------------------------------
head(resSig[ order( -resSig$log2FoldChange ), ])

## ----plotcounts, fig.width=5, fig.height=5--------------------------------------------------------
topGene <- rownames(res)[which.min(res$padj)]
data <- plotCounts(dds, gene=topGene, intgroup=c("dex"), returnData=TRUE)

## ----ggplotcountsdot, fig.height=5----------------------------------------------------------------
ggplot(data, aes(x=dex, y=count, fill=dex)) +
  scale_y_log10() + 
  geom_dotplot(binaxis="y", stackdir="center")

## ----plotma, eval=FALSE---------------------------------------------------------------------------
#  plotMA(res, ylim=c(-5,5))

## ----plotma2, eval=FALSE--------------------------------------------------------------------------
#  plotMA(res, ylim=c(-5,5))
#  with(res[topGene, ], {
#    points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
#    text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
#  })

## ----plotdispests---------------------------------------------------------------------------------
plotDispEsts(dds)

## ----histpvalue-----------------------------------------------------------------------------------
hist(res$pvalue, breaks=20, col="grey50", border="white")

## ----histpvalue2----------------------------------------------------------------------------------
hist(res$pvalue[res$baseMean > 1], breaks=20, col="grey50", border="white")

## -------------------------------------------------------------------------------------------------
library("genefilter")
topVarGenes <- head(order(-rowVars(assay(rld))),35)

## ----genescluster, fig.height=9-------------------------------------------------------------------
colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ rld$dex ]
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
colnames(mat) <- paste0(rld$dex,"-",rld$cell)
heatmap.2(mat, trace="none", col=colors, ColSideColors=sidecols,
          labRow=FALSE, mar=c(10,2), scale="row")

## ----sensitivityovermean, fig.height=4------------------------------------------------------------
# create bins using the quantile function
qs <- c(0, quantile(res$baseMean[res$baseMean > 0], 0:7/7))
# cut the genes into the bins
bins <- cut(res$baseMean, qs)
# rename the levels of the bins using the middle point
levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
# calculate the ratio of $p$ values less than .01 for each bin
ratios <- tapply(res$pvalue, bins, function(p) mean(p < .01, na.rm=TRUE))
# plot these ratios
barplot(ratios, xlab="mean normalized count", ylab="ratio of small p values")

## -------------------------------------------------------------------------------------------------
library(org.Hs.eg.db)

## -------------------------------------------------------------------------------------------------
columns(org.Hs.eg.db)
res$hgnc_symbol <- 
    unname(mapIds(org.Hs.eg.db, rownames(res), "SYMBOL", "ENSEMBL"))
res$entrezgene <- 
    unname(mapIds(org.Hs.eg.db, rownames(res), "ENTREZID", "ENSEMBL"))

## -------------------------------------------------------------------------------------------------
resOrdered <- res[order(res$pvalue),]
head(resOrdered)

## ----eval=FALSE-----------------------------------------------------------------------------------
#  write.csv(as.data.frame(resOrdered), file="results.csv")

## -------------------------------------------------------------------------------------------------
sessionInfo()

