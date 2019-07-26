## ----style, echo = FALSE, results = 'asis'--------------------------------------------------------
options(width=100)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")))

## ----setup, echo=FALSE----------------------------------------------------------------------------
suppressPackageStartupMessages({
    library(airway)
    loadNamespace("matrixStats")
    library(DESeq2)
    library(gplots)
    library(ggplot2)
    library(RColorBrewer)
    library(org.Hs.eg.db)
})

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

## ----library-size---------------------------------------------------------------------------------
colSums(assay(airway))

## ----gene-expression------------------------------------------------------------------------------
means <- rowMeans(assay(airway))
xlim <- range(log(1 + means))
plot(density(log(1 + means)), xlim=xlim)
plot(density(log(1 + means[means > 1])), xlim=xlim)

## ----airway-euclidean-distance--------------------------------------------------------------------
d <- dist(t(log(1 + assay(airway))))
mds <- cmdscale(d)
plot(mds, pch=20, asp=1, cex=2)
plot(mds, pch=20, asp=1, cex=2, col=airway$cell)
plot(mds, pch=20, asp=1, cex=2, col=airway$dex)

## ----mean-var-------------------------------------------------------------------------------------
rowVars <- matrixStats::rowVars
plot(rowVars(1 + assay(airway)) ~ rowMeans(1 + assay(airway)), log="xy")

## -------------------------------------------------------------------------------------------------
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ cell + dex)

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
library(ggplot2)
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

