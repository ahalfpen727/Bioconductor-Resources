## ----opts_chunk, echo=FALSE, results="hide"------------------------------
knitr::opts_chunk$set(cache=TRUE, cache.path="../../tmp/",
                      fig.width=10, fig.height=10, fig.path="../../tmp/")

## ----import_data, message=FALSE------------------------------------------
library(cogena)
data(Psoriasis)
# objects in the Psoriasis dataset.
# Note: label of interest should follow the control label as this
# will affect the direction of gene regulation.
# For instance use factor (c("Normal", "Cancer", "Normal"),
# levels=c("Normal", "Cancer")), instead of factor(c("Normal",
# "Cancer","Normal")) since "Cancer" is the label of interest.



## ----PD_data, echo=TRUE--------------------------------------------------
ls()

## ----gene_list_annotation------------------------------------------------

# KEGG Pathway gene set
annoGMT <- "c2.cp.kegg.v5.0.symbols.gmt.xz"
# GO biological process gene set
annoGObpGMT <- "c5.bp.v5.0.symbols.gmt.xz"
annoGOccGMT <- "c5.cc.v5.0.symbols.gmt.xz"
annoGOmfGMT <- "c5.mf.v5.0.symbols.gmt.xz"

annofile <- system.file("extdata", annoGMT, package="cogena")

# the number of clusters. It can be a vector.
# nClust <- 2:20
nClust <- 10
# Making factor "Psoriasis" behind factor "ct" means Psoriasis Vs Control
# is up-regualted
sampleLabel <- factor(sampleLabel, levels=c("ct", "Psoriasis"))

# the number of cores.
# ncore <- 8
ncore <- 2

# the clustering methods
# clMethods <- c("hierarchical","kmeans","diana","fanny","som","model",
# "sota","pam","clara","agnes") # All the methods can be used together.
clMethods <- c("hierarchical","pam")


# the distance metric
metric <- "correlation"

# the agglomeration method used for hierarchical clustering
# (hierarchical and agnes)
method <- "complete"


## ----cogena, results="hide"----------------------------------------------
# Co-expression Analysis
genecl_result <- coExp(DEexprs, nClust=nClust, clMethods=clMethods,
                       metric=metric, method=method, ncore=ncore)

# Enrichment (Pathway) analysis for the co-expressed genes
clen_res <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel)

## ----cogena_result_analysis----------------------------------------------
summary(clen_res)

## ----enrichment----------------------------------------------------------
# Here we consider the "pam" method and 10 clusters.
# Always make the number as character, please!
enrichment.table <- enrichment(clen_res, "pam", "10")

## ----heatmapCluster, fig.width=10, fig.height=7, fig.cap="Heatmap of expression profiling with clusters"----
# Always make the number as character, please!
heatmapCluster(clen_res, "pam", "10", maintitle="Psoriasis")

## ----heatmapPEI, fig.width=11, fig.height=7, fig.cap="KEGG pathway enrichment"----
# The enrichment score for 10 clusters, together with Down-regulated,
# Up-regulated and All DE genes. The values shown in Figure 2 is the -log2(FDR).
#
# Always make the number as character, please!

heatmapPEI(clen_res, "pam", "10", printGS=FALSE, maintitle="Pathway analysis for Psoriasis")


## ----drp-----------------------------------------------------------------
# A comprehensive way
# cmapDn100_cogena_result <- clEnrich(genecl_result,
# annofile=system.file("extdata", "CmapDn100.gmt.xz", package="cogena"),
# sampleLabel=sampleLabel)

# A quick way
# Based on the pathway analysis results
cmapDn100_cogena_result <- clEnrich_one(genecl_result, "pam", "10",
    annofile=system.file("extdata", "CmapDn100.gmt.xz", package="cogena"),
    sampleLabel=sampleLabel)

## ----drp_figure, fig.width=10, fig.height=7, fig.cap="Drug repositioning"----
heatmapPEI(cmapDn100_cogena_result, "pam", "10", printGS=FALSE,
           orderMethod = "7", maintitle="Drug repositioning for Psoriasis")

# Results based on cluster 5.
# heatmapPEI(cmapDn100_cogena_result, "pam", "10", printGS=FALSE,
#           orderMethod = "5", maintitle="Drug repositioning for Psoriasis")

# Results based on cluster 9, containing down-regulated genes.
# heatmapPEI(cmapUp100_cogena_result, "pam", "10", printGS=FALSE,
#           orderMethod = "9", maintitle="Drug repositioning for Psoriasis")


## ----drp_figure2, fig.width=10, fig.height=7, fig.cap="Drug repositioning (multi-instance merged)"----
heatmapCmap(cmapDn100_cogena_result, "pam", "10", printGS=FALSE,
            orderMethod = "7", maintitle="Drug repositioning for Psoriasis")

## ----geneInCluster-------------------------------------------------------
# Always make the number as character, please!
geneC <- geneInCluster(clen_res, "pam", "10", "4")
head(geneC)


## ----geneExpInCluster----------------------------------------------------
# Always make the number as character, please!
gec <- geneExpInCluster(clen_res, "pam", "10")
gec$clusterGeneExp[1:3, 1:4]
gec$label[1:4]

## ----corInCluster, fig.width=6, fig.height=6, fig.cap="Correlation of genes in a cluster"----
# Always make the number as character, please!
corInCluster(clen_res, "pam", "10", "10")

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

