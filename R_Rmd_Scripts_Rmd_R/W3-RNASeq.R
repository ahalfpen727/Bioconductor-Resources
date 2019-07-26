## ----style, echo = FALSE, results = 'asis'--------------------------------------------------------
BiocStyle::markdown()
options(width=100, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")))

## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE--------------------------------------------
suppressPackageStartupMessages({
    library(airway)
    library(DESeq2)
})

## ----airway---------------------------------------------------------------------------------------
library("airway")
data(airway)
airway

## main components of SummarizedExperiment
head(assay(airway))
colData(airway)
rowRanges(airway)

## e.g., coordinated subset to include dex 'trt'  samples
airway[, airway$dex == "trt"]

## e.g., keep only rows with non-zero counts
airway <- airway[rowSums(assay(airway)) != 0, ]

## ----DESeqDataSet---------------------------------------------------------------------------------
library(DESeq2)
dds <- DESeqDataSet(airway, design = ~ cell + dex)

## ----DESeq-workflow-------------------------------------------------------------------------------
dds <- DESeq(dds)
dds

## ----DESeq-result---------------------------------------------------------------------------------
res <- results(dds)
res

