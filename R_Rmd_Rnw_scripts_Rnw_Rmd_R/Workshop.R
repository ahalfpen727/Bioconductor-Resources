## ----setup, echo=FALSE---------------------------------------------------
knitr::opts_chunk$set(cache=TRUE, autodep=TRUE)
suppressPackageStartupMessages({
    library(ggplot2)
    library(Biostrings)
    library(GenomicRanges)
    library(SummarizedExperiment)

    library(DESeq2)
    library(VariantFiltering)

    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(TxDb.Hsapiens.NCBI.GRCh38)

    library(BSgenome.Hsapiens.UCSC.hg38)
    library(BSgenome.Hsapiens.UCSC.hg19)
   library(BSgenome.Hsapiens.NCBI.GRCh38)
#   library(BSgenome.Hsapiens.1000genomes.hs37d5)
   library(org.Hs.eg.db)
    
    library(airway)
    library(RNAseqData.HNRNPC.bam.chr14)
})

