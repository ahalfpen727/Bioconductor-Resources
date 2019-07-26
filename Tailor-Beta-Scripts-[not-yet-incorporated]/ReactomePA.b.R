## ----style, echo=FALSE, results="asis", message=FALSE----------------------
BiocStyle::markdown()
knitr::opts_chunk$set(tidy = FALSE,
                      warning = FALSE,
                      message = FALSE)

## ----echo=FALSE, results='hide', message=FALSE-----------------------------
library(DOSE)
library(org.Hs.eg.db)

## --------------------------------------------------------------------------
library(DOSE)
data(geneList)
y <- gseDO(geneList,
           nPerm         = 100,
           minGSSize     = 120,
           pvalueCutoff  = 0.2,
           pAdjustMethod = "BH",
           verbose       = FALSE)
head(y, 3)

## --------------------------------------------------------------------------
ncg <- gseNCG(geneList,
              nPerm         = 100,
              minGSSize     = 120,
              pvalueCutoff  = 0.2,
              pAdjustMethod = "BH",
              verbose       = FALSE)
ncg <- setReadable(ncg, 'org.Hs.eg.db')
head(ncg, 3)

## --------------------------------------------------------------------------
dgn <- gseDGN(geneList,
              nPerm         = 100,
              minGSSize     = 120,
              pvalueCutoff  = 0.2,
              pAdjustMethod = "BH",
              verbose       = FALSE)
dgn <- setReadable(dgn, 'org.Hs.eg.db')
head(dgn, 3)

## ----style, echo=FALSE, results="asis", message=FALSE--------------------
knitr::opts_chunk$set(tidy = FALSE,
                      warning = FALSE,
                      message = FALSE)

library(org.Hs.eg.db)
library(DOSE)
library(ReactomePA)
library(affy)
library(estrogen)
library(vsn)
library(voom)
library(genefilter)
source("http://bioconductor.org/biocLite.R")
biocLite(c("vsn","voom")

         datadir <- system.file("extdata", package="estrogen")
         dir(datadir)
         setwd(datadir)

         #Read in the phenotype data and the raw expression 'CEL' files
         pd <- read.AnnotatedDataFrame("estrogen.txt", header=TRUE, sep="", row.names=1)
         a <- ReadAffy(filenames=rownames(pData(pd)), phenoData=pd, verbose=TRUE)

         #Normalise the data
         x <- expresso(
           a,
           bg.correct=FALSE,
           normalize.method="vsn",
           normalize.param=list(subsample=1000),
           pmcorrect.method="pmonly",
           summary.method="medianpolish"
         )

         #NB - if the vsn normalisation does not function, use:
         #x <- expresso(a, bgcorrect.method="rma", normalize.method="constant", pmcorrect.method="pmonly", summary.method="avgdiff")

         #Remove control probes
         controlProbeIdx <- grep("^AFFX", rownames(x))
         x <- x[-controlProbeIdx,]

         #Identify genes of significant effect
         lm.coef <- function(y) lm(y ~ estrogen * time.h)$coefficients
         eff <- esApply(x, 1, lm.coef)
         effectUp <- names(sort(eff[2,], decreasing=TRUE)[1:25])
         effectDown <- names(sort(eff[2,], decreasing=FALSE)[1:25])
         main.effects <- c(effectUp, effectDown)

         #Filter the expression set object to include only genes of significant effect
         estrogenMainEffects <- exprs(x)[main.effects,]
         ## ------------------------------------------------------------------------
         library(ReactomePA)
         data(geneList)
         de <- names(geneList)[abs(geneList) > 1.5]
         head(de)
         x <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
         head(as.data.frame(x))

         ## ----fig.height=6, fig.width=12------------------------------------------
         barplot(x, showCategory=8)

         ## ----fig.height=6, fig.width=12------------------------------------------
         dotplot(x, showCategory=15)

         ## ----fig.height=10, fig.width=10-----------------------------------------
         emapplot(x)

         ## ----fig.height=8, fig.width=8-------------------------------------------
         cnetplot(x, categorySize="pvalue", foldChange=geneList)

         ## ----fig.height=8, fig.width=13, eval=FALSE------------------------------
         require(clusterProfiler)
         #  data(gcSample)
         #  res <- compareCluster(gcSample, fun="enrichPathway")
         #  dotplot(res)

         ## ------------------------------------------------------------------------
         y <- gsePathway(geneList, nPerm=10000,
                         pvalueCutoff=0.2,
                         pAdjustMethod="BH", verbose=FALSE)
         res <- as.data.frame(y)
         head(res)

         ## ----fig.height=8, fig.width=8-------------------------------------------
         emapplot(y, color="pvalue")

         ## ----fig.height=7, fig.width=10------------------------------------------
         gseaplot(y, geneSetID = "R-HSA-69242")

         ## ----fig.height=8, fig.width=8-------------------------------------------
         viewPathway("E2F mediated regulation of DNA replication", readable=TRUE, foldChange=geneList)

