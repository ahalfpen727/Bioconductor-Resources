## ----style, echo=FALSE, results="asis", message=FALSE--------------------
knitr::opts_chunk$set(tidy = FALSE,
                      warning = FALSE,
                      message = FALSE)
data(geneList)
de <- names(geneList)[abs(geneList) > 1.5]
head(de)
myGenes <- getGenes(cuff, sig_geneids)
myGenes

plot(expressionBarplot(myGenes,logMode=T,pseudocount=1))

## [1] "4312"  "8318"  "10874" "55143" "55388" "991"

x <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
head(as.data.frame(x))

## ----echo=FALSE, results='hide', message=FALSE---------------------------
library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(GSEABase)
library(clusterProfiler)

## ------------------------------------------------------------------------
x <- c("GPX3",  "GLRX",   "LBP",   "CRYAB", "DEFB1", "HCLS1",   "SOD2",   "HSPA2",
       "ORM1",  "IGFBP1", "PTHLH", "GPC3",  "IGFBP3","TOB1",    "MITF",   "NDRG1",
       "NR1H4", "FGFR3",  "PVR",   "IL6",   "PTPRM", "ERBB2",   "NID2",   "LAMB1",
       "COMP",  "PLS3",   "MCAM",  "SPP1",  "LAMC1", "COL4A2",  "COL4A1", "MYOC",
       "ANXA4", "TFPI2",  "CST6",  "SLPI",  "TIMP2", "CPM",     "GGT1",   "NNMT",
       "MAL",   "EEF1A2", "HGD",   "TCN2",  "CDA",   "PCCA",    "CRYM",   "PDXK",
       "STC1",  "WARS",  "HMOX1", "FXYD2", "RBP4",   "SLC6A12", "KDELR3", "ITM2B")
eg = bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg)

## ------------------------------------------------------------------------
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

## ------------------------------------------------------------------------
ids <- bitr(x, fromType="SYMBOL", toType=c("UNIPROT", "ENSEMBL"), OrgDb="org.Hs.eg.db")
head(ids)

## ------------------------------------------------------------------------
data(gcSample)
hg <- gcSample[[1]]
head(hg)

eg2np <- bitr_kegg(hg, fromType='kegg', toType='ncbi-proteinid', organism='hsa')
head(eg2np)

## ----eval=FALSE----------------------------------------------------------
#  bitr_kegg("Z5100", fromType="kegg", toType='ncbi-geneid', organism='ece')

## ------------------------------------------------------------------------
bitr_kegg("Z5100", fromType="kegg", toType='ncbi-proteinid', organism='ece')
bitr_kegg("Z5100", fromType="kegg", toType='uniprot', organism='ece')

## ----warning=FALSE-------------------------------------------------------
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
	   	toType = c("ENSEMBL", "SYMBOL"),
		OrgDb = org.Hs.eg.db)
head(gene.df)
ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

head(ggo)

## ------------------------------------------------------------------------
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
		readable      = TRUE)
head(ego)

## ----eval=FALSE----------------------------------------------------------
#  ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
#                  OrgDb         = org.Hs.eg.db,
#                  keyType       = 'ENSEMBL',
#                  ont           = "CC",
#                  pAdjustMethod = "BH",
#                  pvalueCutoff  = 0.01,
#                  qvalueCutoff  = 0.05)

## ----eval=FALSE----------------------------------------------------------
#  ego2 <- setReadable(ego2, OrgDb = org.Hs.eg.db)

## ----eval=FALSE----------------------------------------------------------
#  ego3 <- gseGO(geneList     = geneList,
#                OrgDb        = org.Hs.eg.db,
#                ont          = "CC",
#                nPerm        = 1000,
#                minGSSize    = 100,
#                maxGSSize    = 500,
#                pvalueCutoff = 0.05,
#                verbose      = FALSE)

## ------------------------------------------------------------------------
search_kegg_organism('ece', by='kegg_code')
ecoli <- search_kegg_organism('Escherichia coli', by='scientific_name')
dim(ecoli)
head(ecoli)

## ------------------------------------------------------------------------
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)

## ------------------------------------------------------------------------
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)

## ----eval = FALSE--------------------------------------------------------
#  mkk <- enrichMKEGG(gene = gene,
#                     organism = 'hsa')

## ----eval=FALSE----------------------------------------------------------
#  mkk2 <- gseMKEGG(geneList = geneList,
#                   species = 'hsa')

## ----eval=FALSE----------------------------------------------------------
#  david <- enrichDAVID(gene = gene,
#                       idType = "ENTREZ_GENE_ID",
#                       listType = "Gene",
#                       annotation = "KEGG_PATHWAY",
#                       david.user = "clusterProfiler@hku.hk")

## ------------------------------------------------------------------------
gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)

egmt <- enricher(gene, TERM2GENE=c5)
head(egmt)

egmt2 <- GSEA(geneList, TERM2GENE=c5, verbose=FALSE)
head(egmt2)

## ----fig.height=5, fig.width=9-------------------------------------------
barplot(ggo, drop=TRUE, showCategory=12)

## ----fig.height=5, fig.width=8-------------------------------------------
barplot(ego, showCategory=8)

## ------------------------------------------------------------------------
dotplot(ego)

## ----fig.cap="enrichment map of enrichment result", fig.align="center", fig.height=8, fig.width=8----
emapplot(ego)

## ----fig.height=8, fig.width=8-------------------------------------------
## categorySize can be scaled by 'pvalue' or 'geneNum'
cnetplot(ego, categorySize="pvalue", foldChange=geneList)

## ----fig.height=12, fig.width=8, eval=FALSE------------------------------
#  goplot(ego)

## ----fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8----
gseaplot(kk2, geneSetID = "hsa04145")

## ----eval=FALSE----------------------------------------------------------
#  browseKEGG(kk, 'hsa04110')

## ----eval=FALSE----------------------------------------------------------
#  library("pathview")
#  hsa04110 <- pathview(gene.data  = geneList,
#                       pathway.id = "hsa04110",
#                       species    = "hsa",
#                       limit      = list(gene=max(abs(geneList)), cpd=1))

## ------------------------------------------------------------------------
data(gcSample)
lapply(gcSample, head)

## ------------------------------------------------------------------------
ck <- compareCluster(geneCluster = gcSample, fun = "enrichKEGG")
head(as.data.frame(ck))

## ------------------------------------------------------------------------
mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"

formula_res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichKEGG")

head(as.data.frame(formula_res))

## ----fig.height=7, fig.width=9-------------------------------------------
dotplot(ck)

## ----fig.height=6, fig.width=10------------------------------------------
dotplot(formula_res)
dotplot(formula_res, x=~group) + ggplot2::facet_grid(~othergroup)

## ----echo=FALSE----------------------------------------------------------
sessionInfo()

