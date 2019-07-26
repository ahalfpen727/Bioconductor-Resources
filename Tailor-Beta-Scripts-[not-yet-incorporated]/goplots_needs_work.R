de <- row.names(sig.gene.fold.list)[abs(sig.gene.fold.list) > 1.5]
head(de)
length(de)

eg = bitr(de, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
sym.eg<-as.factor(eg$SYMBOL)
dim(eg)

symbols<-row.names(sig.gene.fold.list)
head(symbols)
length(symbols)

entrez.syms<-which(row.names(sig.gene.fold.list) %in% sym.eg)
length(entrez.syms)
head(sig.gene.fold.list)
foldchange<-sig.gene.fold.list[entrez.syms,]
SYMBOL<-symbols[entrez.syms]
siggene.fc.eg<-cbind(SYMBOL, foldchange)
head(siggene.fc.eg)
siggene.eg.fc<-merge(eg, siggene.fc.eg, by.x="SYMBOL", by.y="SYMBOL")
SigGeneList<-as.data.frame(siggene.eg.fc[,"foldchange"], row.names=siggene.eg.fc$ENTREZID)
colnames(SigGeneList)<-"log2.fold.change"
head(SigGeneList)
de<-row.names(SigGeneList)
head(de)
x <- enrichPathway(gene=SigGeneList,pvalueCutoff=0.05, readable=T)
head(as.data.frame(x))

## ----fig.height=6, fig.width=12------------------------------------------
barplot(x, showCategory=8)
dotplot(x, showCategory=15)
cnetplot(x, categorySize="pvalue", foldChange=fold.change)
SigGeneList<-SigGeneList[order(SigGeneList$log2.fold.change, decreasing = T),]
## ------------------------------------------------------------------------
y <- gsePathway(SigGeneList, nPerm=10000,
                pvalueCutoff=0.2,
                pAdjustMethod="BH", verbose=FALSE)
res <- as.data.frame(y); head(res)

## ----fig.height=8, fig.width=8-------------------------------------------
#emapplot(y, color="pvalue")
#gseaplot(y, geneSetID = "R-HSA-69242")
#viewPathway("E2F mediated regulation of DNA replication", readable=TRUE, foldChange=geneList)
gene.df <- bitr(names(geneList), fromType = "SYMBOL",toType = c("ENTREZID", "SYMBOL"),OrgDb = org.Hs.eg.db)

gene.df <- bitr(gene_diff_data$gene_id, fromType = "SYMBOL",toType = c("ENTREZID", "SYMBOL"),OrgDb = org.Hs.eg.db)
head(gene.df)
ggoCC <- groupGO(gene     = de,
                 OrgDb    = org.Hs.eg.db,
                 ont      = "CC",
                 level    = 3,
                 readable = TRUE)

ggoBP <- groupGO(gene     = de,
                 OrgDb    = org.Hs.eg.db,
                 ont      = "BP",
                 level    = 3,
                 readable = TRUE)

ggoMF <- groupGO(gene     = de,
                 OrgDb    = org.Hs.eg.db,
                 ont      = "MF",
                 level    = 3,
                 readable = TRUE)
head(ggo)
## ----fig.height=5, fig.width=9-------------------------------------------
barplot(ggoCC, drop=TRUE, showCategory=12)
barplot(ggoMF, drop=TRUE, showCategory=12)
barplot(ggoBP, drop=TRUE, showCategory=12)

barplot(ego, showCategory=8)
dotplot(ego)

## ----fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8----
gseaplot(kk2, geneSetID = "hsa04145")

## ------------------------------------------------------------------------
ego <- enrichGO(gene          = de,
                universe      = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
cnetplot(ego, categorySize="pvalue")

## categorySize can be scaled by 'pvalue' or 'geneNum'
cnetplot(ego, categorySize="pvalue", foldChange=geneList)

## ----fig.height=12, fig.width=8, eval=FALSE------------------------------
goplot(ego)

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
egoCC <- gseGO(geneList     = de,
               OrgDb        = org.Hs.eg.db,
               ont          = "CC",
               nPerm        = 1000,
               minGSSize    = 100,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
egoBP <- gseGO(geneList     = de,
               OrgDb        = org.Hs.eg.db,
               ont          = "BP",
               nPerm        = 1000,
               minGSSize    = 100,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
egoMF <- gseGO(geneList     = de,
               OrgDb        = org.Hs.eg.db,
               ont          = "MF",
               nPerm        = 1000,
               minGSSize    = 100,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
## ------------------------------------------------------------------------
kk <- enrichKEGG(gene         = de,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)

## ------------------------------------------------------------------------
kk2 <- gseKEGG(geneList     = de,
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
search_kegg_organism('ece', by='kegg_code')
ecoli <- search_kegg_organism('Escherichia coli', by='scientific_name')
dim(ecoli)
head(ecoli)

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
gseaplot(kk, geneSetID = "hsa04145")

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
dotplot(formula_res)
dotplot(formula_res, x=~group) + ggplot2::facet_grid(~othergroup)

