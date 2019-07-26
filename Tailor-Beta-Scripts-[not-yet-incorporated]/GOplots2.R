###############################################################################
# GOplots for Significantly differentially expressed features: method 1
############################################################################
# Libraries
###########################################################################
#source("http://bioconductor.org/biocLite.R")
#biocLite("topGO");biocLite("org.Hs.eg.db");biocLite("ReactomePA")

library(cummeRbund)
library(gage);library(pathview)
library(biomaRt);library(org.Hs.eg.db)
library(reactome.db);library(ReactomePA)
library(GO.db);library(GOstats);library(topGO)
library(KEGG.db);library(KEGGgraph);library(keggorthology)
library(clusterProfiler);library(DOSE)
library(limma)
library(GSEABase); library(STRINGdb)
library(igraph);source("http://igraph.sf.net")
## ----KEGG Download------------------------------------------------------------
hsa.kegg.code<-search_kegg_organism('hsa', by='kegg_code')
hsa.kegg.code
#go.abund
HSA.kegg <- search_kegg_organism('Homo sapiens', by='scientific_name')
dim(HSA.kegg)
head(HSA.kegg)
hsa.kegg<-download_KEGG(species="hsa", keggType = "KEGG", keyType = "kegg")
hsa.kegg$KEGGPATHID2EXTID[1:10,1]
hsa.kegg$KEGGPATHID2EXTID[1:10,2]
length(hsa.kegg)
head(hsa.kegg)

###############################################################################
# GOplots for Significantly differentially expressed features: method 1
############################################################################
## source('openGraphSaveGraph.r');
options(error=traceback) # causes a traceback to appear if there is an error, and the traceback has the line number
# prefixed by #

### original  R code from Drew'
## options(echo=TRUE)
## args = commandArgs(trailingOnly = T)
args = commandArgs(TRUE)

for (i in 1:length(args)) {
  eval(parse(text=args[[i]]))
}

###################################################
### code chunk number 1: init, 2: loadLib, 3: Read
###################################################
options(width=65)
library(cummeRbund)
#refgtf = Sys.getenv("CUMMERBUND_INPUT_GTF")
#genome_path = Sys.getenv("GENOME_PATH")
#gen_v = Sys.getenv("R_GENOME")
#cuff<-readCufflinks(gtfFile=refgtf,genome=genome_path,rebuild=T)
cuff<-readCufflinks()
cuff

###################################################
### for testing
###################################################
# default
#cuff<-readCufflinks(gtfFile="GRCh38_genes.gtf",genome="genome.fa",rebuild=F)
#cuff
over="LUTS"
under="CTRL"

###############################################################################
# Features, Counts, and all inclusive tables
# get gene symbols for CDS TSS Isoforms features
################################################################################

runInfo(cuff);
repdata<-replicates(cuff)
samps<-repdata$rep_name

genes_exp.diff<-diffData(isoforms(cuff))
head(genes_exp.diff)
g.rep.matrix<-repFpkmMatrix(isoforms(cuff))
head(g.rep.matrix)
gene.xloc.matrix<-featureNames(isoforms(cuff))
head(gene.xloc.matrix)
gene.list<-getGenes(cuff,geneId = gene.xloc.matrix)
gene.list
gene_annotation_data<-featureNames(gene.list@isoforms)
dim(gene_annotation_data)
head(gene_annotation_data)

reps<-colnames(g.rep.matrix)
reps

mySigGenes<-getSig(cuff,x=over,y=under,alpha=.05,level='genes')
head(mySigGenes)
length(mySigGenes)
sigGenes<-getGenes(cuff, mySigGenes)
length(sigGenes)

sig_genes_exp.diff<-diffData(sigGenes)

sig_genes_annot<-sigGenes@annotation
head(sig_genes_exp.diff)
head(sig_genes_annot)

sig_genes_exp.diff["gene_id"]<-sig_genes_annot["gene_short_name"]
head(sig_genes_exp.diff)
#genes_exp.diff<-as.data.frame(cbind(genes=sig_genes_exp.diff$gene_id,logFC=sig_genes_exp.diff$log2_fold_change,p_value=sig_genes_exp.diff$p_value,
#                                    q_value=sig_genes_exp.diff$q_value))
q_value=sig_genes_exp.diff$q_value
genes_exp.diff<-as.data.frame(q_value,row.names = sig_genes_exp.diff$gene_id)
genes_exp.diff<-as.data.frame(cbind(logFC=sig_genes_exp.diff$log2_fold_change,p_value=sig_genes_exp.diff$p_value,
                                    q_value=sig_genes_exp.diff$q_value),row.names=sig_genes_exp.diff$gene_id)
head(genes_exp.diff)
dim(genes_exp.diff)
rownames(genes_exp.diff)
#gene_exp.diff<-unique(rownames(genes_exp.diff))
#gene.exp.diff<-genes_exp.diff[gene_exp.diff,]

over.group<-grep(pattern = over, x = colnames(g.rep.matrix),ignore.case = T)
under.group<-grep(pattern = under, x = colnames(g.rep.matrix),ignore.case = T)

geneList<-mapIds(x = org.Hs.eg.db,
                 keys =  rownames(genes_exp.diff),
                 column = "ENTREZID",
                 keytype = "SYMBOL",
                 multiVals="first")

xx <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
topDiffGenes <- function(allScore) {
  return(allScore < 0.01)}

# Significantly Differentially Expressed by the Under group
# i.e Luts over Ctrl --> Ctrl = Under group, Under group highly
# expresses and No expression in the Over group
loOVERbpGOdata <- new("topGOdata",ontology = "BP",allGenes = genes_exp.diff, nodeSize = 5,
                      annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                      ID = "symbol")

loOVERbpGOdata <- new("topGOdata",ontology = "BP",allGenes = Qval.hi_exprUNDER, nodeSize = 5,
                      annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                      ID = "symbol")
HIunder.MFGOdata <- new("topGOdata",ontology = "MF",allGenes = Qval.hi_exprUNDER, nodeSize = 5,
                        annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                        ID = "symbol")
HIunder.CCGOdata <- new("topGOdata",ontology = "CC",allGenes = Qval.hi_exprUNDER, nodeSize = 5,
                        annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                        ID = "symbol")

LOunder.BPGOdata <- new("topGOdata",ontology = "BP",allGenes = Qval.hi_exprOVER, nodeSize = 5,
                        annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                        ID = "symbol")
hiOVERmfGOdata <- new("topGOdata",ontology = "MF",allGenes = Qval.hi_exprOVER, nodeSize = 5,
                      annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                      ID = "symbol")

hiOVERccGOdata <- new("topGOdata",ontology = "CC",allGenes = Qval.hi_exprOVER, nodeSize = 5,
                      annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                      ID = "symbol")

LOover.BPtKS <- runTest(loOVERbpGOdata, algorithm = "classic", statistic = "ks")
LOover.BPFisher <- runTest(loOVERbpGOdata, algorithm = "classic", statistic = "fisher")
LOover.BPtKS.elim <- runTest(loOVERbpGOdata, algorithm = "elim", statistic = "ks")

LOunder.BPtKS <- runTest(LOunder.BPGOdata, algorithm = "classic", statistic = "ks")
LOunder.BPFisher <- runTest(LOunder.BPGOdata, algorithm = "classic", statistic = "fisher")
LOunder.BPtKS.elim <- runTest(LOunder.BPGOdata, algorithm = "elim", statistic = "ks")

HIover.MFtKS <- runTest(hiOVERmfGOdata, algorithm = "classic", statistic = "ks")
HIover.MFFisher <- runTest(hiOVERmfGOdata, algorithm = "classic", statistic = "fisher")
HIoverMFtKS.elim <- runTest(hiOVERmfGOdata, algorithm = "elim", statistic = "ks")

HIunder.MFtKS <- runTest(HIunder.MFGOdata, algorithm = "classic", statistic = "ks")
HIunder.MFFisher <- runTest(HIunder.MFGOdata, algorithm = "classic", statistic = "fisher")
HIunder.MFtKS.elim <- runTest(HIunder.MFGOdata, algorithm = "elim", statistic = "ks")

HIover.CCtKS <- runTest(hiOVERccGOdata, algorithm = "classic", statistic = "ks")
HIover.CCFisher <- runTest(hiOVERccGOdata, algorithm = "classic", statistic = "fisher")
HIover.CCtKS.elim <- runTest(hiOVERccGOdata, algorithm = "elim", statistic = "ks")

HIunder.CCtKS <- runTest(HIunder.CCGOdata, algorithm = "classic", statistic = "ks")
HIunder.CCFisher <- runTest(HIunder.CCGOdata, algorithm = "classic", statistic = "fisher")
HIunder.CCtKS.elim <- runTest(HIunder.CCGOdata, algorithm = "elim", statistic = "ks")

pdf("GOplots_grch38.CTRL_highly_expressed.pdf")

showSigOfNodes(loOVERbpGOdata, score(LOover.BPFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(loOVERbpGOdata, score(LOover.BPtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(loOVERbpGOdata, score(LOover.BPtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(loOVERbpGOdata, LOover.BPFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(loOVERbpGOdata, LOover.BPtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(loOVERbpGOdata, LOover.BPtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(HIunder.MFGOdata, score(HIunder.MFFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(HIunder.MFGOdata, score(HIunder.MFtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(HIunder.MFGOdata, score(HIunder.MFtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(HIunder.MFGOdata, HIunder.MFFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(HIunder.MFGOdata, HIunder.MFtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(HIunder.MFGOdata, HIunder.MFtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(HIunder.CCGOdata, score(HIunder.CCFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(HIunder.CCGOdata, score(HIunder.CCtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(HIunder.CCGOdata, score(HIunder.CCtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(HIunder.CCGOdata, HIunder.CCFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(HIunder.CCGOdata, HIunder.CCtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(HIunder.CCGOdata, HIunder.CCtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

dev.off()
pdf("GOplots_grch38.LUTS_highly_expressed.pdf")

showSigOfNodes(LOunder.BPGOdata, score(LOunder.BPFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(LOunder.BPGOdata, score(LOunder.BPtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(LOunder.BPGOdata, score(LOunder.BPtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(LOunder.BPGOdata, LOunder.BPFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(LOunder.BPGOdata, LOunder.BPtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(LOunder.BPGOdata, LOunder.BPtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(hiOVERmfGOdata, score(HIover.MFFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(hiOVERmfGOdata, score(HIover.MFtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(hiOVERmfGOdata, score(HIoverMFtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(hiOVERmfGOdata, HIover.MFFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(hiOVERmfGOdata, HIoverMFtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(hiOVERmfGOdata, HIover.MFtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(hiOVERccGOdata, score(HIover.CCFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(hiOVERccGOdata, score(HIover.CCtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(hiOVERccGOdata, score(HIover.CCtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(hiOVERccGOdata, HIover.CCFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(hiOVERccGOdata, HIover.CCtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(hiOVERccGOdata, HIover.CCtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

dev.off()

## ----eval = FALSE--------------------------------------------------------
library(EnrichmentBrowser); library(gage)
hsa.pwys <- download.kegg.pathways("hsa")
hsa.kegg.gs <- get.kegg.genesets("hsa")
head(hsa.pwys);head(hsa.kegg.gs);data(gse16873)
hsa.kegg.sets<-kegg.gsets(species = "hsa",id.type = "kegg")
hsa.kegg.sets
## ----eval = FALSE--------------------------------------------------------
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
columns(org.Hs.eg.db)

#ENSEMBL ENTREZID GO
eg2np <- bitr_kegg(siggene_exp.diff$gene_id, fromType='ncbi-geneid', toType='kegg', organism='hsa')
#bitr_kegg("Z5100", fromType="kegg", toType='ncbi-proteinid', organism='ece')
#bitr_kegg("Z5100", fromType="kegg", toType='uniprot', organism='ece')

siggene.df <- bitr(siggene_exp.diff$gene_id, fromType = "SYMBOL",toType ="ENTREZID",OrgDb = org.Hs.eg.db)
dim(siggene.df)
sigENTREZgenes<-siggene.df$ENTREZID

go.len <- goana(siggene.df, geneid = sigENTREZgenes,FDR = 0.05, trend = FALSE,species = "Hs", plot=T)
topGO(go.len)

sig.kegga.gene <- kegga(siggene.df$ENTREZID, species.KEGG="hsa") # equivalent to previous
topKEGG(sig.kegga.gene)

plotGOgraph(go.len)

## ----fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8----
gseaplot(sig.kegga.gene, geneSetID = "hsa")


sigGOcc <- groupGO(gene=sigENTREZgenes, OrgDb=org.Hs.eg.db, ont="CC",
                   level    = 3,readable = TRUE)
sigGOmf <- groupGO(gene=sigENTREZgenes, OrgDb=org.Hs.eg.db, ont="MF",
                   level    = 3,readable = TRUE)
sigGObp <- groupGO(gene=sigENTREZgenes, OrgDb=org.Hs.eg.db, ont="BP",
                   level    = 3,readable = TRUE)

barplot(sigGOcc, drop=TRUE, showCategory=12,args.legend="Sig Diff Expr Genes Grouping by GO-CC")
barplot(sigGOmf, drop=TRUE, showCategory=12,args.legend="Sig Diff Expr Genes Grouping by GO-MF")
barplot(sigGObp, drop=TRUE, showCategory=12,args.legend="Sig Diff Expr Genes Grouping by GO-BP")

sigKEGGgenes <- enrichKEGG(gene = sigENTREZgenes,organism='hsa',pvalueCutoff = 0.05)
head(sigKEGGgenes)
barplot(sigKEGGgenes, drop=TRUE, showCategory=12)
dotplot(sigKEGGgenes)
cnetplot(sigKEGGgenes, categorySize="pvalue",foldChange) # ,wt1021.vs.wt1021B.kegg
enrichMap(sigKEGGgenes)

sigKEGG.df <- bitr(siggene_exp.diff$gene_id, fromType = "SYMBOL",toType ="PATH",OrgDb = org.Hs.eg.db)

kk2 <- gseKEGG(geneList     = sigENTREZgenes,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
## ------------------------------------------------------------------------
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)

## ------------------------------------------------------------------------
#kk3 <- gseKEGG(geneList     = gene,organism     = 'hsa',nPerm        = 1000,
#               minGSSize    = 120, pvalueCutoff = 0.05, verbose      = FALSE)

#kk3 <- gseKEGG(geneList = geneList, organism = 'hsa',nPerm = 1000,
#               minGSSize = 120,pvalueCutoff = 0.05,verbose = FALSE)
#head(kk2)

## ----eval = FALSE--------------------------------------------------------
mkk <- enrichMKEGG(gene = gene,
                   organism = 'hsa')
barplot(mkk, drop=TRUE, showCategory=12)

## ----eval=FALSE----------------------------------------------------------
mkk2 <- gseMKEGG(geneList = geneList,
                 species = 'hsa')

## ----eval=FALSE----------------------------------------------------------
## david <- enrichDAVID(gene = gene,
##                      idType = "ENTREZ_GENE_ID",
##                      listType = "Gene",
##                      annotation = "KEGG_PATHWAY",
##                      david.user = "clusterProfiler@hku.hk")

## ----fig.height=5, fig.width=9-------------------------------------------
barplot(ggo, drop=TRUE, showCategory=12)

## ----fig.height=5, fig.width=8-------------------------------------------
barplot(mkk, showCategory=8)
## ------------------------------------------------------------------------
dotplot(mkk)
## ----fig.cap="enrichment map of enrichment result", fig.align="center", fig.height=16, fig.width=16, eval=FALSE----
enrichMap(mkk)
## ## categorySize can be scaled by 'pvalue' or 'geneNum'
cnetplot(mkk, categorySize="pvalue", foldChange=geneList)

## ----fig.height=12, fig.width=8------------------------------------------
plotGOgraph(ego)

## ----fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8----
gseaplot(kk, geneSetID = "sme")


## ----eval=FALSE----------------------------------------------------------
library("pathview")
hsa04110 <- pathview(gene.data  = geneList,
                     species    = "sme",
                     limit      = list(gene=max(abs(geneList)), cpd=1))

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
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"
# [8] "EVIDENCE"     "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"        "IPI"          "MAP"
#[15] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"         "PROSITE"
#[22] "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIGENE"      "UNIPROT"
eg = bitr(glist, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg)
uniprot_ids <- bitr(glist, fromType="SYMBOL", toType=c("UNIPROT"), OrgDb="org.Hs.eg.db")
head(uniprot_ids)
refseq_ids <- bitr(glist, fromType="SYMBOL", toType=c("REFSEQ"), OrgDb="org.Hs.eg.db")
head(refseq_ids)
go_ids <- bitr(glist, fromType="SYMBOL", toType=c("UCSCKG"), OrgDb="org.Hs.eg.db")
head(go_ids)
go_ids <- bitr(glist, fromType="SYMBOL", toType=c("GOALL"), OrgDb="org.Hs.eg.db")
head(go_ids)
#eg2np <- bitr_kegg(glist, fromType='ncbi-geneid', toType='kegg', organism='hsa')
#bitr_kegg("Z5100", fromType="kegg", toType='ncbi-proteinid', organism='ece')
#bitr_kegg("Z5100", fromType="kegg", toType='uniprot', organism='ece')
library(DOSE)
gene.df <- bitr(AB.vs.wt1021.DE$GeneSymbol, fromType = "SYMBOL",toType = c("ENTREZID", "KEGG"),OrgDb=Org.Hs.egOMIM2EG@datacache)
gene.df <- bitr(glist, fromType = "SYMBOL",toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
str(gene.df)
entrezgenes<-gene.df[,"ENTREZID"]
ggo <- groupGO(gene=entrezgenes, OrgDb=org.Hs.eg.db, ont="CC",
               level    = 3,readable = TRUE)
head(ggo)
kk <- enrichKEGG(gene = entrezgenes,organism='hsa',pvalueCutoff = 0.05)
head(kk)

gene.df <- bitr(AB.vs.wt1021.DE$GeneSymbol, fromType = "SYMBOL",toType = c("ENTREZID", "KEGG"),OrgDb=Org.Hs.egOMIM2EG@datacache)

S.me1021 <- enrichKEGG(gene = geneList,organism='sme',pvalueCutoff = 0.05)
head(S.me1021)

#############################################################################


AB.vs.wt1021.npid <- bitr_kegg(AB.vs.wt1021.DE$GeneSymbol, fromType='kegg', toType='ncbi-proteinid', organism='sme',drop=T)
head(AB.vs.wt1021.npid)
dim(AB.vs.wt1021.npid)

AB.vs.wt1021.geneid <- bitr_kegg(AB.vs.wt1021.DE$GeneSymbol, fromType='kegg', toType='ncbi-geneid', organism='sme')
dim(AB.vs.wt1021.geneid)

## ------------------------------------------------------------------------
ego <- enrichGO(gene=entrezgenes,
                universe=names(geneList),
                OrgDb= org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

## ----eval=FALSE----------------------------------------------------------
## ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
##                 OrgDb         = org.Hs.eg.db,
## 		keytype       = 'ENSEMBL',
##                 ont           = "CC",
##                 pAdjustMethod = "BH",
##                 pvalueCutoff  = 0.01,
##                 qvalueCutoff  = 0.05)

## ----eval=FALSE----------------------------------------------------------
## ego2 <- setReadable(ego2, OrgDb = org.Hs.eg.db)

## ----eval=FALSE----------------------------------------------------------
## ego3 <- gseGO(geneList     = geneList,
##               OrgDb        = org.Hs.eg.db,
##               ont          = "CC",
##               nPerm        = 1000,
##               minGSSize    = 100,
##               maxGSSize    = 500,
##               pvalueCutoff = 0.05,
##               verbose      = FALSE)



## ------------------------------------------------------------------------
barcodeplot(AB.vs.wt1021B.DE[,8], index = AB.vs.wt1021B.DE[,7],index2 = AB.vs.wt1021B.DE[,8], col.bars = "dodgerblue",alpha=.01,
            labels = "LogFoldChange",xlab="FoldChange")

barcodeplot(A.vs.wt1021.DE[,8], index = A.vs.wt1021.DE[,7],index2 = A.vs.wt1021.DE[,8], col.bars = "dodgerblue",alpha=.01,
            labels = "LogFoldChange",xlab="FoldChange")

barcodeplot(wt1021B.vs.wt1021.DE[,8], index = wt1021B.vs.wt1021.DE[,7],index2 = wt1021B.vs.wt1021.DE[,8], col.bars = "dodgerblue",alpha=.01,
            labels = "LogFoldChange",xlab="FoldChange")

barcodeplot(A.vs.AB.DE[,8], index = A.vs.AB.DE[,7],index2 = A.vs.AB.DE[,8], col.bars = "dodgerblue",alpha=.01,
            labels = "LogFoldChange",xlab="FoldChange")

barcodeplot(AB.vs.wt1021.DE[,8], index = AB.vs.wt1021.DE[,7],index2 = AB.vs.wt1021.DE[,8], col.bars = "dodgerblue",alpha=.01,
            labels = "LogFoldChange",xlab="FoldChange")

## ------------------------------------------------------------------------
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)

## ----eval = FALSE--------------------------------------------------------
sme.genes <- enrichMKEGG(gene = geneList,
                         organism = 'sme')

## ----eval=FALSE----------------------------------------------------------
## mkk2 <- gseMKEGG(geneList = geneList,
##                  species = 'hsa')

## ----eval=FALSE----------------------------------------------------------
## david <- enrichDAVID(gene = gene,
##                      idType = "ENTREZ_GENE_ID",
##                      listType = "Gene",
##                      annotation = "KEGG_PATHWAY",
##                      david.user = "clusterProfiler@hku.hk")

## ------------------------------------------------------------------------
gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)

egmt <- enricher(gene, TERM2GENE=c5)
head(egmt)

egmt2 <- GSEA(geneList, TERM2GENE=c5, verbose=FALSE)
head(egmt2)

## ----fig.height=5, fig.width=9-------------------------------------------
barplot(mkk, drop=TRUE, showCategory=12)

## ----fig.height=5, fig.width=8-------------------------------------------
barplot(ego, showCategory=8)

## ------------------------------------------------------------------------
dotplot(mkk)

## ----
#fig.cap="enrichment map of enrichment result", fig.align="center", fig.height=16, fig.width=16, eval=FALSE----
# enrichMap(ego)

## ----fig.height=14, fig.width=14, eval=FALSE-----------------------------
## ## categorySize can be scaled by 'pvalue' or 'geneNum'
## cnetplot(ego, categorySize="pvalue", foldChange=geneList)

## ----fig.height=12, fig.width=8------------------------------------------
plotGOgraph(ego)

## ----fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8----
gseaplot(kk2, geneSetID = "hsa04145")

## ----eval=FALSE----------------------------------------------------------
## browseKEGG(kk, 'hsa04110')

## ----eval=FALSE----------------------------------------------------------
## library("pathview")
## hsa04110 <- pathview(gene.data  = geneList,
##                      pathway.id = "hsa04110",
##                      species    = "hsa",
##                      limit      = list(gene=max(abs(geneList)), cpd=1))

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

library(package = affyLib, character.only = TRUE)
## the distribution of the adjusted p-values
hist(geneList, 100)

## how many differentially expressed genes are:
sum(topDiffGenes(geneList))

## build the topGOdata class
GOdata <- new("topGOdata",ontology = "BP",
              allGenes = geneList,geneSel = topDiffGenes,
              annot = annFUN.db,affylib = affyLib)

## display the GOdata object
GOdata

##########################################################
## Examples on how to use the methods
##########################################################

## description of the experiment
description(GOdata)

## obtain the genes that will be used in the analysis
a <- genes(GOdata)
str(a)
numGenes(GOdata)

## obtain the score (p-value) of the genes
selGenes <- names(geneList)[sample(1:length(geneList), 10)]
gs <- geneScore(GOdata, whichGenes = selGenes)
print(gs)

## if we want an unnamed vector containing all the feasible genes
gs <- geneScore(GOdata, use.names = FALSE)
str(gs)

## the list of significant genes
sg <- sigGenes(GOdata)
str(sg)
numSigGenes(GOdata)

## to update the gene list
.geneList <- geneScore(GOdata, use.names = TRUE)
GOdata ## more available genes
GOdata <- updateGenes(GOdata, .geneList, topDiffGenes)
GOdata ## the available genes are now the feasible genes

## the available GO terms (all the nodes in the graph)
go <- usedGO(GOdata)
length(go)

## to list the genes annotated to a set of specified GO terms
sel.terms <- sample(go, 10)
ann.genes <- genesInTerm(GOdata, sel.terms)
str(ann.genes)

## the score for these genes
ann.score <- scoresInTerm(GOdata, sel.terms)
str(ann.score)

## to see the number of annotated genes
num.ann.genes <- countGenesInTerm(GOdata)
str(num.ann.genes)

## to summarise the statistics
termStat(GOdata, sel.terms)

#install te gprofiler package if it is not installed
install.packages("gProfileR")
sig_expr.genes<-sig_genes_exp.diff[which(sig_genes_exp.diff$q_value < 0.05 & sig_genes_exp.diff$log2_fold_change > 2),]
dim(sig_expr.genes)

#write.table(topgenes_qvalue005_mesenchymal,
#            "ImmunovsRest_immuno_significantgenes.txt",
#            col.names=FALSE, sep="\t", row.names=FALSE, quote=FALSE)
#

#Significantly up-regulated genes in immunoreactive samples have negative logFC and t-values.
sig_reactive.genes<-sig_genes_exp.diff[which(sig_genes_exp.diff$log2_fold_change < 0 & sig_genes_exp.diff$test_stat < 0),]
dim(sig_reactive.genes)

#write.table(topgenes_qvalue005_immunoreactive,
#            "ImmunovsRest_rest_significantgenes.txt",
#            col.names=FALSE, sep="\t", row.names=FALSE, quote=FALSE)

1.2.7 Create GSEA input list
20b. Create a rank file for GSEA. To run GSEA in pre-ranked mode, you need a two column RNK file with gene/protein/probe name (column 1) and the associated score (column 2). The first column should contain the same type of gene IDs used in the pathway gene-set (GMT) file. GSEA looks for enrichment in the top and bottom parts of the list, ranking the file using the t-statistic. The t-statistic indicates the strength of differential expression and is used in the p-value calculation. Other scores indicating the strength of differential expression may be used as well. GSEA ranks the most up-regulated genes at the top of the list and the most down-regulated at the bottom of the list. Genes at the top of the list are more highly expressed in class A compared to class B, while genes at the bottom of the list are higher in class B. In this workflow, a positive t-value means a higher expression of a gene in the Mesenchymal samples compared to the Immunoreactive samples (variable constrastnm). The following commands create a data frame with gene IDs and t-statistics, remove lines with missing gene IDs, and store the result as a RNK file. An additional step is usually required in analysis of Affymetrix microarray data as genes are represented with multiple probesets. The most significant probeset or average probeset score may be considered for every gene.

sig_reactive.genes<-sig_genes_exp.diff[which(sig_genes_exp.diff$log2_fold_change < 0 & sig_genes_exp.diff$test_stat < 0),]
dim(sig_reactive.genes)
ranks <- data.frame(geneID=sig_reactive.genes[,"gene_id"],t_stat=sig_reactive.genes[,"test_stat"],
                    stringsAsFactors=F)
ranks <- ranks[which(ranks[,"geneID"] != ""),]
head(ranks)

#1.2.8 Create expression file
#Create an expression file for the enrichment map and save files to the home folder of the analysis. The expression file contains the gene IDs as the first column gene description as the second column and the expression values for each sample as the additional columns. Gene IDs should correspond to the first column of the rank file. The text files will be saved on your computer in the directory specified at the beginning of the script using setwd(). The .rnk, .cls and .txt are all tab delimited files that can be viewed in spreadsheet or in a text editor.

library(biomaRt)

mart = useMart(biomart="ensembl") #, dataset="hsapiens_gene_ensembl")

genes = getBM(attributes = c( 'hgnc_symbol'), filters='hgnc_symbol',
              values=ranks$geneID, mart=mart);
genes$description = gsub("\\[Source.*", "", genes$description);

EM_expressionFile <- merge(genes,expressionMatrix,  all.y=TRUE,by.x=1, by.y=0)
colnames(EM_expressionFile)[1] <- "Name"
colnames(EM_expressionFile)[2] <- "Description"
write.table(EM_expressionFile, "TCGA_OV_expression.txt",
            col.name=TRUE, sep="\t", row.names=FALSE, quote=FALSE)

g.under.matrix<-g.rep.matrix[,under.group]
head(g.under.matrix)
dim(g.under.matrix)

g.over.matrix<-g.rep.matrix[,over.group]
head(g.over.matrix)
dim(g.over.matrix)

genes.reps.df<-cbind(gene_annotation_data, g.rep.matrix)
head(genes.reps.df)
genes.reps.df<-na.omit(genes.reps.df[unique(genes.reps.df$gene_short_name),])
head(genes.reps.df)

rownames(g.under.matrix)<-gene_annotation_data$tracking_id
head(g.under.matrix)
rownames(g.over.matrix)<-gene_annotation_data$tracking_id
head(g.over.matrix)


genes<-genes.reps.df$gene_short_name
rownames(genes_exp.diff)
foldchange<-genes_exp.diff[,"logFC"]
qval<-genes_exp.diff[,"q_value"]

under.inf = which(foldchange == "-Inf" & qval < 0.05)
over.inf = which(foldchange == "Inf" & qval < 0.05)
over.in= which(foldchange > 0 & qval < 0.05)
under.in = which(foldchange < 0 & qval < 0.05)

HIexp.inOVER<-as.data.frame(rbind(genes_exp.diff[over.inf,],genes_exp.diff[over.in,]))
head(HIexp.inOVER)
dim(HIexp.inOVER)
Qval.hi_exprOVER<-HIexp.inOVER[,"q_value"]
names(Qval.hi_exprOVER)<-rownames(HIexp.inOVER) #[,"genes"]

HIexp.inUNDER<-as.data.frame(rbind(genes_exp.diff[under.inf,],genes_exp.diff[under.in,]))
head(HIexp.inUNDER)
dim(HIexp.inUNDER)
Qval.hi_exprUNDER<-HIexp.inUNDER[,"q_value"]
names(Qval.hi_exprUNDER)<-rownames(HIexp.inUNDER) #[,"genes"]

#QvalnoOVER #QvalnoUNDER
ENTREZQval.hi_exprOVER<-mapIds(x = org.Hs.eg.db,
                               keys =  names(Qval.hi_exprOVER),
                               column = "ENTREZID",
                               keytype = "SYMBOL",
                               multiVals="first")

ENTREZQval.hi_exprUNDER<-mapIds(x = org.Hs.eg.db,
                                keys =  names(Qval.hi_exprUNDER),
                                column = "ENTREZID",
                                keytype = "SYMBOL",
                                multiVals="first")

ENTREZsiggenes<-mapIds(x = org.Hs.eg.db,
                       keys = sig_genes_exp.diff$gene_id,
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals="first")

NOexp.inUNDER<-as.data.frame(genes_exp.diff[over.inf,])
NOexp.inOVER<-as.data.frame(genes_exp.diff[under.inf,])
QvalnoOVER<-NOexp.inOVER[,"q_value"]
names(QvalnoOVER)<-rownames(NOexp.inOVER) #[,"gene_short_name"]
QvalnoUNDER<-NOexp.inUNDER[,"q_value"]
names(QvalnoUNDER)<-rownames(NOexp.inUNDER) #[,"gene_short_name"]

#LOexp.inUNDER<-as.data.frame(gene_exp.diff[over.in,])
#LOexp.inOVER<-as.data.frame(gene_exp.diff[under.in,])
OVERupexp<-as.data.frame(genes_exp.diff[over.in,])
QvalOVERhi<-OVERupexp[,"q_value"]
names(QvalOVERhi)<-rownames(OVERupexp) #[,"gene_short_name"]
UNDERupexp<-as.data.frame(genes_exp.diff[under.in,])
QvalUNDERhi<-UNDERupexp[,"q_value"]
names(QvalUNDERhi)<-rownames(UNDERupexp) #[,"gene_short_name"]

#QvalOVERhi #QvalUNDERhi
ENTREZQvalUNDERhi<-mapIds(x = org.Hs.eg.db,
                          keys = names(QvalUNDERhi),
                          column = "ENTREZID",
                          keytype = "SYMBOL",
                          multiVals="first")

ENTREZQvalOVERhi<-mapIds(x = org.Hs.eg.db,
                         keys = names(QvalOVERhi),
                         column = "ENTREZID",
                         keytype = "SYMBOL",
                         multiVals="first")

under.hilog<-sig_genes_exp.diff[,"log2_fold_change"]
names(under.hilog)<-ENTREZsiggenes
under.hiqval<-sig_genes_exp.diff[,"q_value"]
names(under.hiqval)<-ENTREZQvalUNDERhi
under.hipval<-(sig_genes_exp.diff[,"p_value"])
names(under.hipval) <- ENTREZQvalUNDERhi

over.hilog<-c(sig_genes_exp.diff[,"log2_fold_change"])
names(over.hilog)<-(ENTREZQvalOVERhi)
over.hiqval<-sig_genes_exp.diff[,"q_value"]
names(over.hiqval)<-ENTREZQvalOVERhi
over.hipval<-sig_genes_exp.diff[,"p_value"]
names(over.hipval)<-ENTREZQvalOVERhi

