library(cummeRbund);library(ReactomePA); library(pathview)
#library(org.Mm.eg.db)
library(limma);library(DOSE);library(GO.db);library(org.Hs.eg.db)
library(topGO);library(GSEABase);library(clusterProfiler)
library(biomaRt);library(KEGG.db);library(reactome.db)
###############################################################################
# Features, Counts, and all inclusive tables 
# (CDS-p_ids, TSS-TSS ids, isoforms- XLOC and TCONS ids, and Genes XLOC ids)
# Significantly differentially expressed features: method 1
############################################################################
## ---CummeRbund and GSEA ------------------------------------------------------------------------------------------------
#source('~/R_Py_scripts/dcEnrichmentGO.R', echo=TRUE)
cuff<-readCufflinks(gtfFile="cuffcmp.combined.gtf",genome="genome.fa",rebuild=F)
cuff
refgtf<-read.table("cuffcmp.combined.gtf", sep=c("\t"), header=F)
head(refgtf)
genome="genome.fa"
over="LUTS"
under="CTRL"
###############################################################################
# Features, Counts, and all inclusive tables 
# get gene symbols for CDS TSS Isoforms features
################################################################################
runInfo(cuff);
repdata<-replicates(cuff)
reps<-repdata$rep_name
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

genes_exp.diff<-diffData(isoforms(cuff))
dim(genes_exp.diff)
g.rep.matrix<-repFpkmMatrix(isoforms(cuff))
head(g.rep.matrix)
gene.xloc.matrix<-featureNames(isoforms(cuff))
gene.list<-getGenes(cuff,geneId = gene.xloc.matrix)
gene.list
gene_annotation_data<-featureNames(gene.list@isoforms)
dim(gene_annotation_data)
gene.rep.matrix<-cbind2(gene_annotation_data, g.rep.matrix)
gene_exp.diff<-cbind2(gene_annotation_data, genes_exp.diff)
genes_exp.diff<-as.data.frame(cbind(gene_annotation_data, 
								   logFC=gene_exp.diff$log2_fold_change,
								   pvalue=gene_exp.diff$p_value,
								   qvalue=gene_exp.diff$q_value))

gene_exp.diff<-genes_exp.diff[unique(genes_exp.diff$gene_short_name),]
repdata<-replicates(cuff)
reps<-repdata$rep_name
samples<-repdata$sample_name
colnames(g.rep.matrix)<-samples
groups<-factor(g.rep.matrix,levels=0:1,labels=c(over, under))
g.over.matrix<-g.rep.matrix[samples == over]
g.under.matrix<-g.rep.matrix[samples == under]

colnames(g.rep.matrix)<-reps
rownames(g.rep.matrix)<-gene_annotation_data$tracking_id
genes.reps.df<-cbind(gene_annotation_data$gene_short_name, g.rep.matrix)
genes.reps.df<-gene.rep.matrix[unique(gene.rep.matrix$gene_short_name),]
genes<-genes.reps.df$gene_short_name


foldchange<-sig_genes_exp.diff[,"log2_fold_change"]
under.hi<-which(foldchange == "-Inf")
over.hi = which(foldchange == "Inf")
over.up = which(foldchange > 0, foldchange == "Inf", foldchange != "-Inf")
under.up = which(foldchange < 0,  foldchange == "-Inf")

NOexp.inUNDER<-as.data.frame(sig_genes_exp.diff[over.hi,])
NOexp.inOVER<-as.data.frame(sig_genes_exp.diff[under.hi,])

LOexp.inUNDER<-as.data.frame(sig_genes_exp.diff[over.up,])
LOexp.inOVER<-as.data.frame(sig_genes_exp.diff[under.up,])

QvalnoOVER<-NOexp.inOVER[,"q_value"]
names(QvalnoOVER)<-NOexp.inOVER[,"gene_id"]

QvalnoUNDER<-NOexp.inUNDER[,"q_value"]
names(QvalnoUNDER)<-NOexp.inUNDER[,"gene_id"]

OVERupexp<-as.data.frame(sig_genes_exp.diff[over.up,])
QvalOVERhi<-OVERupexp[,"q_value"]
names(QvalOVERhi)<-OVERupexp[,"gene_id"]

UNDERupexp<-as.data.frame(sig_genes_exp.diff[under.up,])
QvalUNDERhi<-UNDERupexp[,"q_value"]
names(QvalUNDERhi)<-UNDERupexp[,"gene_id"]

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

ENTREZsiggenes<-mapIds(x = org.Hs.eg.db,
                          keys = sig_genes_exp.diff$gene_id,
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

ENTREZQvalUNDER<-mapIds(x = org.Hs.eg.db,
                          keys =NOexp.inUNDER[,"gene_id"],
                          column = "ENTREZID",
                          keytype = "SYMBOL",
						multiVals="first")

ENTREZQvalOVER<-mapIds(x = org.Hs.eg.db,
                          keys =NOexp.inOVER[,"gene_id"],
                          column = "ENTREZID",
                          keytype = "SYMBOL",
					       multiVals="first" )

ENTREZids<-mapIds(x = org.Hs.eg.db,
                          keys = gene.list@isoforms@annotation$gene_short_name,
                          column = "ENTREZID",
                          keytype = "SYMBOL",
				          multiVals="first") 

sig_genes_exp<-cbind(ENTREZsiggenes,sig_genes_exp.diff[,"log2_fold_change"])
head(sig_genes_exp)

siggenelog<-sig_genes_exp.diff[,"log2_fold_change"]
names(siggenelog)<-ENTREZsiggenes
siggeneqval<-sig_genes_exp.diff[,"q_value"]
q.val<-sig_genes_exp.diff[,"q_value"]
names(siggeneqval)<-q.val
names(siggeneqval)<-ENTREZsiggenes
siggenepval<-sig_genes_exp.diff[,"p_value"]
names(siggenepval)<-ENTREZsiggenes

xx <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
topDiffGenes <- function(allScore) {
return(allScore < 0.01)}

noOVERbpGOdata <- new("topGOdata",ontology = "BP",allGenes = QvalnoOVER, nodeSize = 5,
			  annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
			  ID = "symbol")
noUNDERbpGOdata <- new("topGOdata",ontology = "BP",allGenes = QvalnoUNDER, nodeSize = 5,
			  annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
			  ID = "symbol")
under.upBPGOdata <- new("topGOdata",ontology = "BP",allGenes = QvalUNDERhi, nodeSize = 5,
			  annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
			  ID = "symbol")
over.upMFGOdata <- new("topGOdata",ontology = "BP",allGenes = QvalOVERhi, nodeSize = 5,
			  annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
			  ID = "symbol")
noOVERMFGOdata <- new("topGOdata",ontology = "MF",allGenes = QvalnoOVER, nodeSize = 5,
			  annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
			  ID = "symbol")
noUNDERMFGOdata <- new("topGOdata",ontology = "MF",allGenes = QvalnoUNDER, nodeSize = 5,
			  annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
			  ID = "symbol")
under.upCCGOdata <- new("topGOdata",ontology = "CC",allGenes = QvalUNDERhi, nodeSize = 5,
			  annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
			  ID = "symbol")
over.upCCGOdata <- new("topGOdata",ontology = "CC",allGenes = QvalOVERhi, nodeSize = 5,
			  annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
			  ID = "symbol")
noOVERCCGOdata <- new("topGOdata",ontology = "CC",allGenes = QvalnoOVER, nodeSize = 5,
			  annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
			  ID = "symbol")
noUNDERCCGOdata <- new("topGOdata",ontology = "CC",allGenes = QvalnoUNDER, nodeSize = 5,
			  annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
			  ID = "symbol")

noUNDERbptKS <- runTest(noOVERbpGOdata, algorithm = "classic", statistic = "ks")
noUNDERbpFisher <- runTest(noOVERbpGOdata, algorithm = "classic", statistic = "fisher")
noUNDERbptKS <- runTest(noUNDERbpGOdata, algorithm = "classic", statistic = "ks")
noUNDERbpweight01 <- runTest(noUNDERbpGOdata, algorithm = "classic", statistic = "fisher")
noUNDERbptKS.elim <- runTest(noUNDERbpGOdata, algorithm = "elim", statistic = "ks")

noUNDERcctKS <- runTest(noUNDERCCGOdata, algorithm = "classic", statistic = "ks")
noUNDERccFisher <- runTest(noUNDERCCGOdata, algorithm = "classic", statistic = "fisher")
noUNDERcctKS <- runTest(noUNDERCCGOdata, algorithm = "classic", statistic = "ks")
noUNDERccweight01 <- runTest(noUNDERCCGOdata, algorithm = "classic", statistic = "fisher")
noUNDERcctKS.elim <- runTest(noUNDERCCGOdata, algorithm = "elim", statistic = "ks")

noUNDERmftKS <- runTest(noUNDERMFGOdata, algorithm = "classic", statistic = "ks")
noUNDERmfFisher <- runTest(noUNDERMFGOdata, algorithm = "classic", statistic = "fisher")
noUNDERmftKS <- runTest(noUNDERMFGOdata, algorithm = "classic", statistic = "ks")
noUNDERmfweight01 <- runTest(noUNDERMFGOdata, algorithm = "classic", statistic = "fisher")
noUNDERmftKS.elim <- runTest(noUNDERMFGOdata, algorithm = "elim", statistic = "ks")

noOVERbptKS <- runTest(noOVERbpGOdata, algorithm = "classic", statistic = "ks")
noOVERbpFisher <- runTest(noOVERbpGOdata, algorithm = "classic", statistic = "fisher")
noOVERbptKS <- runTest(noOVERbpGOdata, algorithm = "classic", statistic = "ks")
noOVERbpweight01 <- runTest(noOVERbpGOdata, algorithm = "classic", statistic = "fisher")
noOVERbptKS.elim <- runTest(noOVERbpGOdata, algorithm = "elim", statistic = "ks")

noOVERcctKS <- runTest(noOVERCCGOdata, algorithm = "classic", statistic = "ks")
noOVERccFisher <- runTest(noOVERCCGOdata, algorithm = "classic", statistic = "fisher")
noOVERcctKS <- runTest(noOVERCCGOdata, algorithm = "classic", statistic = "ks")
noOVERccweight01 <- runTest(noOVERCCGOdata, algorithm = "classic", statistic = "fisher")
noOVERcctKS.elim <- runTest(noOVERCCGOdata, algorithm = "elim", statistic = "ks")

noOVERmftKS <- runTest(noOVERMFGOdata, algorithm = "classic", statistic = "ks")
noOVERmfFisher <- runTest(noOVERMFGOdata, algorithm = "classic", statistic = "fisher")
noOVERmftKS <- runTest(noOVERMFGOdata, algorithm = "classic", statistic = "ks")
noOVERmfweight01 <- runTest(noOVERMFGOdata, algorithm = "classic", statistic = "fisher")
noOVERmftKS.elim <- runTest(noOVERMFGOdata, algorithm = "elim", statistic = "ks")

over.upMFtKS <- runTest(over.upMFGOdata, algorithm = "classic", statistic = "ks")
over.upMFFisher <- runTest(over.upMFGOdata, algorithm = "classic", statistic = "fisher")
over.upMFtKS <- runTest(over.upMFGOdata, algorithm = "classic", statistic = "ks")
over.upMFweight01 <- runTest(over.upMFGOdata, algorithm = "classic", statistic = "fisher")
over.upMFtKS.elim <- runTest(over.upMFGOdata, algorithm = "elim", statistic = "ks")

over.upCCtKS <- runTest(over.upCCGOdata, algorithm = "classic", statistic = "ks")
over.upCCFisher <- runTest(over.upCCGOdata, algorithm = "classic", statistic = "fisher")
over.upCCtKS <- runTest(over.upCCGOdata, algorithm = "classic", statistic = "ks")
over.upCCweight01 <- runTest(over.upCCGOdata, algorithm = "classic", statistic = "fisher")
over.upCCtKS.elim <- runTest(over.upCCGOdata, algorithm = "elim", statistic = "ks")

under.upMFtKS <- runTest(under.upMFGOdata, algorithm = "classic", statistic = "ks")
under.upMFFisher <- runTest(under.upMFGOdata, algorithm = "classic", statistic = "fisher")
under.upMFtKS <- runTest(under.upMFGOdata, algorithm = "classic", statistic = "ks")
under.upMFweight01 <- runTest(under.upMFGOdata, algorithm = "classic", statistic = "fisher")
under.upMFtKS.elim <- runTest(under.upMFGOdata, algorithm = "elim", statistic = "ks")

under.upCCtKS <- runTest(under.upCCGOdata, algorithm = "classic", statistic = "ks")
under.upCCFisher <- runTest(under.upCCGOdata, algorithm = "classic", statistic = "fisher")
under.upCCtKS <- runTest(under.upCCGOdata, algorithm = "classic", statistic = "ks")
under.upCCweight01 <- runTest(under.upCCGOdata, algorithm = "classic", statistic = "fisher")
under.upCCtKS.elim <- runTest(under.upCCGOdata, algorithm = "elim", statistic = "ks")

under.upBPtKS <- runTest(under.upBPGOdata, algorithm = "classic", statistic = "ks")
under.upBPFisher <- runTest(under.upBPGOdata, algorithm = "classic", statistic = "fisher")
under.upBPtKS <- runTest(under.upBPGOdata, algorithm = "classic", statistic = "ks")
under.upBPweight01 <- runTest(under.upBPGOdata, algorithm = "classic", statistic = "fisher")
under.upBPtKS.elim <- runTest(under.upBPGOdata, algorithm = "elim", statistic = "ks")

pdf("GOplots_LUTS_vx_CTRL.pdf")

showSigOfNodes(under.upBPGOdata, score(under.upBPFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(under.upBPGOdata, score(under.upBPweight01), firstSigNodes = 5, useInfo = "def" )
printGraph(under.upBPGOdata, under.upBPFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(under.upBPGOdata, under.upBPweight01, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(over.upMFGOdata, score(over.upMFFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(over.upMFGOdata, score(over.upMFweight01), firstSigNodes = 5, useInfo = "def" )
printGraph(over.upMFGOdata, over.upMFFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(over.upMFGOdata, over.upMFweight01, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(noOVERMFGOdata, score(noOVERmfFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(noOVERMFGOdata, score(noOVERmfweight01), firstSigNodes = 5, useInfo = "def" )
printGraph(noOVERCCGOdata, noOVERccFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(noOVERCCGOdata, noOVERccweight01, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(noUNDERMFGOdata, score(noUNDERmfFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(noUNDERMFGOdata, score(noUNDERmfweight01), firstSigNodes = 5, useInfo = "def" )
printGraph(noUNDERCCGOdata, noUNDERccFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(noUNDERCCGOdata, noUNDERccweight01, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
showSigOfNodes(noUNDERbpGOdata, score(noUNDERbpFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(noUNDERbpGOdata, score(noUNDERbpweight01), firstSigNodes = 5, useInfo = "def" )
printGraph(noUNDERbpGOdata, noUNDERbpFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(noUNDERbpGOdata, noUNDERbpFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
showSigOfNodes(noOVERMFGOdata, score(noOVERmfFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(noOVERMFGOdata, score(noOVERmfweight01), firstSigNodes = 5, useInfo = "def" )
printGraph(over.upCCGOdata, over.upCCFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(noUNDERbpGOdata, noUNDERbpFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "def", pdfSW = TRUE)
printGraph(noOVERbpGOdata, noOVERbpweight01, firstSigNodes = 10, noOVERbpweight01, fn.prefix = "tGO", useInfo = "def")
printGraph(noOVERMFGOdata, noOVERmfweight01, firstSigNodes = 10, noOVERmfFisher, fn.prefix = "tGO", useInfo = "def")
printGraph(noOVERCCGOdata, noOVERcctKS.elim, firstSigNodes = 15, noOVERccFisher, fn.prefix = "tGO", useInfo = "all")		   

sig.tab <- GenTable(noOVERbpGOdata, Fis = noOVERbpFisher, topNodes = 20)
## results of both test
sig.tab.CC.under.up <- GenTable(under.upCCGOdata, under.upCCFisher, under.upCCtKS, topNodes = 20)
## results of both test with specified names
sig.tab.BP.noover.<- GenTable(noOVERbpGOdata, Fis = noOVERbpFisher, KS = noOVERbptKS, topNodes = 20)
## results of both test with specified names and specified ordering
sig.tab.MF.over.up <- GenTable(over.upMFGOdata, Fis = over.upMFFisher, KS = over.upMFtKS, orderBy = "KS", ranksOf = "Fis", topNodes = 20)

goID <- sig.tab.MF.over.up[1, "GO.ID"]
print(showGroupDensity(over.upMFGOdata, goID, ranks = TRUE))
goID <- sig.tab.BP.noover.[1, "GO.ID"]
print(showGroupDensity(noOVERbpGOdata, goID, ranks = TRUE))
goID <- sig.tab.CC.under.up[1, "GO.ID"]
print(showGroupDensity(under.upCCGOdata, goID, ranks = TRUE))

sigpaths<- enrichPathway(gene=names(siggeneqval),pvalueCutoff=0.05, readable=T)
enrichMap(sigpaths,fixed=T,n=6, vertex.label.cex =.7,vertex.label.font = 2)
UNDERhinopaths<- enrichPathway(gene=ENTREZQvalUNDER,pvalueCutoff=0.05, readable=T,organism="human")
enrichMap(UNDERhinopaths,fixed=T,n=6, vertex.label.cex =.7,vertex.label.font = 2)
OVERhinopaths<- enrichPathway(gene=ENTREZQvalOVER,pvalueCutoff=0.05, readable=T)
enrichMap(OVERhinopaths,fixed=T,n=6, vertex.label.cex =.7,vertex.label.font = 2)
UNDERnopaths<- enrichPathway(gene=ENTREZQvalOVERhi,pvalueCutoff=0.05, readable=T,organism="human")
enrichMap(UNDERnopaths,fixed=T,n=6, vertex.label.cex =.7,vertex.label.font = 2)
OVERnopaths<- enrichPathway(gene=ENTREZQvalUNDERhi,pvalueCutoff=0.05, readable=T)
enrichMap(OVERnopaths,fixed=T,n=6, vertex.label.cex =.7,vertex.label.font = 2)

dotplot(UNDERnopaths, showCategory=50,title="UNDER paths",font.size = 4)
dotplot(OVERnopaths, showCategory=50,title="OVER paths",font.size = 4)
dotplot(UNDERhinopaths, showCategory=50,title="UNDER paths",font.size = 4)

barplot(UNDERnopaths, drop=F, showCategory=50,title="UNDER paths",font.size = 3)
barplot(OVERnopaths,drop=F,showCategory = 50,title="OVER paths", font.size = 3)
barplot(UNDERhinopaths, drop=F, showCategory=50,title="UNDER highly expressed",font.size = 3)

enrichMap(sigpaths,fixed=T,n=6, vertex.label.cex =.7,vertex.label.font = 2)
cnetplot(sigpaths, categorySize="pvalue", foldChange=siggenelog)
cnetplot(UNDERnopaths, categorySize="pvalue", foldChange=siggenelog)
cnetplot(OVERnopaths, categorySize="pvalue", foldChange=siggenelog)
cnetplot(UNDERhinopaths, categorySize="pvalue", foldChange=siggenelog)

egoCC <- enrichGO(gene = names(siggeneqval),
                universe = ENTREZids, OrgDb = org.Hs.eg.db, ont = "CC",
                pAdjustMethod = "BH", pvalueCutoff  = 0.01,
				  qvalueCutoff  = 0.05,readable = TRUE)
egoMF <- enrichGO(gene= names(siggeneqval),
                universe  = ENTREZids , OrgDb  = org.Hs.eg.db,
                ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05, readable = TRUE)
egoBP <- enrichGO(gene = names(siggeneqval), universe = ENTREZids,
				  OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,readable = TRUE)
noUNDERegoBP <- enrichGO(gene = ENTREZQvalUNDER,
                universe  = ENTREZids  ,  #ENTREZids,,
                OrgDb = org.Hs.eg.db, ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05, readable = TRUE)
noOVERegoBP <- enrichGO(gene = ENTREZQvalOVER,
                universe  = ENTREZids  ,  #ENTREZids,,
                OrgDb = org.Hs.eg.db, ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05, readable = TRUE)
noUNDERegoCC <- enrichGO(gene = ENTREZQvalUNDER,
                universe  = ENTREZids  ,  #ENTREZids,,
                OrgDb = org.Hs.eg.db, ont = "CC",
                pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05, readable = TRUE)
noOVERegoCC <- enrichGO(gene = ENTREZQvalOVER,
                universe  = ENTREZids  ,  #ENTREZids,,
                OrgDb = org.Hs.eg.db, ont = "CC",
                pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05, readable = TRUE)
noUNDERegoMF <- enrichGO(gene = ENTREZQvalUNDER,
                universe  = ENTREZids  ,  #ENTREZids,,
                OrgDb = org.Hs.eg.db, ont = "MF",
                pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05, readable = TRUE)
noOVERegoMF <- enrichGO(gene = ENTREZQvalOVER,
                universe  = ENTREZids  ,  #ENTREZids,,
                OrgDb = org.Hs.eg.db, ont = "MF",
                pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05, readable = TRUE)

kegg_enrichsig<- enrichKEGG(gene = ENTREZsiggenes, organism = 'human')
head(summary(kegg_enrichsig))
kegg_enrichqover <- enrichKEGG(gene = ENTREZQvalOVER, organism = 'human')
head(summary(kegg_enrichqover))
kegg_enrichHI <- enrichKEGG(gene = ENTREZQvalUNDER, organism = 'human')
head(summary(kegg_enrichHI))

barplot(kegg_enrichsig,drop = TRUE, showCategory = 50,
		title = "KEGG Enrichment ALL Pathways", font.size = 4)
barplot(kegg_enrichqover,drop = TRUE, showCategory = 50,
		title = "KEGG Enrichment OVER expr Pathways", font.size = 4)
barplot(kegg_enrichHI,drop = TRUE, showCategory = 50,
		title = "KEGG Enrichment UNDER exprPathways", font.size = 4)

ccgo <- groupGO(gene = names(siggeneqval), OrgDb = org.Hs.eg.db,
               ont = "CC", level    = 3, readable = TRUE)
bpgo <- groupGO(gene = names(siggeneqval), OrgDb = org.Hs.eg.db, ont = "BP",
               level    = 3, readable = TRUE)
mfgo <- groupGO(gene = names(siggeneqval), OrgDb = org.Hs.eg.db,
               ont = "MF",level  = 3, readable = TRUE)

cnetplot(egoBP, categorySize="pvalue", foldChange=siggenelog)
cnetplot(egoMF, categorySize="pvalue", foldChange=siggenelog)
cnetplot(egoCC, categorySize="pvalue", foldChange=siggenelog)
cnetplot(noOVERegoBP, categorySize="pvalue", foldChange=siggenelog)
cnetplot(noUNDERegoCC, categorySize="pvalue", foldChange=siggenelog)
cnetplot(noOVERegoMF, categorySize="pvalue", foldChange=siggenelog)

plotGOgraph(egoBP, firstSigNodes = 50, useInfo = "all", sigForAll = TRUE, useFullNames = F)
plotGOgraph(egoMF, firstSigNodes = 50, useInfo = "all", sigForAll = TRUE, useFullNames = F)
plotGOgraph(egoCC, firstSigNodes = 50, useInfo = "all", sigForAll = TRUE, useFullNames = F)


barplot(ccgo, drop=TRUE, showCategory=50,title="Cellular Component",font.size = 4)
barplot(mfgo, drop = TRUE, showCategory=50,title="Molecular Function",font.size = 4)
barplot(bpgo,drop = TRUE, showCategory=50, title = "Biological Process",font.size = 4) 

dotplot(sigpaths, title = "Significantly DE paths",font.size = 4)
dotplot(egoBP,showCategory=50, title = "Biological Process",font.size = 4)
dotplot(egoMF,showCategory=50, title = "Molecular Function",font.size = 4)
dotplot(egoCC,showCategory=50, title = "Biological Process",font.size = 4)

enrichMap(egoBP,vertex.label.font=.1)
enrichMap(egoMF,vertex.label.font=.1)
enrichMap(egoCC,vertex.label.font=.1,)

cnetplot(kegg_enrich,,showCategory=5,categorySize="geneNum", foldChange=siggenelog)
cnetplot(ccgo,showCategory=5,categorySize="geneNum", foldChange=siggenelog)
cnetplot(mfgo,showCategory=5,categorySize="geneNum", foldChange=siggenelog)
cnetplot(bpgo,showCategory=5,categorySize="geneNum", foldChange=siggenelog)

dev.off()
