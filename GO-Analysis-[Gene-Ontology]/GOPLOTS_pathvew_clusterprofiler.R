###############################################################################
# GOplots for Significantly differentially expressed features: method 1
############################################################################
# Libraries
###########################################################################
#source("http://bioconductor.org/biocLite.R")
#biocLite("topGO");biocLite("org.Hs.eg.db");biocLite("ReactomePA")
#biocLite("pathview");biocLite("reactome.db");biocLite("clusterProfiler")
#biocLite("GSEABase");biocLite("DOSE")
#biocLite("GO.db");biocLite("KEGG.db")
## ----echo=FALSE, results='hide', message=FALSE---------------------------
library(DOSE)
library(GO.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)
library(cummeRbund);library(biomaRt)
library(reactome.db);library(ReactomePA)
library(org.Mm.eg.db);library(org.Hs.eg.db)
library(KEGG.db);library(pathview)
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
refgtf = Sys.getenv("CUMMERBUND_INPUT_GTF")
genome_path = Sys.getenv("GENOME_PATH")
gen_v = Sys.getenv("R_GENOME")
cuff<-readCufflinks(gtfFile=refgtf,genome=genome_path,rebuild=T)
cuff
###################################################
### for testing
###################################################
#cuff<-readCufflinks(gtfFile="GRCh38_genes.gtf",genome="genome.fa",rebuild=F)
cuff<-readCufflinks(gtfFile="cuffcmp.combined.gtf",genome="genome.fa",rebuild=F)
cuff
refgtf="cuffcmp.combined.gtf"
genome_path="genome.fa"
over="LUTS"
under="CTRL"
###########################################################################
### Chunk_2: Significantly differentially expressed features: method 1
############################################################################
mySigGenes<-getSig(cuff,x=over,y=under,alpha=.05,level='genes')
head(mySigGenes)
length(mySigGenes)
mySigGeneTable<-getSigTable(cuff,alpha=0.05,level='genes')
head(mySigGeneTable,20)
length(mySigGeneTable)
sigGenes<-getGenes(cuff, mySigGenes)
head(sigGenes)
length(sigGenes)

sig_genes_exp.diff<-diffData(sigGenes)
sig_genes_annot<-sigGenes@annotation
head(sig_genes_exp.diff)
head(sig_genes_annot)
sig_genes_exp.diff["gene_id"]<-sig_genes_annot["gene_short_name"]
head(sig_genes_exp.diff)


mySigCDS<-getSig(cuff,x=over,y=under,alpha=.05,level='CDS')
head(mySigCDS)
length(mySigCDS)

sigCDS<-getGenes(cuff, mySigCDS)
sigCDS
length(sigCDS)

sig_cds_exp.diff<-diffData(sigCDS)
sig_cds_annot<-sigCDS@annotation
head(sig_cds_exp.diff)
head(sig_cds_annot)
sig_cds_exp.diff["gene_id"]<-sig_cds_annot["gene_short_name"]
head(sig_cds_exp.diff)


mySigIsos<-getSig(cuff,x=over,y=under,alpha=.05,level='isoforms')
head(mySigIsos)
length(mySigIsos)

sigIsos<-getGenes(cuff, mySigIsos)
sigIsos
length(sigIsos)

sig_isos_exp.diff<-diffData(sigIsos)
sig_isos_annot<-sigIsos@annotation
head(sig_isos_exp.diff)
head(sig_isos_annot)
sig_isos_exp.diff["gene_id"]<-sig_isos_annot["gene_short_name"]
head(sig_isos_exp.diff)

mySigTSS<-getSig(cuff,x=over,y=under,alpha=.05,level='TSS')
head(mySigTSS)
length(mySigTSS)

sigTSS<-getGenes(cuff, mySigTSS)
sigTSS
length(sigTSS)

sig_tss_exp.diff<-diffData(sigTSS)
sig_tss_annot<-sigTSS@annotation
head(sig_tss_exp.diff)
head(sig_tss_annot)
sig_tss_exp.diff["gene_id"]<-sig_tss_annot["gene_short_name"]
head(sig_tss_exp.diff)
dim(sig_tss_exp.diff)

###################################################
### Chunk_3: Differentially Loaded Features
###################################################

mySigrelCDS<-getSig(cuff,x=over,y=under,alpha=.05,level='relCDS')
head(mySigrelCDS)
length(mySigrelCDS)

relcds.diff<-distValues(relCDS(cuff))
relcds.sigdiff<-subset(relcds.diff, significant=="yes")
head(relcds.sigdiff)
dim(relcds.sigdiff)

mySigSplices<-getSig(cuff,x=over,y=under,alpha=.05,level='splicing')
head(mySigSplices)
length(mySigSplices)
splice.diff<-distValues(splicing(cuff))
splice.sigdiff<-subset(splice.diff, significant=="yes")
head(splice.sigdiff)
dim(splice.sigdiff)

mySigPromoters<-getSig(cuff,x=over,y=under,alpha=.05,level='promoters')
head(mySigPromoters)
length(mySigPromoters)

promoter.diff<-distValues(cuff@promoters)
promoter.sigdiff<-subset(promoter.diff, significant=="yes")
head(promoter.sigdiff)
dim(promoter.sigdiff)

###################################################
### Chunk_4: Sig Matrix_plots
###################################################
mySigMat<-sigMatrix(cuff,level='genes',alpha=0.05)
print(mySigMat)
sigIsoformIds1<-sigMatrix(cuff, level='isoforms',alpha=0.05)
print(sigIsoformIds1)
sigIsoformIds2<-sigMatrix(cuff, level='TSS',alpha=0.05)
print(sigIsoformIds2)
sigIsoformIds3<-sigMatrix(cuff, level='CDS', alpha=0.05)
print(sigIsoformIds3)

###############################################################################
# Chunk_5: Features, Counts, and all inclusive tables 
# (CDS-p_ids, TSS-TSS ids, isoforms- XLOC and TCONS ids, and Genes XLOC ids)
################################################################################
runInfo(cuff)
replicates(cuff)
conditions(cuff)

isodiff<-diffTable(isoforms(cuff))
head(isodiff)

###################################################################
### Chunk_7: get gene symbols for CDS TSS Isoforms and Genes
#####################################################################

gene.xloc.matrix<-featureNames(genes(cuff))
head(gene.xloc.matrix)
gene.list<-getGenes(cuff,geneId = gene.xloc.matrix)
gene.list
gene_annotation_data<-featureNames(gene.list)
head(gene_annotation_data)

isos.tcons.matrix<-featureNames(isoforms(cuff))
head(isos.tcons.matrix)
iso.list<-getGenes(cuff,geneId = isos.tcons.matrix)
iso.list
iso_annotation_data<-featureNames(iso.list)
head(iso_annotation_data)

cds.p_id.matrix<-featureNames(CDS(cuff))
head(cds.p_id.matrix)
cds.list<-getGenes(cuff,geneId = cds.p_id.matrix)
cds.list
cds_annotation_data<-featureNames(cds.list)
head(cds_annotation_data)

tss.tssgroup.matrix<-featureNames(TSS(cuff))
head(tss.tssgroup.matrix)
tss.list<-getGenes(cuff,geneId = tss.tssgroup.matrix)
tss.list
tss_annotation_data<-featureNames(tss.list)
head(tss_annotation_data)

#####################################################
### code chunk number 8: Write repFPKMMatrix file
#####################################################
g.rep.matrix<-repFpkmMatrix(isoforms(cuff))
gene.xloc.matrix<-featureNames(isoforms(cuff))
gene.list<-getGenes(cuff,geneId = gene.xloc.matrix)
gene.list
gene_annotation_data<-featureNames(gene.list)
head(gene_annotation_data)
gene.rep.matrix<-cbind2(g.rep.matrix, gene.xloc.matrix)
gene_table = file.path("./RepFpkmMatrix.genes.txt")
write.table(gene.rep.matrix,gene_table, sep="  ", row.names = T , col.names = T,quote = F)
head(gene.rep.matrix)
head(gene_annotation_data)
length(gene_annotation_data)
dim(gene_annotation_data)


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
colnames(isodiff)<-c("gene_id","class_code", "nearest_ref_id","gene_short_name",
                            "locus","length", "coverage","seqnames", "start","end",
                            "width" ,"strand", "source","type","score","phase", "isoform_id",
                            "exon_number","gene_name","oId","nearest_ref","TSS_group_id",
                            "CDS_id","contained_in", "LUTS_CTRL_status","LUTS_CTRL_value_1",
                            "LUTS_CTRL_value_2","LUTS_CTRL_log2_fold_change","LUTS_CTRL_test_stat",
                            "LUTS_CTRL_p_value", "LUTS_CTRL_q_value","LUTS_CTRL_significant")

colnames(isodiff)<-c("gene_id","class_code", "nearest_ref_id","gene_short_name",
                            "locus","length", "coverage","seqnames", "start","end",
                            "width" ,"strand", "source","type","score","phase", "isoform_id",
                            "exon_number","gene_name","oId","nearest_ref","TSS_group_id",
                            "CDS_id","contained_in", "status","value_1","value_2",
                            "log2_fold_change","test_stat", "p_value", "q_value","significant")

cuffde<-as.data.frame(genes_exp.diff, col.names=T, rownames=T,header=T, sep="\t", as.is=T)
head(cuffde)
## ------------------------------------------------------------------------
library(ReactomePA)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
## ------------------------------------------------------------------------
ids <- bitr(isodiff$gene_short_name, fromType="SYMBOL", toType=c("ENTREZID","GOALL","ONTOLOGYALL","PATH"), annoDb="org.Hs.eg.db")
head(ids)

Genes_exp.diff<-as.data.frame(genes_exp.diff, rownames=T, col.names=T, 
                              header=T,sep="\t",na.rm=T)

repdata<-replicates(cuff)
reps<-repdata$rep_name
samples<-repdata$sample_name

g.rep.matrix<-repFpkmMatrix(isoforms(cuff))
head(g.rep.matrix)
colnames(g.rep.matrix)<-samples

groups<-factor(g.rep.matrix,levels=0:1,labels=c(over, under))
g.over.rep<-g.rep.matrix[samples == over]
g.under.rep<-g.rep.matrix[samples == under]

groups.de<-factor(cuffde,levels=0:1,labels=c(over, under))
g.over.matrix<-cuffde[,"sample_1" == over]
g.under.matrix<-cuffde[,"sample_2" == under]

range(cuffde$log2_fold_change)

colnames(g.rep.matrix)<-reps

Cuff.overde<-cbind(geneID=genes_exp.diff$isoform_id, g.over.matrix)
Cuff.underde<-cbind(geneID=genes_exp.diff$isoform_id, g.under.matrix)

gene.feats<-gene.list@annotation
gene.diff<-gene.list@diff

colnames(gene.feats)
colnames(gene.diff)
dim(gene.feats)
dim(gene.diff)

Cuff.overDE<-cbind(gene_symbol=gene.feats$gene_short_name, geneID=gene.diff$gene_id)
Cuff.underDE<-cbind(gene_symbol=gene.feats$gene_short_name, geneID=gene.diff$gene_id)

head(Cuff.overDE)
head(Cuff.underDE)


genes.under.de<-Cuff.underDE[unique(gene.feats[,"gene_short_name"])]
genes.over.de<-Cuff.overDE[unique(gene.feats[,"gene_short_name"])]

foldchange<-sig_genes_exp.diff[,"log2_fold_change"]
under.inf<-which(foldchange == "-Inf")
over.inf = which(foldchange == "Inf")
over.= which(foldchange > 0, foldchange == "Inf", foldchange != "-Inf")
under. = which(foldchange < 0,  foldchange == "-Inf")

NOexp.inUNDER<-as.data.frame(sig_genes_exp.diff[over.inf,])
NOexp.inOVER<-as.data.frame(sig_genes_exp.diff[under.inf,])

LOexp.inUNDER<-as.data.frame(sig_genes_exp.diff[over.,])
LOexp.inOVER<-as.data.frame(sig_genes_exp.diff[under.,])

################################################################
##  Alternative script for directional expression
################################################################
mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"

################################################################

formula_res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichKEGG")

head(as.data.frame(formula_res))
dotplot(ck)
dotplot(formula_res)
dotplot(formula_res, x=~group) + ggplot2::facet_grid(~othergroup)

################################################################



QvalnoOVER<-NOexp.inOVER[,"q_value"]
names(QvalnoOVER)<-NOexp.inOVER[,"gene_id"]

QvalnoUNDER<-NOexp.inUNDER[,"q_value"]
names(QvalnoUNDER)<-NOexp.inUNDER[,"gene_id"]

OVERupexp<-as.data.frame(sig_genes_exp.diff[over.,])
QvalOVERhi<-OVERupexp[,"q_value"]
names(QvalOVERhi)<-OVERupexp[,"gene_id"]

UNDERupexp<-as.data.frame(sig_genes_exp.diff[under.,])
QvalUNDERhi<-UNDERupexp[,"q_value"]
names(QvalUNDERhi)<-UNDERupexp[,"gene_id"]

#QvalOVERhi #QvalUNDERhi
ENTREZQvalexp<-mapIds(x = org.Hs.eg.db,
                          keys = gene.list@annotation$gene_short_name,
                          column = "ENTREZID",
                          keytype = "SYMBOL",
			                    multiVals="first")

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
                          keys = gene.list@annotation$gene_short_name,
                          column = "ENTREZID",
                          keytype = "SYMBOL",
			                    multiVals="first")

under.hilog<-Genes_exp.diff[,"log2_fold_change"]
names(under.hilog)<-ENTREZsiggenes
under.hiqval<-Genes_exp.diff[,"q_value"]
names(under.hiqval)<-ENTREZQvalUNDERhi
under.hipval<-(Genes_exp.diff[,"p_value"])
names(under.hipval) <- ENTREZQvalUNDERhi

over.hilog<-c(Genes_exp.diff[,"log2_fold_change"])
names(over.hilog)<-(ENTREZQvalOVERhi)
over.hiqval<-Genes_exp.diff[,"q_value"]
names(over.hiqval)<-ENTREZQvalOVERhi
over.hipval<-Genes_exp.diff[,"p_value"]
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
                          keys = gene.list@annotation$gene_short_name,
                          column = "ENTREZID",
                          keytype = "SYMBOL",
					          multiVals="first") 
entrezIDS<-mapIds(x = org.Hs.eg.db,
                  keys = gene.list@annotation$gene_short_name,
                  column = "ENTREZID",
                  keytype = "SYMBOL",
                  multiVals="first") 

sig_Genes_exp<-cbind(ENTREZsiggenes,log2fold=gene.diff[,"log2_fold_change"])
sig_genes_exp<-cbind(ENTREZsiggenes,log2fold=gene.diff[,"log2_fold_change"])
head(sig_genes_exp)
dim(sig_genes_exp)
names(sig_genes_exp) = rownames(sig_genes_exp)

siggenelog<-Genes_exp.diff[,"log2_fold_change"]
names(siggenelog)<-ENTREZsiggenes
siggeneqval<-Genes_exp.diff[,"q_value"]
q.val<-Genes_exp.diff[,"q_value"]
names(siggeneqval)<-q.val
names(siggeneqval)<-ENTREZsiggenes
siggenepval<-Genes_exp.diff[,"p_value"]
names(siggenepval)<-ENTREZsiggenes

xx <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
topDiffGenes <- function(allScore) {
return(allScore < 0.01)}

# Significantly Differentially Expressed by the Under group
# i.e Luts over Ctrl --> Ctrl = Under group, Under group highly 
# expresses and No expression in the Over group

noOVERbpGOdata <- new("topGOdata",ontology = "BP",allGenes = QvalnoOVER, nodeSize = 5,
	       	  			     annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
					     	       ID = "symbol")
under.upBPGOdata <- new("topGOdata",ontology = "BP",allGenes = QvalUNDERhi, nodeSize = 5,
		    			       annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
					       	         ID = "symbol")
noOVERmfGOdata <- new("topGOdata",ontology = "MF",allGenes = QvalnoOVER, nodeSize = 5,
	       	  			     annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
					     	       ID = "symbol")
under.upMFGOdata <- new("topGOdata",ontology = "MF",allGenes = QvalUNDERhi, nodeSize = 5,
                        annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                        ID = "symbol")
noOVERccGOdata <- new("topGOdata",ontology = "CC",allGenes = QvalnoOVER, nodeSize = 5,
                      annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                      ID = "symbol")
under.upCCGOdata <- new("topGOdata",ontology = "CC",allGenes = QvalUNDERhi, nodeSize = 5,
                        annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                        ID = "symbol")


under.upBPtKS <- runTest(under.upBPGOdata, algorithm = "classic", statistic = "ks")
under.upBPFisher <- runTest(under.upBPGOdata, algorithm = "classic", statistic = "fisher")
under.upBPtKS.elim <- runTest(under.upBPGOdata, algorithm = "elim", statistic = "ks")

NOover.BPtKS <- runTest(noOVERbpGOdata, algorithm = "classic", statistic = "ks")
NOover.BPFisher <- runTest(noOVERbpGOdata, algorithm = "classic", statistic = "fisher")
NOover.BPtKS.elim <- runTest(noOVERbpGOdata, algorithm = "elim", statistic = "ks")

NOover.MFtKS <- runTest(noOVERmfGOdata, algorithm = "classic", statistic = "ks")
NOover.MFFisher <- runTest(noOVERmfGOdata, algorithm = "classic", statistic = "fisher")
NOover.MFtKS.elim <- runTest(noOVERmfGOdata, algorithm = "elim", statistic = "ks")

NOover.CCtKS <- runTest(noOVERccGOdata, algorithm = "classic", statistic = "ks")
NOover.CCFisher <- runTest(noOVERccGOdata, algorithm = "classic", statistic = "fisher")
NOover.CCtKS.elim <- runTest(noOVERccGOdata, algorithm = "elim", statistic = "ks")

underUP.MFtKS <- runTest(under.upMFGOdata, algorithm = "classic", statistic = "ks")
underUP.MFFisher <- runTest(under.upMFGOdata, algorithm = "classic", statistic = "fisher")
underUP.MFtKS.elim <- runTest(under.upMFGOdata, algorithm = "elim", statistic = "ks")

underUP.CCtKS <- runTest(under.upCCGOdata, algorithm = "classic", statistic = "ks")
underUP.CCFisher <- runTest(under.upCCGOdata, algorithm = "classic", statistic = "fisher")
underUP.CCtKS.elim <- runTest(under.upCCGOdata, algorithm = "elim", statistic = "ks")

pdf("GOplots_OVER_^.pdf")

showSigOfNodes(under.upMFGOdata, score(underUP.MFFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(under.upMFGOdata, score(underUP.MFtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(under.upMFGOdata, score(underUP.MFtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(under.upMFGOdata, underUP.MFFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(under.upMFGOdata, underUP.MFtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(under.upMFGOdata, underUP.MFtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(under.upCCGOdata, score(underUP.CCFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(under.upCCGOdata, score(underUP.CCtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(under.upCCGOdata, score(underUP.CCtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(under.upCCGOdata, underUP.CCFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(under.upCCGOdata, underUP.CCtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(under.upCCGOdata, underUP.CCtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(under.upBPGOdata, score(under.upBPFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(under.upBPGOdata, score(under.upBPtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(under.upBPGOdata, score(under.upBPtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(under.upBPGOdata, under.upBPFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(under.upBPGOdata, under.upBPtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(under.upBPGOdata, under.upBPtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(noOVERccGOdata, score(NOover.CCFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(noOVERccGOdata, score(NOover.CCtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(noOVERccGOdata, score(NOover.CCtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(noOVERccGOdata, NOover.CCFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(noOVERccGOdata, NOover.CCtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(noOVERccGOdata, NOover.CCtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(noOVERbpGOdata, score(NOover.BPFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(noOVERbpGOdata, score(NOover.BPtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(noOVERbpGOdata, score(NOover.BPtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(noOVERbpGOdata, NOover.BPFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(noOVERbpGOdata, NOover.BPtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(noOVERbpGOdata, NOover.BPtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(noOVERmfGOdata, score(NOover.MFFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(noOVERmfGOdata, score(NOover.MFtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(noOVERmfGOdata, score(NOover.MFtKS.elim), firstSigNodes = 5, useInfo = "def" )


ggo <- groupGO(gene     = ids,
               #OrgDb    = org.Hs.eg.db,
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
  ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
                  OrgDb         = org.Hs.eg.db,
              		keytype       = 'ENSEMBL',
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)

  ego2 <- setReadable(ego2, OrgDb = org.Hs.eg.db)

    ego3 <- gseGO(geneList     = geneList,
               OrgDb        = org.Hs.eg.db,
                ont          = "CC",
                nPerm        = 1000,
                minGSSize    = 100,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                verbose      = FALSE)

## ------------------------------------------------------------------------
search_kegg_organism('hsa', by='kegg_code')
Hsapien <- search_kegg_organism('Homo sapien', by='scientific_name')
dim(Hsapien)
head(Hsapien)

hsa.gene <- enrichKEGG(gene = gene, organism = 'hsa',
                       pvalueCutoff = 0.05)
head(hsa.gene)

kk2 <- gseKEGG(geneList = geneList, organism = 'hsa',
               nPerm = 1000, minGSSizev = 120,
               pvalueCutoff = 0.05,verbose  = FALSE)
head(kk2)

## ----eval = FALSE--------------------------------------------------------
 mkk <- enrichMKEGG(gene = gene, organism = 'hsa')

  mkk2 <- gseMKEGG(geneList = geneList, species = 'hsa')
## ----eval = FALSE--------------------------------------------------------

david <- enrichDAVID(gene = gene,idType = "ENTREZ_GENE_ID",
                     listType = "Gene",
                     annotation = "KEGG_PATHWAY",
                     david.user = "clusterProfiler@hku.hk")

## ------------------------------------------------------------------------
gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt",
                       package="clusterProfiler")
  
c5 <- read.gmt(gmtfile)

egmt <- enricher(gene, TERM2GENE=c5)
head(egmt)

egmt2 <- GSEA(geneList, TERM2GENE=c5, verbose=FALSE)
head(egmt2)

## ----fig.height=5, fig.width=9-------------------------------------------
barplot(ggo, drop=TRUE, showCategory=12)
barplot(ego, showCategory=8)
dotplot(ego)

## ----fig.cap="enrichment map of enrichment result", fig.align="center", fig.height=16, fig.width=16, eval=FALSE----
  enrichMap(ego)

# categorySize can be scaled by 'pvalue' or 'geneNum'
 cnetplot(ego, categorySize="pvalue", foldChange=geneList)

## ----fig.height=12, fig.width=8------------------------------------------
plotGOgraph(ego)

## ----fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8----
gseaplot(kk2, geneSetID = "hsa04145")

## ----eval=FALSE----------------------------------------------------------
  browseKEGG(kk, 'hsa04110')

## ----eval=FALSE----------------------------------------------------------
  library("pathview")
  hsa04110 <- pathview(gene.data  = geneList,hway.id = "hsa04110",
                       species = "hsa",limit = list(gene=max(abs(geneList)),
                                                    cpd=1))
data(gcSample)
lapply(gcSample, head)
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
dotplot(ck)
dotplot(formula_res)
dotplot(formula_res, x=~group) + ggplot2::facet_grid(~othergroup)



