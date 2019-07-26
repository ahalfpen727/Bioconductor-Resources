library(cummeRbund)
library(DOSE);library(GO.db);library(org.Hs.eg.db)
library(topGO);library(GSEABase);library(clusterProfiler)
library(biomaRt);library(KEGG.db);library(reactome.db)
###############################################################################
# Features, Counts, and all inclusive tables 
# (CDS-p_ids, TSS-TSS ids, isoforms- XLOC and TCONS ids, and Genes XLOC ids)
# Significantly differentially expressed features: method 1
############################################################################
## ---CummeRbund and GSEA ------------------------------------------------------------------------------------------------
#source('~/R_Py_scripts/dcEnrichmentGO.R', echo=TRUE)
cuff<-readCufflinks(gtfFile="cuffcmp.combined.gtf",genome="genome.fa",rebuild=T)
cuff
refgtf="cuffcmp.combined.gtf"
genome_path="genome.fa"
over="LUTS"
under="CTRL"
###############################################################################
# Features, Counts, and all inclusive tables 
# get gene symbols for CDS TSS Isoforms features
################################################################################
runInfo(cuff);replicates(cuff);conditions(cuff)
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

g.rep.matrix<-repFpkmMatrix(isoforms(cuff))
gene.xloc.matrix<-featureNames(isoforms(cuff))
gene.list<-getGenes(cuff,geneId = gene.xloc.matrix)
gene.list
gene_annotation_data<-featureNames(isoforms(gene.list))
head(gene_annotation_data)
gene.rep.matrix<-cbind2(gene_annotation_data, g.rep.matrix)
head(gene.xloc.matrix)
length(gene.xloc.matrix)
head(gene.rep.matrix)
#####################################################
# Gene synbol mapping and preparation for go
#####################################################
glist<-as.list(as.character(gene.rep.matrix[,"gene_short_name"]))
length(glist)
dim(glist)
siggenelist<-sig_genes_exp.diff[,"gene_id"]

fpkmatrix<-fpkmMatrix(gene.list)
ctmtx<-countMatrix(isoforms(cuff))
names(ctmtx)
genematrix<-cbind2(glist, ctmtx)
head(genematrix)
isodiff<-diffTable(isoforms(cuff))

foldchange=sig_genes_exp.diff[,"log2_fold_change"]
Lutsup.ctrlnoexp = which(foldchange == "-Inf")
NOexp.inCTRL<-as.data.frame(sig_genes_exp.diff[Lutsup.ctrlnoexp,])
NOexp.inCTRL
ctrl.noexp<-length(which(foldchange == "-Inf"))
ctrl.noexp

lutsup = which(foldchange > 1)
LUTSup<-as.data.frame(sig_genes_exp.diff[lutsup,])
LUTSup
luts.upexp = length(which(foldchange > 1))
luts.upexp

lutsdown = which(foldchange < 1)
LUTS.downexp<-as.data.frame(sig_genes_exp.diff[lutsdown,])
LUTS.downexp
lutsdown = length(which(foldchange < 1, foldchange ))
lutsdown

counts=table(tnc_de[,"gene_name"])
order(sig_genes_exp.diff[,"log2_fold_change"], decreasing=T,)
which(foldchange == "-Inf"))

hist(Lutsup.ctrlnoexp, breaks=50, col="bisque4",) # xlab="Transcripts per gene", main="Distribution of isoform transcript count per gene")
legend_text = c(paste("Genes with one isoform =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
legend("bottomright", legend_text, lty=NULL)
counts=table(tnc_de[,"gene_name"])
c_one = length(which(counts == 1))
c_more_than_one = length(which(counts > 1))
c_max = max(counts)
order(sig_genes_exp.diff[,"log2_fold_change"], decreasing=T,)
