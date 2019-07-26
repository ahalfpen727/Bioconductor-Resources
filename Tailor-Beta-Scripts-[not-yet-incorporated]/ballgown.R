install.packages("devtools",repos="http://cran.us.r-project.org")
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("alyssafrazee/RSkittleBrewer")
library(cummeRbund)
###################################################################
dir("~GRCh38_gtf_only")
cuff<-readCufflinks(gtfFile="cuffcmp.combined.gtf",genome="hg38",rebuild=T)
cuff
refgtf="cuffcmp.combined.gtf"
genome_path="hg38"
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
length(sigGenes)

sig_genes_exp.diff<-diffData(sigGenes)
sig_genes_annot<-sigGenes@annotation
head(sig_genes_exp.diff)
head(sig_genes_annot)
sig_genes_exp.diff<-cbind2(sig_genes_annot["gene_short_name"], sig_genes_exp.diff)
head(sig_genes_exp.diff)
head(sig_genes_exp.diff)
colnames(sig_genes_exp.diff)

mySigCDS<-getSig(cuff,x=over,y=under,alpha=.05,level='CDS')
head(mySigCDS)
length(mySigCDS)
mySigCDSTable<-getSigTable(cuff,alpha=0.05,level='CDS')
head(mySigCDSTable,20)
length(mySigCDSTable)
sigCDS<-getGenes(cuff, mySigCDS)
sigCDS
length(sigCDS)

sig_cds_exp.diff<-diffData(sigCDS)
sig_cds_annot<-sigCDS@annotation
head(sig_cds_exp.diff)
head(sig_cds_annot)
sig_cds_exp.diff<-cbind2(sig_cds_annot["gene_short_name"], sig_cds_exp.diff)
head(sig_cds_exp.diff)
colnames(sig_cds_exp.diff)

mySigIsos<-getSig(cuff,x=over,y=under,alpha=.05,level='isoforms')
head(mySigIsos)
length(mySigIsos)
mySigIsoTable<-getSigTable(cuff,alpha=0.05,level='isoforms')
head(mySigIsoTable,20)
length(mySigIsoTable)
sigIsos<-getGenes(cuff, mySigIsos)
sigIsos
length(sigIsos)

sig_isos_exp.diff<-diffData(sigIsos)
sig_isos_annot<-sigIsos@annotation
head(sig_isos_exp.diff)
head(sig_isos_annot)
sig_isos_exp.diff<-cbind2(sig_isos_annot["gene_short_name"], sig_isos_exp.diff)
head(sig_isos_exp.diff)
colnames(sig_isos_exp.diff)

mySigTSS<-getSig(cuff,x=over,y=under,alpha=.05,level='TSS')
head(mySigTSS)
length(mySigTSS)
mySigTSSTable<-getSigTable(cuff,alpha=0.05,level='TSS')
head(mySigTSSTable,20)
length(mySigTSSTable)
sigTSS<-getGenes(cuff, mySigTSS)
sigTSS
length(sigTSS)

sig_tss_exp.diff<-diffData(sigTSS)
sig_tss_annot<-sigTSS@annotation
head(sig_tss_exp.diff)
head(sig_tss_annot)
sig_tss_exp.diff<-cbind2(sig_tss_annot["gene_short_name"], sig_tss_exp.diff)
head(sig_tss_exp.diff)
colnames(sig_tss_exp.diff)

runInfo(cuff)
idpatient<-replicates(cuff)
colnames(idpatient)
g.rep.matrix<-repFpkmMatrix(genes(cuff))

cufftable<-diffTable(genes(cuff))
head(cufftable)

gene.xloc.matrix<-featureNames(genes(cuff))
gene.list<-getGeneId(cuff,gene.xloc.matrix)
gene.list
gene_annotation_data<-featureNames(gene.list)
head(gene_annotation_data)
gene.rep.matrix<-cbind2(gene_annotation_data, g.rep.matrix)

phenotype=idpatient[,"sample_name"]
ids=idpatient[,"file"]
replicate=idpatient[,"replicate"]
pheno<-cbind2(replicate, ids)
head(pheno)
pheno<-cbind2(pheno, phenotype)
head(pheno)
pheno<-as.data.frame(pheno, colnames=T, header=T, sep="\t")
head(pheno)

######################################################

g.rep.matrix<-repFpkmMatrix(genes(cuff))
head(g.rep.matrix)
colnames(g.rep.matrix)<-idpatient[,"file"]
gene.xloc.matrix<-featureNames(genes(cuff))
gene_annotation_data<-featureNames(gene.list)
head(gene_annotation_data)



library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

bg = ballgown(dataDir=".", samplePattern='')
pData(bg) = data.frame(id=sampleNames(bg), group=rep(c(1,0), each=10))
dataDir = dir

# gtf = file.path(dataDir, 'hg19_genes_small.gtf.gz')
gtf=refgtf
rsemobj = ballgownrsem(dir=dataDir, samples=c('tiny', 'tiny2'), gtf=gtf,
    bamout='none', zipped=TRUE)
rsemobj
ballgownrsem(dir = "", samples, gtf, UCSC = TRUE,
  tfield = "transcript_id", attrsep = "; ", bamout = "transcript",
  pData = NULL, verbose = TRUE, meas = "all", zipped = FALSE)

dataDir = system.file('extdata', package='ballgown')

rsemobj = ballgownrsem(dir=dataDir, samples=c('tiny', 'tiny2'), gtf=gtf,
    bamout='none', zipped=TRUE)
rsemobj
bg_chrX = ballgown(samplePattern = "LUTS", pData=pheno)
# filter out genes with zero-ish counts by removing transcripts
# with a variance less than one
bg_chrX_filt = subset(bg_chrX,
					  "rowVars(texpr(bg_chrX)) >1",
					  genomesubset=TRUE)
# identify DE transcripts using glm
# getFC=TRUE parameter confounder-adjusted fold change 
# between the two groups 
# in this case the groups under observation are m/f (sex)
results_transcripts=stattest(bg_chrX_filt,
							   feature="transcript",
							   covariate="sex",
							   adjustvars=c("population"),
							   getFC=TRUE, meas="FPKM")
# do the same statistical test for gene level features now
results_genes = stattest(bg_chrX_filt, feature="gene",
covariate="sex", adjustvars = c("population"), getFC=TRUE,
meas="FPKM")

# Add gene names + IDs to the results_transcripts data frame

results_transcripts=data.frame(geneNames=ballgown::geneNames(bg_chrX_filt),
							   geneIDs=ballgown::geneIDs(bg_chrX_filt),
							   results_transcripts)
# Sort by smallest P value to the largest:
results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)
# write results to csv
write.csv(results_transcripts, "chrX_transcript_results.csv",
		  row.names=FALSE)
write.csv(results_genes, "chrX_gene_results.csv",
		  row.names=FALSE)
# Identify transcripts and genes with a q value <0.05:
subset(results_transcripts,results_transcripts$qval<0.05)
subset(results_genes,results_genes$qval<0.05)
# plots
tropical= c('darkorange', 'dodgerblue','hotpink', 
			'limegreen', 'yellow')
palette(tropical)
# plot dist with log2 transform 
# add one to prevent undefined terms
fpkm = texpr(bg_chrX,meas="FPKM")
fpkm = log2(fpkm+1)
boxplot(fpkm,col=as.numeric(pheno_data$sex),
		las=2,ylab='log2(FPKM+1)')

#Make plots of individual transcripts across samples
ballgown::transcriptNames(bg_chrX)[12]
ballgown::geneNames(bg_chrX)[12]

plot(fpkm[12,] ~ pheno_data$sex, border=c(1,2),
main=paste(ballgown::geneNames(bg_chrX)[12],' : ',
ballgown::transcriptNames(bg_chrX)[12]),pch=19, xlab="Sex",
ylab='log2(FPKM+1)')

points(fpkm[12,] ~ jitter(as.numeric(pheno_data$sex)),
col=as.numeric(pheno_data$sex))

plotTranscripts(ballgown::geneIDs(bg_
chrX)[1729], bg_chrX, main=c('Gene XIST in
sample ERR188234'), sample=c('ERR188234'))

plotMeans('MSTRG.56', bg_chrX_filt,groupvar="sex",legend=FALSE)








??structure()
data(bg)
gtfPath = system.file('extdata', 'annot.gtf.gz', package='ballgown')
geneoverlaps = getGenes(gtfPath, structure(bg)$trans, UCSC=FALSE)

gene.list<-getGenes(gtf = refgtf, UCSC = T,attribute = "gene_name")
gene.list

gtfPath = system.file('extdata', 'annot.gtf.gz', package='ballgown')
geneoverlaps = getGenes(gtfPath, structure(bg)$trans, UCSC=FALSE)

gtf = file.path(dataDir, 'hg19_genes_small.gtf.gz')






