
#list of GO terms (biological processes) acquired from MSigDB
#http://www.broadinstitute.org/gsea/downloads.jsp
mSigDB<-readLines('c5.bp.v3.0.symbols.gmt')

#the gmt format has the following structure
#each line represents a gene set, which are seperated by a tab
#the first entry is the name
#the second one the web address to mSigDB
#and the rest is a list of genes, in gene symbols
mSigDB[1]

#split each line at the tabs
mSigDB<-strsplit(mSigDB,'\t')
mSigDB[1]

#use the first entry of each gene set as name
names(mSigDB)<-sapply(mSigDB,function(x)x[1])
mSigDB[1]

#remove the first two entries
mSigDB<-sapply(mSigDB,function(x)x[3:length(x)])
mSigDB[1]

#get a unique list of all genes
allGenes<-unique(unlist(mSigDB))

#use biomaRt to map the gene symbols to affmetrix u133 plus2 identifiers
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("biomaRt")

library("biomaRt")
ensembl=useMart("ensembl",dataset="hsapiens_gene_ensembl")
geneMaps<-getBM (attributes =  c('hgnc_symbol','affy_hg_u133_plus_2'),
                 filters='hgnc_symbol',
                 values=allGenes,
                 mart=ensembl)
head(geneMaps)

#remove empty entries
geneMaps<-geneMaps[!geneMaps[,2]=='',]

#vector with all affymetrix u133 plus2 probe ids from our breast cancer dataset
load('breastCancer.RData')
probeIds<-rownames(exprs(breastCancer))
binaryGeneSet<-rep(0,length(probeIds))
names(binaryGeneSet)<-probeIds

#function that maps the list of genes to a binary vector where 1 indicates that
#a gene is present in a particular gene set
mapGeneSet<-function(geneSet,geneMap,binaryGeneSet){
   probesOfGeneSet<-geneMap[geneMap[,1]%in%geneSet,2]
   binaryGeneSet[probesOfGeneSet]<-1
   return(binaryGeneSet)
}
mSigDB<-sapply(mSigDB,mapGeneSet,geneMaps,binaryGeneSet)
mSigDB<-t(mSigDB)

#save the processed gene sets for further use
save(mSigDB,file='mSigDB.RData')

