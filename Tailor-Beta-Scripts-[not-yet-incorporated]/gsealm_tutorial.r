
#acquire the GSEAlm package from bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite("made4")
library(made4)
#load the library
library(GSEAlm)
??GSEAlm
#load the provided gene expression dataset
load('breastCancer.RData')

#meta information in the gene expression dataset
breastCancer

#phenotype data contained in the expression dataset
pData(breastCancer)

#table of the phenotype data shows normal tissue and basal breast cancer
table(pData(breastCancer))

#load GO terms
load('mSigDB.RData')

#run GSEAlm
#nperm determines the number of permutations used to estimate the null distribution
#of the enrichment score
pVals=gsealmPerm(breastCancer,~disease_type,mSigDB,nperm=1000)
head(pVals)

#correction for multiple testing
pVals<-apply(pVals,2,p.adjust,method='BH',n=nrow(pVals))
head(pVals)

#set a significance threshold
THRESHOLD<-0.03

#gene sets that are upregulated in basal breast cancer
sort(pVals[pVals[,1]<THRESHOLD,1])[1:20]

#gene sets that are downregulated in the BLC
pVals[pVals[,2]<THRESHOLD,2]


