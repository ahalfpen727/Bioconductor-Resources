

setwd("/project/umb_triley/cpct/rna-seq/urine1/NetClass")

#load the Adjeacncy Matrix
load("Adj.rda")
#load the Gene Expression data
load("geneData.rda")
library(netClass)
data(EN2SY)
labley=c(1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1)

dk=calc.diffusionKernelp(L=kk,p=3,a=1)



r.stsvm <- cv.stsvm(x=ValuesofGenes,x.mi=NULL,y=labley,folds=5,Gsub=kk,op.method="pt",
                    repeats=1, parallel=FALSE, cores=2,DEBUG=TRUE,pt.pvalue=0.05,op=0.9,
                    aa=5,a=1,p=2,allF=TRUE, seed=1234,Cs=10^(-3:3))
save(r.stsvm,file="Gtf-NetClass.rda")