source("geneExpressionContribution-source.R")
nb.samp = 50

## Single gene, 2 isoforms, simple expression matrix
dat1 = matrix(rnorm(nb.samp*2,c(3,6)),2)

VlsVt(dat1)
VlsVt(dat1,rescale=FALSE)
VlsVt.bootstrap(dat1)
VlsVt.bootstrap(dat1,rescale=FALSE)
VlsVt.bootstrap(dat1,nb.bs=1000)


##
## Several genes into a data.frame
##

## Dummy data construction
dat2 = gTr =  NULL
for(nb.iso in 1:5){
  for(tr in 1:nb.iso){
    dat2 = rbind(dat2,abs(rnorm(nb.samp,runif(1,0,6))))
    gTr = rbind(gTr, data.frame(trId=paste("t",tr,sep=""),geneId=paste("g",nb.iso,sep="")))
  }
}
colnames(dat2) = paste("samp",1:nb.samp,sep="")
dat2 = cbind(gTr,dat2)
dat2[1:6,1:6] ## Have a look at the format

## Vls.Vt computation
VlsVt.comp(dat2)
VlsVt.comp(dat2,min.iso2.RPKM=2,rescale=TRUE,nb.bs=100,prop.bs=.5,verbose=TRUE)
VlsVt.comp(dat2,nb.bs=NULL)
