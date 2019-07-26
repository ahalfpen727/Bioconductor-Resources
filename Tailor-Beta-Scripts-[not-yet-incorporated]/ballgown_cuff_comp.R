source("http://bioconductor.org/biocLite.R")
biocLite("mygene")
biocLite("pathfindr")
biocLite("spliceR")
vignette("pathfindr_vignette"))
library("spliceR")
library(devtools)
library(ballgown)
library(reshape2)
library(limma)
library(usefulstuff)
#install_github('alyssafrazee/usefulstuff')
library(EBSeq)

#Ballgown/cuffdiff data comparison script

transcripts = read.table("./isoform_exp.diff", header=TRUE)
rdata = read.table("./isoforms.read_group_tracking", header=TRUE)
fpkmtable = dcast(rdata, formula=tracking_id~replicate+condition, 
                  value.var='FPKM')
hiexpr = which(rowMeans(fpkmtable[,-1]) > 1)
length(hiexpr) 

hiID = fpkmtable$tracking_id[hiexpr]
hitranscripts = subset(transcripts, test_id %in% hiID)
hist(hitranscripts$p_value[hitranscripts$status == "OK"], xlab='p-values, highly-expressed transcripts, status "OK"', main='Cuffdiff 2.2.1 p-values, transcripts', col='gray')

genes = read.table("./gene_exp.diff", header=TRUE)
hist(genes$p_value[genes$status=="OK"], main='Cuffdiff 2.2.1 p-values, genes', xlab='p-values, status "OK"', col='gray') 

bgtable = fpkmtable[,-1]
rownames(bgtable) = fpkmtable[,1]
system.time(bgresults <- stattest(gowntable=bgtable, 
                                  pData=data.frame(group=as.numeric(grepl('LUTS', names(bgtable)))),
                                  feature='transcript', covariate='group'))
hist(bgresults$pval[bgresults$id %in% hiID], main='p-values, transcripts', xlab='p-values, highly-expressed transcripts', col='gray')


Data = acast(rdata, formula=tracking_id~replicate+condition, value.var='raw_frags')
Conditions = rep(c('CTRL', 'LUTS'), 12)
IsoformNames = rownames(Data)
iso_gene_relationship = read.table("./isoform_exp.diff", 
                                   colClasses=c('character', 'character', rep('NULL', 12)), header=TRUE)
sum(IsoformNames != iso_gene_relationship$test_id) # expect 0

IsosGeneNames = iso_gene_relationship$gene_id
IsoSizes = MedianNorm(Data)
NgList = GetNg(IsoformNames, IsosGeneNames)
IsoNgTrun = NgList$IsoformNgTrun

xy<-system.time(IsoEBOut <- EBTest(Data=Data, NgVector=IsoNgTrun, 
                               Conditions=as.factor(samples), sizeFactors=IsoSizes, maxround=10))

bhist(bgresults$pval[bgresults$id %in% hiID], fill='dodgerblue', alpha=0.6,
      xlab="p-values", ylab='Frequency', main="Luts vs. control")
bhist(hitranscripts$p_value[hitranscripts$status == "OK"], fill='orange', alpha=0.6, 
      add=TRUE)

transcripts = read.table("./isoform_exp.diff", header=TRUE)
rdata = read.table("./isoforms.read_group_tracking", header=TRUE)
fpkmtable = dcast(rdata, formula=tracking_id~replicate+condition, value.var='FPKM')
hiexpr = which(rowMeans(fpkmtable[,-1]) > 1)
length(hiexpr) 

hiID = fpkmtable$tracking_id[hiexpr]
hitranscripts = subset(transcripts, test_id %in% hiID)
hist(hitranscripts$p_value[hitranscripts$status == "OK"], 
     xlab='p-value (highly-expressed, status "OK")', main='Cuffdiff 2.2.1 p-value histogram', col='gray')

genes = read.table("./gene_exp.diff", header=TRUE)
hist(genes$p_value[genes$status=="OK"], xlab='p-value (status "OK")', 
     main="Cuffdiff 2.2.1 p-value histogram (genes)", col='gray')

bgtable = fpkmtable[,-1]
rownames(bgtable) = fpkmtable[,1]
system.time(bgresults <- stattest(gowntable=bgtable, 
                                  pData=data.frame(group=as.numeric(grepl('embryonic', names(bgtable)))),
                                  feature='transcript', covariate='group'))

hist(bgresults$pval[bgresults$id %in% hiID], main='p-value', xlab='p-value (highly-expressed transcripts)', col='gray')

bgres_hi = bgresults[bgresults$id %in% hiID,]
bgres_hi$qval = p.adjust(bgres_hi$pval, 'fdr')
sum(bgres_hi$qval < 0.05)

lib_adj = apply(bgtable, 2, function(x){
  lognz = log2(x[x!=0] + 1)
  q3 = quantile(lognz, 0.75)
  sum(lognz[lognz<q3])
})
y = log2(bgtable+1)
group = as.numeric(grepl('embryonic', names(bgtable)))
x = model.matrix(~group + lib_adj)
fit = lmFit(y, x)
fit = eBayes(fit, trend=TRUE)

limma_p = fit$p.value[,2][names(fit$p.value[,2]) %in% hiID]
hist(limma_p,y = log2(bgtable+1), xlab='p-value (highly-expressed transcripts)', 
     main='Limma p-value histogram', col='gray')


plot(limma_p, bgres_hi$pval, xlab='Limma p-value', ylab='p-value')
hist(abs(limma_p - bgres_hi$pval), main='Absolute differences between Limma and Ballgown p-values', col='gray')
bhist(bgresults$pval[bgresults$id %in% hiID], fill='dodgerblue', 
      alpha=0.6, xlab="p-values", ylab='Frequency', 
      main="Luts vs CTRL")
bhist(hitranscripts$p_value[hitranscripts$status == "OK"], 
      fill='orange', alpha=0.6, add=TRUE)
legend('topright', col=c('dodgerblue', 'orange'), pch=c(15,15), c('Ballgown', 'Cuffdiff'))


