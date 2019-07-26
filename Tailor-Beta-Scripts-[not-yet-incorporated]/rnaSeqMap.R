### R code from vignette source 'rnaSeqMap.Rnw'
### Encoding: UTF-8
source("http://bioconductor.org/workflows.R")
workflowInstall("rnaSeqMap")
library(easyRNASeq)
###################################################
### code chunk number 1: camelRegions (eval = FALSE)
###################################################
         rs <- newSeqReads(ch,st, en, str);
         rs <- getBamData(rs,idx.both, cvd=cvd)
         nd <- getCoverageFromRS(rs, idx.both) 


###################################################
### code chunk number 2: camel1 (eval = FALSE)
###################################################
  idxT <- which(samples$condition=="T")
  idxC <- which(samples$condition=="C")


###################################################
### code chunk number 3: camel2 (eval = FALSE)
###################################################
  regions.gR <- rnaSeqMap:::.fiveCol2GRanges(tmp)


###################################################
### code chunk number 4: camel3 (eval = FALSE)
###################################################
  regionsCamelMeasures <- gRanges2CamelMeasures(regions.gR,samples,idxT,idxC,sums=sums,progress=10)


###################################################
### code chunk number 5: camel4 (eval = FALSE)
###################################################
 idx <- which(regionsCamelMeasures[,"covDensC1"]>10 | regionsCamelMeasures[,"covDensC1"]>10)
regionsCamelMeasures <- regionsCamelMeasures[idx, ]


###################################################
### code chunk number 6: camel5 (eval = FALSE)
###################################################
 o <- order(regionsCamelMeasures[,"QQ.mm"], decreasing=T)
 regionsCamelMeasures <- regionsCamelMeasures [o, ] 


###################################################
### code chunk number 7: Lindell (eval = FALSE)
###################################################
 nd.AL <- findRegionsAsND(nd, 15, minsup=5)


