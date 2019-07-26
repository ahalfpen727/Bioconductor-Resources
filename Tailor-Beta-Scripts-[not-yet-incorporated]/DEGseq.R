### R code from vignette source 'DEGseq.Rnw'

###################################################
### code chunk number 1: DEGseq.Rnw:189-195
###################################################
  library(DEGseq)
  geneExpFile <- system.file("extdata", "GeneExpExample5000.txt", package="DEGseq")
  geneExpMatrix1 <- readGeneExp(file=geneExpFile, geneCol=1, valCol=c(7,9,12,15,18))
  geneExpMatrix2 <- readGeneExp(file=geneExpFile, geneCol=1, valCol=c(8,10,11,13,16))
  write.table(geneExpMatrix1[30:31,],row.names=FALSE)
  write.table(geneExpMatrix2[30:31,],row.names=FALSE)


###################################################
### code chunk number 2: DEGseq.Rnw:203-208
###################################################
layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow=TRUE))
par(mar=c(2, 2, 2, 2))
DEGexp(geneExpMatrix1=geneExpMatrix1, geneCol1=1, expCol1=c(2,3,4,5,6), groupLabel1="kidney",
       geneExpMatrix2=geneExpMatrix2, geneCol2=1, expCol2=c(2,3,4,5,6), groupLabel2="liver",
       method="MARS")

genes<-row.names(g.cnt.matrix)
g.u.cnt.ma<-cbind(genes, g.cnt.matrix[,under.group])
str(g.u.cnt.ma)
g.o.cnt.ma<-cbind(genes, g.cnt.matrix[,over.group])
str(g.o.cnt.ma)
DEGexp(geneExpMatrix1=g.u.cnt.ma,geneCol1=1, expCol1=c(2:9), geneExpMatrix2=g.o.cnt.ma, geneCol2=1, expCol2=c(2:9),pValue=0.05, outputDir="/home/drew/umb_triley/urine1/cuffdiff_results_hg38_default/LUTS-over-CTRL/")

###################################################
### code chunk number 3: DEGseq.Rnw:222-229
###################################################
layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow=TRUE))
par(mar=c(2, 2, 2, 2))
DEGexp(geneExpMatrix1=geneExpMatrix1, expCol1=2, groupLabel1="kidneyR1L1",
       geneExpMatrix2=geneExpMatrix2, expCol2=2, groupLabel2="liverR1L2",
       replicateExpMatrix1=geneExpMatrix1, expColR1=3, replicateLabel1="kidneyR1L3",
       replicateExpMatrix2=geneExpMatrix1, expColR2=4, replicateLabel2="kidneyR1L7",
       method="MATR")


###################################################
### code chunk number 4: DEGseq.Rnw:244-252
###################################################
  kidneyR1L1 <- system.file("extdata", "kidneyChr21.bed.txt", package="DEGseq")
  liverR1L2  <- system.file("extdata", "liverChr21.bed.txt", package="DEGseq")
  refFlat    <- system.file("extdata", "refFlatChr21.txt", package="DEGseq")
  mapResultBatch1 <- c(kidneyR1L1)  ## only use the data from kidneyR1L1 and liverR1L2
  mapResultBatch2 <- c(liverR1L2)
  outputDir <- file.path(tempdir(), "DEGseqExample")
  DEGseq(mapResultBatch1, mapResultBatch2, fileFormat="bed", refFlat=refFlat,
         outputDir=outputDir, method="LRT")


###################################################
### code chunk number 5: DEGseq.Rnw:266-272
###################################################
  kidneyR1L1 <- system.file("extdata", "kidneyChr21.bed.txt", package="DEGseq")
  refFlat    <- system.file("extdata", "refFlatChr21.txt", package="DEGseq")
  mapResultBatch <- c(kidneyR1L1)
  output <- file.path(tempdir(), "kidneyChr21.bed.exp")
  exp <- getGeneExp(mapResultBatch, refFlat=refFlat, output=output)
  write.table(exp[30:32,], row.names=FALSE)


###################################################
### code chunk number 6: DEGseq.Rnw:280-283
###################################################
  geneExpFile <- system.file("extdata", "GeneExpExample1000.txt", package="DEGseq")
  exp <- readGeneExp(file=geneExpFile, geneCol=1, valCol=c(7,9,12,15,18,8,10,11,13,16))
  write.table(exp[30:32,], row.names=FALSE)


