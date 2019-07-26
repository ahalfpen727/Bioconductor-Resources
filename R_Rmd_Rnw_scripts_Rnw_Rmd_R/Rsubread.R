### R code from vignette source 'Rsubread.Rnw'

###################################################
### code chunk number 1: Rsubread.Rnw:67-70
###################################################
library(Rsubread)
ref <- system.file("extdata","reference.fa",package="Rsubread")
buildindex(basename="reference_index",reference=ref)


###################################################
### code chunk number 2: Rsubread.Rnw:87-89
###################################################
reads <- system.file("extdata","reads.txt.gz",package="Rsubread")
align(index="reference_index",readfile1=reads,output_file="alignResults.BAM",phredOffset=64)


###################################################
### code chunk number 3: Rsubread.Rnw:96-100
###################################################
reads1 <- system.file("extdata","reads1.txt.gz",package="Rsubread")
reads2 <- system.file("extdata","reads2.txt.gz",package="Rsubread")
align(index="reference_index",readfile1=reads1,readfile2=reads2,
output_file="alignResultsPE.BAM",phredOffset=64)


###################################################
### code chunk number 4: Rsubread.Rnw:121-131
###################################################
ann <- data.frame(
GeneID=c("gene1","gene1","gene2","gene2"),
Chr="chr_dummy",
Start=c(100,1000,3000,5000),
End=c(500,1800,4000,5500),
Strand=c("+","+","-","-"),
stringsAsFactors=FALSE)
ann
fc_SE <- featureCounts("alignResults.BAM",annot.ext=ann)
fc_SE


###################################################
### code chunk number 5: Rsubread.Rnw:138-140
###################################################
fc_PE <- featureCounts("alignResultsPE.BAM",annot.ext=ann,isPairedEnd=TRUE)
fc_PE


###################################################
### code chunk number 6: Rsubread.Rnw:160-162
###################################################
x <- qualityScores(filename=reads,offset=64,nreads=1000)
x[1:10,1:10]


###################################################
### code chunk number 7: Rsubread.Rnw:187-188
###################################################
propmapped("alignResults.BAM")


