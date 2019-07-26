### R code from vignette source 'seq2pathwaypackage.Rnw'

###################################################
### code chunk number 1: setup
###################################################
library("seq2pathway.data") 
library("seq2pathway") 


###################################################
### code chunk number 2: shoefunction(runseq2pathway)
###################################################
head(runseq2pathway, n=8)


###################################################
### code chunk number 3: Chipseq_Peak_demo
###################################################
data(Chipseq_Peak_demo)
class(Chipseq_Peak_demo)
head(Chipseq_Peak_demo)


###################################################
### code chunk number 4: GRanges_demo
###################################################
data(GRanges_demo)
class(GRanges_demo)
GRanges_demo[1:3,]


###################################################
### code chunk number 5: MsigDB_C5
###################################################
data(MsigDB_C5,package="seq2pathway.data")
class(MsigDB_C5)     


###################################################
### code chunk number 6: Chipseq_Peak_demo
###################################################
data(Chipseq_Peak_demo)
head(Chipseq_Peak_demo)    


###################################################
### code chunk number 7: runseq2gene
###################################################
Chipseq_anno <- runseq2gene(
                  inputfile=Chipseq_Peak_demo, 
                  genome="hg38", adjacent=FALSE, SNP=FALSE, search_radius=1000, 
                  PromoterStop=FALSE,NearestTwoDirection=TRUE)
class(Chipseq_anno)
head(Chipseq_anno[[1]])


###################################################
### code chunk number 8: MsigDB_C5
###################################################
## give the previously defined gene-sets  
data(MsigDB_C5,package="seq2pathway.data")
class(MsigDB_C5)    
## load the gene-level measurements, here is an example of ChIP-seq scores
data(dat_chip)
head(dat_chip)


###################################################
### code chunk number 9: dat_gene2path_chip
###################################################
data(dat_gene2path_chip,package="seq2pathway.data")
names(dat_gene2path_chip)
class(dat_gene2path_chip$gene2pathway_result.2)
names(dat_gene2path_chip$gene2pathway_result.2)
dat_gene2path_chip$gene2pathway_result.2$GO_BP[1:3,]

class(dat_gene2path_chip$gene2pathway_result.FET)
names(dat_gene2path_chip$gene2pathway_result.FET)
colnames(dat_gene2path_chip$gene2pathway_result.FET$GO_BP)
dat_gene2path_chip$gene2pathway_result.FET$GO_BP[1:3,-2]


###################################################
### code chunk number 10: MsigDB_C5
###################################################
data(MsigDB_C5,package="seq2pathway.data")
class(MsigDB_C5)  


###################################################
### code chunk number 11: dat_chip
###################################################
data(dat_chip)
head(dat_chip)


###################################################
### code chunk number 12: FisherTest_GO_BP_MF_CC
###################################################
data(dat_chip)
head(dat_chip)


###################################################
### code chunk number 13: importseq2pathway
###################################################
require(seq2pathway)


###################################################
### code chunk number 14: sessionInfo
###################################################
sessionInfo();


