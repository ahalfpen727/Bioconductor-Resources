### R code from vignette source 'ToPASeq.Rnw'
source("http://bioconductor.org/biocLite.R")
biocLite("ToPASeq")
library(ToPASeq);library(DEGraph)
data(Loi2008_DEGraphVignette)

pathways<-pathways("hsapiens", "kegg")[1:5]
pathways[[1]]
str(pathways[[1]])
pathways<-pathways("hsapiens", "kegg")[1:5]
ls()
 top<-TopologyGSA(exprLoi2008, classLoi2008, pathways, type="MA", perms=200)
res(top) # $results
## #99 node labels mapped to the expression data
## #$errors
## #named list()
deg<-DEGraph(exprLoi2008, classLoi2008, pathways, type="MA")
res(deg)


###################################################
### code chunk number 6: ToPASeq.Rnw:197-226 (eval = FALSE)
###################################################
cli<-clipper( exprLoi2008, classLoi2008, pathways,type="MA", method="mean")
## #99 node labels mapped to the expression data
## #Warning messages:
## #1: In getJunctionTreePaths(graph, root) :
## #  The DAG presents cliques that are not connected.
## #2: In prunePaths(clipped, pruneLevel) : pathSummary is NULL
## #3: In getJunctionTreePaths(graph, root) :
## #  The DAG presents cliques that are not connected.
## #4: In prunePaths(clipped, pruneLevel) : pathSummary is NULL

 res(cli)$results[[1]]
spi<-SPIA(exprLoi2008, classLoi2008,pathways , type="MA", logFC.th=-1)
res(spi)

tap<-TAPPA(exprLoi2008, classLoi2008, pathways, type="MA")
res(tap)

Prs<-PRS( exprLoi2008, classLoi2008, pathways, type="MA",  logFC.th=-1, nperm=100)
res(Prs)

 pwe<-PWEA(exprLoi2008, classLoi2008, pathways,  type="MA", nperm=100)
 res(pwe)
## #$results

library(gageData)
data(hnrnp.cnts)
hnrnp.cnts<-hnrnp.cnts[rowSums(hnrnp.cnts)>0,]
group<-c(rep("sample",4), rep("control",4))
pathways<-pathways("hsapiens", "kegg")
 top<-TopologyGSA(hnrnp.cnts, group, pathways[1:3], type="RNASeq", nperm=1000)
## #528 node labels mapped to the expression data
## #Average coverage 83.16538
## #0 (out of 10) pathways without a mapped node
## #Normalization method was not specified. TMM used as default
## #Acute myeloid leukemia
## #Adherens junction
## #Adipocytokine signaling pathway
res(top)
#data frame with 0 columns and 1 rows
deg<-DEGraph(hnrnp.cnts, group, pathways, type="RNASeq")
res(deg)[[1]][[1]]
 cli<-clipper(hnrnp.cnts, group, pathways, type="RNASeq", method="mean")

 ## #530 node labels mapped to the expression data
## #Average coverage 82.98681 %
## #0 (out of 10) pathways without a mapped node
## #1  pathways were filtered out
## #Analysing pathway:
res(cli)$results[[1]][1:2,]

spi<-SPIA(hnrnp.cnts, group, pathways, type="RNASeq", logFC.th=-1)
res(spi)

tap<-TAPPA(hnrnp.cnts, group, pathways, type="RNASeq")
res(tap)

Prs<-PRS(hnrnp.cnts, group, pathways, type="RNASeq",  logFC.th=-1, nperm=100)
res(Prs)
## pwe<-PWEA(hnrnp.cnts, group, pathways, type="RNASeq", nperm=100)
## #528 node labels mapped to the expression data
## #Average coverage 83.16538
## #0 (out of 10) pathways without a mapped node
## #Acute myeloid leukemia
## #Adherens junction
## #Adipocytokine signaling pathway
## #Adrenergic signaling in cardiomyocytes
## #African trypanosomiasis
## #Alanine, aspartate and glutamate metabolism
## #Alcoholism
## #Aldosterone-regulated sodium reabsorption
## #Allograft rejection
## #alpha-Linolenic acid metabolism
## res(pwe)
## #                                                   ES    p     p.adj
## #Acute myeloid leukemia                      0.3526104 0.29 0.4142857
## #Adherens junction                           0.3829831 1.00 1.0000000
## #Adipocytokine signaling pathway             0.3102945 1.00 1.0000000
## #Adrenergic signaling in cardiomyocytes      0.3611207 0.20 0.3333333
## #African trypanosomiasis                     0.3272899 0.20 0.3333333
## #Alanine, aspartate and glutamate metabolism 0.2720946 0.20 0.3333333
## #Alcoholism                                  0.4708293 0.86 1.0000000
## #Aldosterone-regulated sodium reabsorption   0.3951037 0.20 0.3333333
## #Allograft rejection                         0.9421248 0.03 0.3000000
## #alpha-Linolenic acid metabolism             0.6587026 0.20 0.3333333


###################################################
### code chunk number 19: plot1 (eval = FALSE)
###################################################
## #Fails during check
## library(gageData)
## data(hnrnp.cnts)
## 
## group<-c(rep("sample",4), rep("control",4))
## hnrnp.cnts<-hnrnp.cnts[rowSums(hnrnp.cnts)>0,]
## res<-clipper(hnrnp.cnts, group, pathways[1:2], type="RNASeq", testCliques=TRUE)
## 
## plot(res,1, pathways)
## 


###################################################
### code chunk number 20: plot2 (eval = FALSE)
###################################################
## library(gageData)
## data(hnrnp.cnts)
## 
## group<-c(rep("sample",4), rep("control",4))
## hnrnp.cnts<-hnrnp.cnts[rowSums(hnrnp.cnts)>0,]
## pathways<-pathways("hsapiens", "kegg")[50:55]
## spi<-SPIA(hnrnp.cnts, group, pathways, type="RNASeq", logFC.th=-1)
## plot(spi,"Complement and coagulation cascades", pathways, fontsize=50)
## 


