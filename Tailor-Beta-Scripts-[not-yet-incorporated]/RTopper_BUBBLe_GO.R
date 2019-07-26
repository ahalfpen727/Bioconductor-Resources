### R code from vignette source 'RTopper.Rnw'

###################################################
### code chunk number 1: start
###################################################
options(width=50)
rm(list=ls())

## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("RTopper")

###################################################
### code chunk number 2: loadData
###################################################
library(limma)
library(GO.db)
browseVignettes("RTopper")
library(KEGG.db)
library(org.Hs.eg.db)
library(RTopper)
data(exampleData)
ls()
class(dat)
names(dat)
sapply(dat,class)
sapply(dat,dim)
dim(pheno)
str(pheno)

###data structure
lapply(dat,function(x) head(x)[,1:3])
sum(rownames(dat[[1]])%in%rownames(dat[[2]]))
sum(rownames(dat[[2]])%in%rownames(dat[[3]]))
###################################################
# list the contents
org.Hs.eg()
KEGG()
GO()

###################################################
kegg <- as.list(org.Hs.egPATH2EG)
go <- as.list(org.Hs.egGO2ALLEGS)
length(kegg)
length(go)
str(kegg[1:5])
names(kegg)[1:5]
str(go[1:5])
names(go)[1:5]


###################################################
### code chunk number 6: convertIDs
###################################################
kegg <- lapply(kegg[sample(1:length(kegg),5)],function(x) unique(unlist(mget(x,org.Hs.egSYMBOL))))
go <- lapply(go[sample(1:length(go),5)],function(x) unique(unlist(mget(x,org.Hs.egSYMBOL))))
str(kegg)
str(go)

###################################################
### code chunk number 7: annotateFGS
###################################################

names(kegg) <- paste(names(kegg),unlist(mget(names(kegg),KEGGPATHID2NAME)),sep=".")
names(kegg)

names(go) <- paste(names(go),Term(names(go)),sep=".")
names(go)

###################################################
### code chunk number 9: listFGS
###################################################
fgsList <- list(go=go,kegg=kegg)
str(fgsList)

###################################################
### code chunk number 10: convertToDr
###################################################
dataDr <- convertToDr(dat, pheno, 4)
class(dataDr)
length(dataDr)
names(dataDr)[1:5]
str(dataDr[1:2])


###################################################
### code chunk number 11: integratedScore
###################################################
bicStatInt <- computeDrStat(dataDr, columns = c(1:4), method="bic", integrate = TRUE)
names(bicStatInt)
str(bicStatInt)


###################################################
### code chunk number 12: separateScore
###################################################
bicStatSep <- computeDrStat(dataDr, columns = c(1:4), method="bic", integrate = FALSE)
names(bicStatSep)
str(bicStatSep)


###################################################
### code chunk number 13: runGSEbatchArgs
###################################################
args(runBatchGSE)


###################################################
### code chunk number 14: runBatchGSE.int1
###################################################
gseABS.int <- runBatchGSE(dataList=bicStatInt, fgsList=fgsList)
gseABS.int <- runBatchGSE(dataList=bicStatInt, fgsList=fgsList,
				 absolute=TRUE, type="f", alternative="mixed")


###################################################
### code chunk number 15: runBatchGSE.int2
###################################################
gseUP.int <- runBatchGSE(dataList=bicStatInt, fgsList=fgsList,
				 absolute=FALSE, type="t", alternative="up")
gseDW.int <- runBatchGSE(dataList=bicStatInt, fgsList=fgsList,
				 absolute=FALSE, type="t", alternative="down")
gseBOTH.int <- runBatchGSE(dataList=bicStatInt, fgsList=fgsList,
				 absolute=FALSE, type="t", alternative="either")


###################################################
### code chunk number 16: runBatchGSE.int3
###################################################
gseABSsim.int <- runBatchGSE(dataList=bicStatInt, fgsList=fgsList,
				    absolute=TRUE, type="f", alternative="mixed",
				    ranks.only=FALSE, nsim=1000)
gseUPsim.int <- runBatchGSE(dataList=bicStatInt, fgsList=fgsList,
				    absolute=FALSE, type="t", alternative="up",
				    ranks.only=FALSE, nsim=1000)

str(gseUP.int)
gseABSsim.int

gseUP.int.2 <- runBatchGSE(dataList=bicStatInt, fgsList=fgsList,
				 absolute=FALSE, gseFunc=wilcoxGST, alternative="up")

str(gseUP.int.2)
all(gseUP.int.2$go==gseUP.int$go)


###################################################
### code chunk number 20: runBatchGSE.altFunc2
###################################################
gseFunc <- function (selected, statistics, threshold) {
	diffExpGenes <- statistics > threshold
	tab <- table(diffExpGenes, selected)
	pVal <- fisher.test(tab)[["p.value"]]
	}
gseUP.int.3 <- runBatchGSE(dataList=bicStatInt, fgsList=fgsList,
				 absolute=FALSE, gseFunc=gseFunc, threshold=7.5)


###################################################
### code chunk number 21: runBatchGSE.format3
###################################################
str(gseUP.int.3)
data.frame(fisher=gseUP.int.3$integrated$kegg,wilcoxon=gseUP.int$integrated$kegg)

gseABS.sep <- runBatchGSE(dataList=bicStatSep, fgsList=fgsList)
gseABS.geoMean.sep <- combineGSE(gseABS.sep, method="geometricMean")
gseABS.max.sep <- combineGSE(gseABS.sep, method="max")

names(gseABS.sep)
str(gseABS.sep)
str(gseABS.geoMean.sep)
gseABS.geoMean.sep


###################################################
### code chunk number 25: adjustP
###################################################
gseABS.int.BH <- adjustPvalGSE(gseABS.int)
gseABS.int.holm <- adjustPvalGSE(gseABS.int, proc = "Holm")

names(gseABS.int.BH)
names(gseABS.int.holm)
str(gseABS.int.BH)
str(gseABS.int.holm)

#install.packages("googleVis")  ## If you need to install the package
library(googleVis)
M<-gvisMotionChart(Fruits, "Fruit", "Year")
plot(M)
cat(M$html$chart, file = "tmp.html")

 

library(clusterProfiler)
data(geneList)
str(geneList)	
gene = names(geneList)[geneList > 1]
str(gene)	
yy = enrichDO(gene, pvalueCutoff=0.05)
	summary(yy)
data(gcSample)
gcSample
class(gcSample)
str(gcSample)
head(gcSample$X1)
x=compareCluster(yy)
x=compareCluster(gcSample, fun='enrichDO')
str(x)

# still some kinks not sure 
#enricher()
#enrichGO()
#str(ClusterList)
#CompareGO_BP=compareCluster(ClusterList, fun="enrichGO", pvalueCutoff=0.01, pAdjustMethod="BH", OrgDb=org.Hs.eg.db,ont="BP",readable=T)
#dotplot(CompareGO_BP, showCategory=10, title="GO - Biological Process")

 enrichDO(gene, ont = "DO", pvalueCutoff = 0.05, pAdjustMethod = "BH",
  universe, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2,
  readable = FALSE)

 # (gene -- vector of entrex ids),
 #  (pAdjustMethod --"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")

 dotplot(x, showCategory=5, includeAll=FALSE)
dotplot(x, showCategory=5)
x=compareCluster(gcSample, fun='enrichDO')
str(x)
dotplot(x, showCategory=5, includeAll=FALSE)
dotplot(x, showCategory=5)


