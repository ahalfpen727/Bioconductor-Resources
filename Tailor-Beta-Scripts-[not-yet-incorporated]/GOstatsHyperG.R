### R code from vignette source 'GOstatsForUnsupportedOrganisms.Rnw'

###################################################
### code chunk number 1: available Schemas
###################################################
library("AnnotationForge")
available.dbschemas()


###################################################
### code chunk number 2: Acquire annotation data
###################################################
library("org.Hs.eg.db")
frame = toTable(org.Hs.egGO)
goframeData = data.frame(frame$go_id, frame$Evidence, frame$gene_id)
head(goframeData)


###################################################
### code chunk number 3: transformGOFrame
###################################################
goFrame=GOFrame(goframeData,organism="Homo sapiens")
goAllFrame=GOAllFrame(goFrame)


###################################################
### code chunk number 4: Make GSC
###################################################
library("GSEABase")
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())


###################################################
### code chunk number 5: <make parameter
###################################################
library("GOstats")
universe = Lkeys(org.Hs.egGO)
genes = universe[1:500]
params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params", 
                             geneSetCollection=gsc, 
                             geneIds = genes, 
                             universeGeneIds = universe, 
                             ontology = "MF", 
                             pvalueCutoff = 0.05, 
                             conditional = FALSE, 
                             testDirection = "over")


###################################################
### code chunk number 6: call HyperGTest
###################################################
Over <- hyperGTest(params)
head(summary(Over))


###################################################
### code chunk number 7: KEGGFrame object
###################################################
frame = toTable(org.Hs.egPATH)
keggframeData = data.frame(frame$path_id, frame$gene_id)
head(keggframeData)
keggFrame=KEGGFrame(keggframeData,organism="Homo sapiens")


###################################################
### code chunk number 8: KEGG Parameters
###################################################
gsc <- GeneSetCollection(keggFrame, setType = KEGGCollection())
universe = Lkeys(org.Hs.egGO)
genes = universe[1:500]
kparams <- GSEAKEGGHyperGParams(name="My Custom GSEA based annot Params", 
                                geneSetCollection=gsc, 
                                geneIds = genes, 
                                universeGeneIds = universe,  
                                pvalueCutoff = 0.05, 
                                testDirection = "over")
kOver <- hyperGTest(params)
head(summary(kOver))


###################################################
### code chunk number 9: info
###################################################
toLatex(sessionInfo())



### R code from vignette source 'GOstatsHyperG.Rnw'

###################################################
### code chunk number 1: Setup
###################################################
library("ALL")
library("hgu95av2.db")
library("GO.db")
library("annotate")
library("genefilter")
library("GOstats")
library("RColorBrewer")
library("xtable")
library("Rgraphviz")


###################################################
### code chunk number 2: universeSizeVsPvalue
###################################################
hg_tester <- function(size) {
    numFound <- 10
    numDrawn <- 400
    numAtCat <- 40
    numNotAtCat <- size - numAtCat
    phyper(numFound-1, numAtCat, numNotAtCat, numDrawn, lower.tail=FALSE)
}
pv1000 <- hg_tester(1000)
pv5000 <- hg_tester(5000)


###################################################
### code chunk number 3: bcrAblOrNegSubset
###################################################
data(ALL, package="ALL")
## For this data we can have ALL1/AF4 or BCR/ABL
subsetType <- "ALL1/AF4"

## Subset of interest: 37BRC/ABL + 42NEG = 79 samples
Bcell <- grep("^B", as.character(ALL$BT))
bcrAblOrNegIdx <- which(as.character(ALL$mol) %in% c("NEG", subsetType))

bcrAblOrNeg <- ALL[, intersect(Bcell, bcrAblOrNegIdx)]
bcrAblOrNeg$mol.biol = factor(bcrAblOrNeg$mol.biol)


###################################################
### code chunk number 4: nsFiltering-noEntrez
###################################################
## Remove genes that have no entrezGene id
entrezIds <- mget(featureNames(bcrAblOrNeg), envir=hgu95av2ENTREZID)
haveEntrezId <- names(entrezIds)[sapply(entrezIds, function(x) !is.na(x))]
numNoEntrezId <- length(featureNames(bcrAblOrNeg)) - length(haveEntrezId)
bcrAblOrNeg <- bcrAblOrNeg[haveEntrezId, ]


###################################################
### code chunk number 5: nsFiltering-noGO
###################################################
## Remove genes with no GO mapping
haveGo <- sapply(mget(featureNames(bcrAblOrNeg), hgu95av2GO),
                 function(x) {
                     if (length(x) == 1 && is.na(x)) 
                       FALSE 
                     else TRUE
                 })
numNoGO <- sum(!haveGo)
bcrAblOrNeg <- bcrAblOrNeg[haveGo, ]


###################################################
### code chunk number 6: nsFiltering-IQR
###################################################
## Non-specific filtering based on IQR
iqrCutoff <- 0.5
bcrAblOrNegIqr <- apply(exprs(bcrAblOrNeg), 1, IQR)
selected <- bcrAblOrNegIqr > iqrCutoff

## Drop those that are on the Y chromosome
## because there is an imbalance of men and women by group
chrN <- mget(featureNames(bcrAblOrNeg), envir=hgu95av2CHR)
onY <- sapply(chrN, function(x) any(x=="Y"))
onY[is.na(onY)] <- FALSE
selected <- selected & !onY

nsFiltered <- bcrAblOrNeg[selected, ]


###################################################
### code chunk number 7: nsFiltering-unique
###################################################
numNsWithDups <- length(featureNames(nsFiltered))
nsFilteredIqr <- bcrAblOrNegIqr[selected]
uniqGenes <- findLargest(featureNames(nsFiltered), nsFilteredIqr, 
                         "hgu95av2")
nsFiltered <- nsFiltered[uniqGenes, ]
numSelected <- length(featureNames(nsFiltered))

##set up some colors
BCRcols = ifelse(nsFiltered$mol == subsetType, "goldenrod", "skyblue")
cols = brewer.pal(10, "RdBu")


###################################################
### code chunk number 8: defineGeneUniverse
###################################################
## Define gene universe based on results of non-specific filtering
affyUniverse <- featureNames(nsFiltered)
entrezUniverse <- unlist(mget(affyUniverse, hgu95av2ENTREZID))
if (any(duplicated(entrezUniverse)))
  stop("error in gene universe: can't have duplicate Entrez Gene Ids")

## Also define an alternate universe based on the entire chip
chipAffyUniverse <- featureNames(bcrAblOrNeg)
chipEntrezUniverse <- mget(chipAffyUniverse, hgu95av2ENTREZID)
chipEntrezUniverse <- unique(unlist(chipEntrezUniverse)) 


###################################################
### code chunk number 9: parametric1
###################################################
ttestCutoff <- 0.05
ttests = rowttests(nsFiltered, "mol.biol")

smPV = ttests$p.value < ttestCutoff

pvalFiltered <- nsFiltered[smPV, ]
selectedEntrezIds <- unlist(mget(featureNames(pvalFiltered),
                                 hgu95av2ENTREZID))


###################################################
### code chunk number 10: withYourData1 (eval = FALSE)
###################################################
## entrezUniverse <- unlist(mget(featureNames(yourData), 
##                               hgu95av2ENTREZID))
## if (any(duplicated(entrezUniverse)))
##   stop("error in gene universe: can't have duplicate Entrez Gene Ids")
## 
## pvalFiltered <- yourData[hasSmallPvalue, ]
## selectedEntrezIds <- unlist(mget(featureNames(pvalFiltered),
##                                  hgu95av2ENTREZID))


###################################################
### code chunk number 11: standardHyperGeo
###################################################
hgCutoff <- 0.001
params <- new("GOHyperGParams",
              geneIds=selectedEntrezIds,
              universeGeneIds=entrezUniverse,
              annotation="hgu95av2.db",
              ontology="BP",
              pvalueCutoff=hgCutoff,
              conditional=FALSE,
              testDirection="over")



###################################################
### code chunk number 12: condHyperGeo
###################################################
paramsCond <- params
conditional(paramsCond) <- TRUE


###################################################
### code chunk number 13: standardHGTEST
###################################################
hgOver <- hyperGTest(params)
hgCondOver <- hyperGTest(paramsCond)



###################################################
### code chunk number 14: HGTestAns
###################################################
hgOver
hgCondOver


###################################################
### code chunk number 15: summaryEx
###################################################
df <- summary(hgOver)
names(df)                               # the columns
dim(summary(hgOver, pvalue=0.1))
dim(summary(hgOver, categorySize=10))


###################################################
### code chunk number 16: resultAccessors
###################################################
pvalues(hgOver)[1:3]

oddsRatios(hgOver)[1:3]

expectedCounts(hgOver)[1:3]

geneCounts(hgOver)[1:3]
universeCounts(hgOver)[1:3]

length(geneIds(hgOver))
length(geneIdUniverse(hgOver)[[3]])

## GOHyperGResult _only_
## (NB: edges go from parent to child)
goDag(hgOver)

geneMappedCount(hgOver)
universeMappedCount(hgOver)

conditional(hgOver)
testDirection(hgOver)
testName(hgOver)



###################################################
### code chunk number 17: htmlReportExample
###################################################
htmlReport(hgCondOver, file="ALL_hgco.html")


###################################################
### code chunk number 18: helperFunc
###################################################
sigCategories <- function(res, p) {
    if (missing(p))
      p <- pvalueCutoff(res)
    pv <- pvalues(res)
    goIds <- names(pv[pv < p])
    goIds
}


###################################################
### code chunk number 19: plotFuns
###################################################
coloredGoPlot <- function(ccMaxGraph, hgOver, hgCondOver) {
    nodeColors <- sapply(nodes(ccMaxGraph),
                         function(n) {
                             if (n %in% sigCategories(hgCondOver))
                               "dark red"
                             else if (n %in% sigCategories(hgOver))
                               "pink"
                             else 
                               "gray"
                         })
    nattr <- makeNodeAttrs(ccMaxGraph,
                           label=nodes(ccMaxGraph),
                           shape="ellipse",
                           fillcolor=nodeColors,
                           fixedsize=FALSE)
    plot(ccMaxGraph, nodeAttrs=nattr)
}

getMaxConnCompGraph <- function(hgOver, hgCondOver) {
    ##uGoDagRev <- ugraph(goDag(hgOver))
    sigNodes <- sigCategories(hgOver)
    ##adjNodes <- unlist(adj(uGoDagRev, sigNodes))
    adjNodes <- unlist(adj(goDag(hgOver), sigNodes))
    displayNodes <- unique(c(sigNodes, adjNodes))
    displayGraph <- subGraph(displayNodes, goDag(hgOver))
    cc <- connComp(displayGraph)
    ccSizes <- listLen(cc)
    ccMaxIdx <- which(ccSizes == max(ccSizes))
    ccMaxGraph <- subGraph(cc[[ccMaxIdx]], displayGraph)
    ccMaxGraph 
}



###################################################
### code chunk number 20: info
###################################################
toLatex(sessionInfo())


