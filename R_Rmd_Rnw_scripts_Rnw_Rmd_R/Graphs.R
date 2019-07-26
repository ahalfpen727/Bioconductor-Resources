###################################################
### chunk number 1: setup1
###################################################
library("RbcBook1")
library("Rgraphviz")
library("Biobase")
stopifnot(package.version("RbcBook1") >= package_version("0.2.1"))
stopifnot(package.version("Rgraphviz") >= package_version("1.5.6"))
stopifnot(package.version("CoCiteStats") >= package_version("0.5.3"))
humanLLMappings

###################################################
### chunk number 2: dataIMCA
###################################################
library("graph")
data("integrinMediatedCellAdhesion")
class(IMCAGraph)


###################################################
### chunk number 3: SOS
###################################################
s <- acc(IMCAGraph, "SOS")


###################################################
### chunk number 4: degree
###################################################
deg <- degree(IMCAGraph)$outDegree
deg[which.max(deg)]


###################################################
### chunk number 5: tfG1
###################################################
library("GO")
library("GOstats")


###################################################
### chunk number 6: tfG2
###################################################
GOTERM$"GO:0003700"


###################################################
### chunk number 7: tfG3
###################################################
tfG   <- GOGraph("GO:0003700", GOMFPARENTS)


###################################################
### chunk number 8: tfG4
###################################################
library("Rgraphviz")
nL    <- nodes(tfG)
ggt   <- unlist(getGOTerm(nL))

labs <- character(length(nL))
names(labs) <- names(nL)
labs[sub("MF.", "", names(ggt))] <- ggt
labs["top"] <- "GO"

nattr <- makeNodeAttrs(tfG,
                       label     = labs[nodes(tfG)],
                       shape     = "ellipse",
                       fillcolor = "#f2f2f2",
                       fixedsize = FALSE)


###################################################
### chunk number 9: tfGplot
###################################################
plot(tfG, nodeAttrs=nattr)


###################################################
### chunk number 10: tf.children
###################################################
tfch <- GOMFCHILDREN$"GO:0003700"
tfchild <- mget(tfch, GOTERM)


###################################################
### chunk number 11: literatureGraphCalc
###################################################
library("humanLLMappings")

g1   <- "1736"        ##DKC1: 1736
paps <- humanLLMappingsLL2PMID[[g1]]

## throw out the Strausberg one
longpaps <- c("15302935", "12477932")
paps <- paps[!(paps %in% longpaps)]

genes <- lapply(paps, function(p) humanLLMappingsPMID2LL[[p]])
names(genes) <- paps

LLstring <- function(i) paste("LL", i, sep=":")
PMstring <- function(i) paste("PM", i, sep=":")

nd <- c(LLstring(unique(unlist(genes))),
        PMstring(paps))
ed        <- lapply(nd, function(z) list(edges=integer(0)))
names(ed) <- nd

for(i in 1:length(genes)) {
  p <- PMstring(names(genes)[i])
  g <- LLstring(genes[[i]])
  ed[[p]] <- list(edges=match(g, nd))
}
g <- new("graphNEL", nodes=nd, edgeL=ed, edgemode="directed")

nt <- match(substr(nodes(g), 1, 2), c("LL", "PM"))
nattr <-  makeNodeAttrs(g,
                        fillcolor=c("#8da0cb", "#e78ac3")[nt],
                        shape=c("ellipse", "rect")[nt],
                        fixedsize=FALSE)


###################################################
### chunk number 12: literatureGraphPlot
###################################################
plot(g, "neato", nodeAttrs=nattr)


###################################################
### chunk number 1: loadlibs
###################################################
library("RbcBook1")
library("Rgraphviz")
library("RBGL")
library("graph")
library("RColorBrewer")
library("geneplotter")
library("GOstats")
library("GO")

histcolor <- brewer.pal(7, "Pastel1")[5]
twocolors <- c(brewer.pal(11,"RdYlGn")[7], brewer.pal(11,"RdYlBu")[7])
sixcolors <- brewer.pal(6, "Dark2")
hundredcolors <- colorRampPalette(
      c(brewer.pal(9, "Pastel1"), brewer.pal(8, "Pastel2")))(100)

graphVizAttrs <- getDefaultAttrs()
graphVizAttrs$node$fillcolor <- twocolors[1]


###################################################
### chunk number 2: sgraph
###################################################
myNodes <- c('s', 'p', 'q', 'r')
myEdges <- list(s=list(edges=c('p', 'q')),
                p=list(edges=c('p', 'q')),
                q=list(edges=c('p', 'r')),
                r=list(edges=c('s')))
g <- new('graphNEL', nodes=myNodes, edgeL=myEdges, edgemode='directed')
g


###################################################
### chunk number 3: sgraph
###################################################
plot(g, attrs=graphVizAttrs)


###################################################
### chunk number 4: prepareSetOps
###################################################
V <- letters[1:4]
set.seed(4713)
ug1 <- randomGraph(V, 1, .55)
ug2 <- randomGraph(V, 1, .55)

myPlot <- function(g, i, tit) {
  plot(g, nodeAttrs=makeNodeAttrs(g, fillcolor=sixcolors[i]), main=tit)
}


###################################################
### chunk number 5: setOps1
###################################################
par(cex.main=2)
myPlot(ug1, 1, "ug1")


###################################################
### chunk number 6: setOps2
###################################################
par(cex.main=2)
myPlot(complement(ug1), 2, "complement(ug1)")


###################################################
### chunk number 7: setOps3
###################################################
par(cex.main=2)
myPlot(ug2, 3, "ug2")


###################################################
### chunk number 8: setOps4
###################################################
par(cex.main=2)
myPlot(complement(ug2), 4, "complement(ug2)")


###################################################
### chunk number 9: setOps5
###################################################
par(cex.main=2)
myPlot(intersection(ug1,ug2), 5, "intersection(ug1, ug2)")


###################################################
### chunk number 10: setOps6
###################################################
par(cex.main=2)
myPlot(union(ug1,ug2), 6, "union(ug1, ug2)")

###################################################
### chunk: load0
###################################################
library("BiocCaseStudies")
#Installing regular R packages:
install.packages(“vegan”)
#Installing Bioconductor packages:
source("http://bioconductor.org/biocLite.R")
biocLite(c("BiocCaseStudies","RCytoscape", "RpsiXML","yeastExpData",
           "ppiStats","apComplex","EnrichmentBrowser","featureCounts",
           "rsbml","paxtoolsr","MLInterfaces","bioDist","EGSEA123"))
###################################################
### chunk: loadLibs
###################################################
library("Biobase");library("graph");library("Rgraphviz");library("RColorBrewer")
library("RBGL")
library(EGSEA123)
library("yeastExpData")
library(RpsiXML)
library("ppiStats")
library("Biobase")
library("RCytoscape")
library("bioDist")
library("genefilter")
library("class")
library("MLInterfaces")
library("hgu95av2.db")
library("annotate")
library("randomForest")
library("apComplex")
library("GSEABase")
library(EnrichmentBrowser); library(featureCounts); library(rsbml);library(paxtoolsr)
###################################################
### chunk: theData
###################################################
data("ccyclered")
data("litG")


###################################################
### chunk: litGNodes
###################################################
nodes(litG)[1:5]


###################################################
### chunk: exA
###################################################
class(ccyclered)
str(ccyclered)
dim(ccyclered)
names(ccyclered)


###################################################
### chunk: connComp
###################################################
cc = connectedComp(litG)
length(cc)
cclens = sapply(cc, length)
table(cclens)


###################################################
### chunk: cc7
###################################################
cc[[7]]


###################################################
### chunk: createSubG
###################################################
ord = order(cclens, decreasing=TRUE)
sg1 = subGraph(cc[[ord[1]]], litG)
sg2 = subGraph(cc[[ord[2]]], litG)


###################################################
### chunk: layoutSubG
###################################################
lsg1 = layoutGraph(sg1, layoutType="neato")
lsg2 = layoutGraph(sg2, layoutType="neato")


###################################################
### chunk: sg1
###################################################
renderGraph(lsg1)


###################################################
### chunk: sg2
###################################################
renderGraph(lsg2)


###################################################
### chunk: extraLayout
###################################################
lay12neato = layoutGraph(sg1, layoutType="dot")
renderGraph(lay12neato,
            graph.pars=list(nodes=list(fillcolor="steelblue2")))
lay12twopi = layoutGraph(sg2, layoutType="twopi")
renderGraph(lay12twopi,
            graph.pars=list(nodes=list(fillcolor="steelblue2")))


###################################################
### chunk: sps
###################################################
sps = sp.between(sg1, "YHR129C", "YOL039W")
pth = sps[[1]]$path_detail
pth


###################################################
### chunk: sglayout
###################################################
fill = rep("steelblue2", length(pth))
names(fill) = pth
nodeRenderInfo(lsg1) = list(fill=fill)
edges = paste(pth[-length(pth)], pth[-1], sep="~")
lwd = rep(5, length(edges))
col = rep("steelblue2", length(edges))
names(lwd) = names(col) = edges
edgeRenderInfo(lsg1) = list(col=col, lwd=lwd)


###################################################
### chunk: sg
###################################################
renderGraph(lsg1)


###################################################
### chunk: diam88sg
###################################################
allp = johnson.all.pairs.sp(sg1)


###################################################
### chunk: diamC
###################################################
max(allp)


###################################################
### chunk: howManyDiam
###################################################
sum(allp == max(allp))


###################################################
### chunk: geneToClusterList
###################################################
clusts = with(ccyclered, split(Y.name, Cluster))


###################################################
### chunk: makeClusterGraph
###################################################
cg = new("clusterGraph", clusters = clusts)


###################################################
### chunk: cgConnectedComp
###################################################
ccClust = connectedComp(cg)
length(ccClust)


###################################################
### chunk: intersect
###################################################
commonG = intersection(cg, litG)


###################################################
### chunk: numEdges
###################################################
numEdges(commonG)


###################################################
### chunk: defNodePerm
###################################################
nodePerm = function (g1, g2, B=1000) {
  n1 = nodes(g1)
  sapply(1:B, function(i) {
    nodes(g1) = sample(n1)
    numEdges(intersection(g1, g2))
  })
}


###################################################
### chunk: nodePermDo
###################################################
data("nPdist")
summary(nPdist)


###################################################
### chunk: nPdist
###################################################
barplot(table(nPdist))


###################################################
### chunk: loadlibs
###################################################
library("RpsiXML")
library("ppiStats")
library("apComplex")
library("xtable")


###################################################
### chunk: psi25int
###################################################
fold <- system.file("/extdata/psi25files", package="RpsiXML")
fn <- file.path(fold, "intact_2008_test.xml")
eg <- parsePsimi25Interaction(psimi25file=fn,
                              psimi25source=INTACT.PSIMI25,
                              verbose=FALSE)


###################################################
### chunk: slotEntries
###################################################
slotNames(eg)


###################################################
### chunk: simpleSlots
###################################################
organismName(eg)
taxId(eg)
releaseDate(eg)


###################################################
### chunk: interactionSlot
###################################################
length(interactions(eg))
class(interactions(eg)[[1]])
interactions(eg)[[1]]
slotNames(interactions(eg)[[1]])


###################################################
### chunk: getBaitPrey
###################################################
egbait = sapply(interactions(eg), bait)
egprey = sapply(interactions(eg), prey)


###################################################
### chunk: baitVec
###################################################
egbait
egprey


###################################################
### chunk: interactors
###################################################
interactors(eg)[1:2]
bt <- egbait[4]
annBt <- interactors(eg)[[bt]]
annBt
xref(annBt)


###################################################
### chunk: psi25complex
###################################################
fn2 = file.path(fold, "intact_complexSample.xml")
comps = parsePsimi25Complex(fn2, INTACT.PSIMI25)
slotNames(comps)


###################################################
### chunk: complex1
###################################################
length(complexes(comps))
class(complexes(comps)[[1]])
slotNames(complexes(comps)[[1]])


###################################################
### chunk: showComplex1
###################################################
complexes(comps)[[1]]


###################################################
### chunk: intactgraph
###################################################
s1 = file.path(fold, "human_stelzl-2005-1_01.xml")
s2 = file.path(fold, "human_stelzl-2005-1_02.xml")
stelzl = separateXMLDataByExpt(xmlFiles=c(s1,s2),
                               psimi25source=INTACT.PSIMI25,
                               type="direct", directed=TRUE,
                               abstract=TRUE, verbose=FALSE)
class(stelzl[[1]])
stelzl1 = removeSelfLoops(stelzl[[1]])


###################################################
### chunk: abs
###################################################
abstract(stelzl1)


###################################################
### chunk: complexHG
###################################################
compXML <- file.path(fold, "intact_complexSample.xml")
pcHg <- buildPCHypergraph(compXML, INTACT.PSIMI25,
                          split.by = "organismName")
pcHg[[1]]


###################################################
### chunk: egHE
###################################################
complexes(pcHg[[1]])


###################################################
### chunk: activeBait
###################################################
deg = degree(stelzl1)
activeBait = names(which(deg$outDegree > 10 &
                           deg$outDegree<15))
proteins = union(activeBait, unlist(adj(stelzl1,
                                        activeBait)))
stelzlSG = subGraph(proteins, stelzl1)


###################################################
### chunk: graphAtt
###################################################
graph.par(list(nodes=list(fill="steelblue", label="")))
baitCol = rep("yellow", length(activeBait))
names(baitCol) = activeBait
nodeRenderInfo(stelzlSG) <- list(fill=baitCol)


###################################################
### chunk: graph-stelzlSG1
###################################################
stelzlSG <- layoutGraph(stelzlSG, layoutType="neato")
renderGraph(stelzlSG)


###################################################
### chunk: graph-stelzlSG2
###################################################
stelzlSG <- layoutGraph(stelzlSG, layoutType="twopi")
renderGraph(stelzlSG)



