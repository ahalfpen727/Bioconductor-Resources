### R code from vignette source 'biocGraph.Rnw'

###################################################
### code chunk number 1: createGraph1
###################################################
library("biocGraph")
set.seed(123)
V <- letters[1:10]
M <- 1:4


###################################################
### code chunk number 2: createGraph2
###################################################
g1 <- randomGraph(V, M, 0.2)


###################################################
### code chunk number 3: plotDot
###################################################
plot(g1)


###################################################
### code chunk number 4: plotNeato
###################################################
plot(g1, "neato")


###################################################
### code chunk number 5: plotTwopi
###################################################
plot(g1, "twopi")


###################################################
### code chunk number 6: rEG
###################################################
rEG <- new("graphNEL", nodes=c("A", "B"), edgemode="directed")
rEG <- addEdge("A", "B", rEG, 1)
rEG <- addEdge("B", "A", rEG, 1)


###################################################
### code chunk number 7: recipEdgesComb
###################################################
plot(rEG)


###################################################
### code chunk number 8: recipEdgesDistinct
###################################################
plot(rEG, recipEdges="distinct")


###################################################
### code chunk number 9: removedEdges
###################################################
removedEdges(g1)


###################################################
### code chunk number 10: getSubgraphs
###################################################
sg1 <- subGraph(c("a","d","j","i"), g1)
sg1
sg2 <- subGraph(c("b","e","h"), g1)
sg3 <- subGraph(c("c","f","g"), g1)


###################################################
### code chunk number 11: subGplot
###################################################
subGList <- vector(mode="list", length=3)
subGList[[1]] <- list(graph=sg1)
subGList[[2]] <- list(graph=sg2, cluster=FALSE)
subGList[[3]] <- list(graph=sg3)
plot(g1, subGList=subGList)


###################################################
### code chunk number 12: subGPlot2
###################################################
sg1 <- subGraph(c("a","c","d","e","j"), g1)
sg2 <- subGraph(c("f","h","i"), g1)
plot(g1, subGList=list(list(graph=sg1), list(graph=sg2)))


###################################################
### code chunk number 13: edgeNames
###################################################
edgeNames(g1)
edgeNames(g1, recipEdges="distinct")


###################################################
### code chunk number 14: defAttrs
###################################################
defAttrs <- getDefaultAttrs()


###################################################
### code chunk number 15: defAttrs2
###################################################
plot(g1, attrs=list(node=list(label="foo", fillcolor="lightgreen"),
           edge=list(color="cyan"),
           graph=list(rankdir="LR")))


###################################################
### code chunk number 16: baseLists
###################################################
nAttrs <- list()
eAttrs <- list()


###################################################
### code chunk number 17: makeLabels1
###################################################
z <- strsplit(packageDescription("Rgraphviz")$Description, " ")[[1]]
z <- z[1:numNodes(g1)]
names(z) = nodes(g1)
nAttrs$label <- z


###################################################
### code chunk number 18: makeLabels2
###################################################
eAttrs$label <- c("a~h"="Label 1", "c~h"="Label 2")


###################################################
### code chunk number 19: makeLabels3
###################################################
attrs <- list(node=list(shape="ellipse", fixedsize=FALSE))


###################################################
### code chunk number 20: figLabels
###################################################
plot(g1, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs)


###################################################
### code chunk number 21: edgeWeights
###################################################
ew <- as.character(unlist(edgeWeights(g1)))
ew <- ew[setdiff(seq(along=ew), removedEdges(g1))]
names(ew) <- edgeNames(g1)
eAttrs$label <- ew
## FIXME (wh 17.6.06): This does not work - see bug report:
## attrs$edge$labelfontsize="27"


###################################################
### code chunk number 22: edgeWeightLabels
###################################################
plot(g1, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs)


###################################################
### code chunk number 23: colors
###################################################
## Specify node drawing color
nAttrs$color <- c(a="red", b="red", g="green", d="blue")

## Specify edge drawing color
eAttrs$color <- c("a~d"="blue", "c~h"="purple")

## Specify node fill color
nAttrs$fillcolor <- c(j="yellow")

## label color
nAttrs$fontcolor <- c(e="green", f="red")
eAttrs$fontcolor <- c("a~h"="green", "c~h"="brown")

nAttrs
eAttrs


###################################################
### code chunk number 24: figColors
###################################################
plot(g1, nodeAttrs=nAttrs, attrs=attrs)


###################################################
### code chunk number 25: nodeShapes
###################################################
attrs$node$shape <- "ellipse"
nAttrs$shape <- c(g="box", f="circle", j="box", a="plaintext")


###################################################
### code chunk number 26: figNodeShapes
###################################################
plot(g1, attrs=attrs, nodeAttrs=nAttrs)


###################################################
### code chunk number 27: getLists1
###################################################
nodes <- buildNodeList(g1)
edges <- buildEdgeList(g1)


###################################################
### code chunk number 28: getLists2
###################################################
nodes[[1]]
edges[[1]]


###################################################
### code chunk number 29: buildwithAttrs
###################################################
nodes <- buildNodeList(g1, nodeAttrs=nAttrs, defAttrs=defAttrs$node)
edges <- buildEdgeList(g1, edgeAttrs=eAttrs, defAttrs=defAttrs$edge)
nodes[[1]]
edges[[1]]


###################################################
### code chunk number 30: arrowheads
###################################################
for(j in c("a~e", "a~h"))
  edges[[j]]@attrs$arrowhead <- "open"


###################################################
### code chunk number 31: plotbuild
###################################################
vv <- agopen(name="foo", nodes=nodes, edges=edges, attrs=attrs, 
  edgeMode="undirected")
plot(vv)


###################################################
### code chunk number 32: graph17
###################################################
data(graphExamples)
z <- graphExamples[[8]]
nNodes <- length(nodes(z))

nA <- list()
nA$fixedSize<-rep(FALSE, nNodes)
nA$height <- nA$width <- rep("1", nNodes)
nA$label <- rep("z", nNodes)
nA$color <- rep("green", nNodes)
nA$fillcolor <- rep("orange", nNodes)
nA$shape <- rep("circle", nNodes)
nA$fontcolor <- rep("blue", nNodes)
nA$fontsize <- rep(14, nNodes)
nA <- lapply(nA, function(x) { names(x) <- nodes(z); x})


###################################################
### code chunk number 33: graph17
###################################################
plot(z, nodeAttrs=nA)


###################################################
### code chunk number 34: pieChartCalc
###################################################
set.seed(123)
counts = matrix(rexp(numNodes(g1)*4), ncol=4)

g1layout <- agopen(g1, name="foo")

makeNodeDrawFunction <- function(x) {
  force(x)
  function(node, ur, attrs, radConv) {
    nc <- getNodeCenter(node)
    pieGlyph(x, 
             xpos=getX(nc),
             ypos=getY(nc),
             radius=getNodeRW(node),
             col=rainbow(4))
    text(getX(nc), getY(nc), paste(signif(sum(x), 2)), 
      cex=0.5, col="white", font=2)
  }
}

drawFuns <- apply(counts, 1, makeNodeDrawFunction)


###################################################
### code chunk number 35: pieChartGraph
###################################################
plot(g1layout, drawNode=drawFuns, main="Example Pie Chart Plot")


###################################################
### code chunk number 36: clusterGraph1
###################################################
cG <- new("clusterGraph", clusters=list(a=c(1:10), b=c(11:13),
                            c=c(14:20), d=c(21, 22)))


###################################################
### code chunk number 37: cGdot
###################################################
plot(cG, main="dot")


###################################################
### code chunk number 38: cGtwopi
###################################################
par(bg="#e0e0e0")
plot(cG, "twopi", main="twopi")


###################################################
### code chunk number 39: cGneato
###################################################
plot(cG, "neato", main="neato")


###################################################
### code chunk number 40: bipartite1
###################################################
set.seed(123)
nodes1 <- paste(0:7)
nodes2 <- letters[1:10]

ft <- cbind(sample(nodes1, 24, replace=TRUE), 
            sample(nodes2, 24, replace=TRUE))
ft <- ft[!duplicated(apply(ft, 1, paste, collapse="")),]

g <-  ftM2graphNEL(ft, edgemode='directed')
g


###################################################
### code chunk number 41: bipartitelayout
###################################################
twocolors <- c("#D9EF8B", "#E0F3F8")
nodeType <- 1 + (nodes(g) %in% nodes1)
nA = makeNodeAttrs(g, fillcolor=twocolors[nodeType])

sg1 = subGraph(nodes1, g)
sgL = list(list(graph=sg1, cluster = FALSE, attrs = c(rank="sink")))

att = list(graph = list(rankdir = "LR", rank = ""))


###################################################
### code chunk number 42: figbipartite
###################################################
plot(g, attrs = att, nodeAttrs=nA, subGList = sgL)


###################################################
### code chunk number 43: agopenSimpleDemo
###################################################
library("graph")
g1_gz <- gzfile(system.file("GXL/graphExample-01.gxl.gz",package="graph"), open="rb")
g11_gz <- gzfile(system.file("GXL/graphExample-11.gxl.gz",package="graph"), open="rb")
g1 <- fromGXL(g1_gz)
g11 <- fromGXL(g11_gz)
g1_11 <- join(g1, g11)
sgl <- vector(mode="list", length=2)
sgl[[1]] <- list(graph=g1, cluster=FALSE)
sgl[[2]] <- list(graph=g11, cluster=TRUE)
ng <- agopenSimple(g1_11, "tmpsg", subGList=sgl)


###################################################
### code chunk number 44: DataDefaultsDemo1
###################################################
graphDataDefaults(ng)
nodeDataDefaults(ng)
edgeDataDefaults(ng)


###################################################
### code chunk number 45: DataDefaultsDemo2
###################################################
graphDataDefaults(ng, c("size", "bgcolor")) <- c("1", "yellow")
nodeDataDefaults(ng, c("fontcolor", "width")) <- c("blue", 0.5)
edgeDataDefaults(ng, c("color", "style")) <- c("green", "dotted")


###################################################
### code chunk number 46: DataDemo1
###################################################
graphData(ng, "bgcolor")
nodeData(ng, "a", c("fontcolor", "width"))
edgeData(ng, "f", "h", c("color", "arrowhead"))



###################################################
### code chunk number 47: DataDemo2
###################################################
graphData(ng, "bgcolor") <- "orange"
clusterData(ng, 2, "bgcolor") <- "red"
nodeData(ng, "a", c("fontcolor", "width")) <- c("red", "0.8")
edgeData(ng, "f", "h", c("color", "style")) <- c("blue", "solid")


###################################################
### code chunk number 48: layoutRenderDemo1
###################################################
plot(ng, "neato")
plot(ng, "circo")


###################################################
### code chunk number 49: layoutRenderDemo1
###################################################
toFile(ng, layoutType="dot", filename="test_dot.svg", fileType="svg")
toFile(ng, layoutType="circo", filename="test_circo.ps", fileType="ps")
toFile(ng, layoutType="twopi", filename="test_twopi.dot", fileType="dot")


###################################################
### code chunk number 50: integrinMediatedCellAdhesion
###################################################
data("integrinMediatedCellAdhesion")


###################################################
### code chunk number 51: dummyPlots
###################################################
outdir=tempdir()
nd = nodes(IMCAGraph)
plotFiles = paste(seq(along=nd), 'png', sep='.')
for(i in seq(along=nd)) {
  png(file.path(outdir, plotFiles[i]), width=400, height=600)
  plot(cumsum(rnorm(100)), type='l', col='blue', lwd=2, main=nd[i])
  dev.off()
}


###################################################
### code chunk number 52: indexhtml
###################################################
fhtml = file.path(outdir, "index.html")
con = file(fhtml, open="wt")
cat("<HTML><HEAD><TITLE>",
"Integrin Mediated Cell Adhesion graph with tooltips and hyperlinks, please click on the nodes.",
"</TITLE></HEAD>",
"<FRAMESET COLS=\"3*,2*\" BORDER=0>",
"  <FRAME SRC=\"graph.html\">",
"  <FRAME NAME=\"nodedata\">",
"</FRAMESET></HTML>", sep="\n", file=con)
close(con)


###################################################
### code chunk number 53: makegraph (eval = FALSE)
###################################################
## width = 600
## height = 512
## imgname = "graph.png"
## png(file.path(outdir, imgname), width=width, height=height)
## par(mai=rep(0,4))
## 
## lg = agopen(IMCAGraph, name="foo", 
##   attrs = list(graph=list(rankdir="LR", rank=""), node=list(fixedsize=FALSE)), 
##   nodeAttrs = makeNodeAttrs(IMCAGraph), 
##   subGList = IMCAAttrs$subGList)
## plot(lg)    
##  
## con = file(file.path(outdir, "graph.html"), open="wt")
## writeLines("<html><body>\n", con)
## imageMap(lg, con=con,
##          tags=list(HREF=plotFiles,
##            TITLE = nd,
##            TARGET = rep("nodedata", length(nd))),
##          imgname=imgname, width=width, height=height)
## writeLines("</body></html>", con)
## close(con)
## dev.off()


###################################################
### code chunk number 54: browseurl
###################################################
fhtml
if(interactive())
  browseURL(fhtml)


###################################################
### code chunk number 55: biocGraph.Rnw:987-988
###################################################
sessionInfo()


