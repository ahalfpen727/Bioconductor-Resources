### R code from vignette source 'layingOutPathways.Rnw'

###################################################
### code chunk number 1: setup
###################################################
library("Rgraphviz")
data("integrinMediatedCellAdhesion")
IMCAGraph


###################################################
### code chunk number 2: initial
###################################################
attrs <- list(graph=list(rankdir="LR"))
IMCAGraph <- layoutGraph(IMCAGraph, attrs=attrs)
renderGraph(IMCAGraph)


###################################################
### code chunk number 3: longLabel
###################################################
n <- nodes(IMCAGraph)
names(labels) <- labels <- n
nc <- nchar(labels)
table(nc)
long <- labels[order(nc, decreasing=TRUE)][1:4]
long


###################################################
### code chunk number 4: linefeed
###################################################
labels[long] <- c(paste("Phosphatidyl-\ninositol\n",
    "signaling\nsystem", sep=""), "cell\nproliferation", 
    "cell\nmaintenance", "cell\nmotility") 


###################################################
### code chunk number 5: setLabel
###################################################
nodeRenderInfo(IMCAGraph) <- list(label=labels)
renderGraph(IMCAGraph)


###################################################
### code chunk number 6: redoLayout1
###################################################
attrs$node <- list(fixedsize=FALSE)
width <- c(2.5, 1.5, 1.5, 1.5)
height <- c(1.5, 1.5, 1.5, 1.5)
names(width) <- names(height) <- long
nodeAttrs <- list(width=width, height=height)
IMCAGraph <- layoutGraph(IMCAGraph, attrs=attrs, 
                         nodeAttrs=nodeAttrs)
renderGraph(IMCAGraph)


###################################################
### code chunk number 7: redoLayout2
###################################################
shape <- rep("rectangle", length(n))
names(shape) <- n
shape[long[1]] <- "ellipse" 
shape[c(long[2:4], "F-actin")] <- "plaintext"
nodeRenderInfo(IMCAGraph) <- list(shape=shape)
IMCAGraph <- layoutGraph(IMCAGraph, attrs=attrs, 
                         nodeAttrs=nodeAttrs)
renderGraph(IMCAGraph)


###################################################
### code chunk number 8: colorPlot
###################################################
colors <- rep("lightgreen", length(n))
names(colors) <- n
transp <- c("ITGB", "ITGA", "MYO", "ACTN", "JNK", "p110", 
            "Phosphatidylinositol signaling system", 
            "PI5K", "MYO-P", "cell maintenance", "cell motility", 
            "F-actin", "cell proliferation")
colors[transp] <- "transparent"
nodeRenderInfo(IMCAGraph) <- list(fill=colors)
renderGraph(IMCAGraph)


###################################################
### code chunk number 9: subgraphs
###################################################
sg1 <- subGraph(c("ITGA", "ITGB", "ILK", "CAV"), IMCAGraph)
sg2 <- subGraph(c("cell maintenance", "cell motility", 
                  "F-actin", "cell proliferation"), IMCAGraph)
sg3 <- subGraph(c("ACTN", "VCL", "TLN", "PXN", "TNS", "VASP"), 
                IMCAGraph)
sg4 <- subGraph(setdiff(n, c(nodes(sg1), nodes(sg2), nodes(sg3))), 
                IMCAGraph)


###################################################
### code chunk number 10: subGList
###################################################
subGList <- vector(mode="list", length=4)
subGList[[1]] <- list(graph=sg1, attrs=c(rank="source"))
subGList[[2]] <- list(graph=sg2, attrs=c(rank="sink"))
subGList[[3]] <- list(graph=sg3, cluster=TRUE)
subGList[[4]] <- list(graph=sg3, cluster=TRUE)


###################################################
### code chunk number 11: plotSubgraph
###################################################
IMCAGraph <- layoutGraph(IMCAGraph, attrs=attrs, 
                         nodeAttrs=nodeAttrs, subGList=subGList)
renderGraph(IMCAGraph)


###################################################
### code chunk number 12: expressionGraph
###################################################
require("geneplotter")
require("fibroEset")
require("hgu95av2.db")
data("fibroEset")
plotExpressionGraph(IMCAGraph, IMCAAttrs$LocusLink,
                    exprs(fibroEset)[,1], hgu95av2ENTREZID, 
                    attrs=attrs,
                    subGList=subGList, nodeAttr=nodeAttrs)


###################################################
### code chunk number 13: VJCGraph
###################################################
z <- IMCAGraph
nodeRenderInfo(z) <- list(shape="plaintext", fontsize=100)
nag <- layoutGraph(z, attrs=list(edge=list(arrowsize=2.8, minlen=3)))
renderGraph(nag)


###################################################
### code chunk number 14: layingOutPathways.Rnw:313-314
###################################################
sessionInfo()


###################################################
### code chunk number 15: layingOutPathways.Rnw:318-319
###################################################
graphvizVersion()


