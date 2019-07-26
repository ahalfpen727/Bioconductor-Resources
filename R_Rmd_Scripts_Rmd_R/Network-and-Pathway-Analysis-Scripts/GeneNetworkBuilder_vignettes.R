## ---- echo=FALSE, results="hide", warning=FALSE----------------------------
suppressPackageStartupMessages({
  library(GeneNetworkBuilder)
  library(Rgraphviz)
  library(XML)
})
knitr::opts_chunk$set(warning=FALSE, message=FALSE)

## ----quickStart------------------------------------------------------------
library(GeneNetworkBuilder)
##load C. elegans miRNA ID lists
data("ce.miRNA.map")
##load GNB-embedded regulatory network of C. elegans.
data("ce.interactionmap")
##load data required
data("example.data")
##build the network by binding list and interaction map
sifNetwork<-buildNetwork(TFbindingTable=example.data$ce.bind,
                        interactionmap=ce.interactionmap, level=2)
##filter the network by expression data
cifNetwork<-filterNetwork(rootgene="WBGene00000912", sifNetwork=sifNetwork,
                    exprsData=uniqueExprsData(example.data$ce.exprData),
                    mergeBy="symbols",
                    miRNAlist=as.character(ce.miRNA.map[ , 1]),
                    remove_miRNA=FALSE, tolerance=1)
##generate graphNEL object for the network
gR<-polishNetwork(cifNetwork=cifNetwork,
                  nodecolor=colorRampPalette(c("green",
                                               "yellow",
                                               "red"))(5))
##show network
(html <- browseNetwork(gR))
##save network
exportNetwork(html, "network.html")

## ----ce_example------------------------------------------------------------
library(GeneNetworkBuilder)
data("example.data")
##Initialize a binding matrix by TF and the related gene lists of TFBDs.
##For example, TF is daf-16, and the ChIP-chip result indicates that it can bind to
##upstream regions of gene "zip-2", "zip-4", "nhr-3" and "nhr-66".
bind<-cbind(from="daf-16", to=c("zip-2", "zip-4", "nhr-3", "nhr-66"))
##For same gene, there are multple gene alias. In order to eliminate the possibility of
##missing any interactions, convert the gene symbols to unique gene ids is important.
data("ce.IDsMap")
bind<-convertID(toupper(bind), IDsMap=ce.IDsMap, ByName=c("from", "to"))
##build the network by binding list and interaction map
data("ce.interactionmap")
sifNetwork<-buildNetwork(TFbindingTable=example.data$ce.bind,
                        interactionmap=ce.interactionmap, level=2)
##filter the network by expression data
##For each gene id, it should have only single record for expression change.
unique.ce.microarrayData<-uniqueExprsData(example.data$ce.exprData,
                        method="Max", condenseName='logFC')
data("ce.miRNA.map")
cifNetwork<-filterNetwork(rootgene="WBGene00000912", sifNetwork=sifNetwork,
                    exprsData=unique.ce.microarrayData, mergeBy="symbols",
                    miRNAlist=as.character(ce.miRNA.map[ , 1]),
                    tolerance=1, cutoffPVal=0.01, cutoffLFC=1)
##convert the unique gene ids back to gene symbols
data("ce.mapIDs")
cifNetwork<-convertID(cifNetwork, ce.mapIDs, ByName=c("from","to"))
##generate graphNEL object for the network
gR<-polishNetwork(cifNetwork, nodecolor=colorRampPalette(c("green", "yellow", "red"))(10))
##plot the figure
browseNetwork(gR)
## or plot by Rgraphviz
library(Rgraphviz)
plotNetwork<-function(gR, layouttype="dot", ...){
    if(!is(gR,"graphNEL")) stop("gR must be a graphNEL object")
    if(!(GeneNetworkBuilder:::inList(layouttype, c("dot", "neato", "twopi", "circo", "fdp")))){
        stop("layouttype must be dot, neato, twopi, circo or fdp")
    }
    g1<-Rgraphviz::layoutGraph(gR, layoutType=layouttype, ...)
    nodeRenderInfo(g1)$col <- nodeRenderInfo(gR)$col
    nodeRenderInfo(g1)$fill <- nodeRenderInfo(gR)$fill
    renderGraph(g1)
}
plotNetwork(gR)
##output file for cytoscape
exportNetwork(browseNetwork(gR),
              file="network.xgmml",
              format = "XGMML")
##output the GXL file
library("XML")
xml<-saveXML(toGXL(gR)$value())
z<-textConnection(xml)
cat(readLines(z, 8), sep="\n")
##calculate shortest path, ...
library(RBGL)
sp.between(gR,"daf-16","lam-2")


## ----hs_example------------------------------------------------------------
library(GeneNetworkBuilder)
data("hs.interactionmap")
data("hs.miRNA.map")
data("hs.IDsMap")
data("hs.mapIDs")
data("example.data")
rootgene<-"6657"
sifNetwork<-buildNetwork(example.data$hs.bind, hs.interactionmap, level=3)
##example.data$ce$exprData is the combination of gene/miRNA expression profile
##note, here should set the miRNAtol to TRUE
cifNetwork<-filterNetwork(rootgene=rootgene, sifNetwork=sifNetwork,
                   exprsData=example.data$hs.exprData, mergeBy="symbols",
                   miRNAlist=as.character(hs.miRNA.map[,1]),
                   tolerance=0, miRNAtol=TRUE)
cifNetwork<-convertID(cifNetwork, hs.mapIDs, ByName=c("from","to"))
gR<-polishNetwork(cifNetwork)
##plot the figure
browseNetwork(gR)
## or plot by Rcyjs
# library(RCyjs)
# rcy <- RCyjs(portRange=9047:9067, title='sox2', graph=gR)
# setBrowserWindowTitle(rcy, "sox2")
# setNodeLabelRule(rcy, "label")
# setNodeLabelAlignment(rcy, "center", "center")
# setDefaultNodeColor(rcy, "white")
# setDefaultNodeBorderColor(rcy, "black")
# setDefaultNodeBorderWidth(rcy, 1)
# logFC.range <- range(cifNetwork$logFC)
# logFC.range <- max(abs(logFC.range))
# setNodeColorRule(rcy, "logFC", c(-logFC.range, 0, logFC.range), c("green", "yellow", "red"), mode="interpolate")
# size.range <- range(noa(gR, "size"))
# setNodeSizeRule(rcy, "size", size.range, size.range)
# layout(rcy, "cose")
# redraw(rcy)

## ----sessionInfo-----------------------------------------------------------
sessionInfo()

