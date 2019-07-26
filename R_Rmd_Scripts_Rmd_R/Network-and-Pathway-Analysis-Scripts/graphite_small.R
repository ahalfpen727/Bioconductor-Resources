### R code from vignette source 'graphite.Rnw'
library(graph)
library(graphite)
pathwayDatabases()

humanReactome <- pathways("hsapiens", "reactome")
names(humanReactome)[1:10]
p <- humanReactome[["ABC-family proteins mediated transport"]]
p

p <- humanReactome[[21]]
p@title

head(nodes(p))
head(edges(p))
pathwayDatabases()
g <- pathwayGraph(p)
g


###################################################
### code chunk number 9: graph2
###################################################
edgeData(g)[1]


###################################################
### code chunk number 10: ident1
###################################################
pSymbol <- convertIdentifiers(p, "SYMBOL")
pSymbol
head(nodes(pSymbol))


###################################################
### code chunk number 11: ident2
###################################################
reactomeSymbol <- convertIdentifiers(humanReactome[1:5], "SYMBOL")


###################################################
### code chunk number 12: spia1
###################################################
library(SPIA)
data(colorectalcancer)

library(hgu133plus2.db)
x <- hgu133plus2ENTREZID
top$ENTREZ <- unlist(as.list(x[top$ID]))
top <- top[!is.na(top$ENTREZ), ]
top <- top[!duplicated(top$ENTREZ), ]
tg1 <- top[top$adj.P.Val < 0.05, ]

DE_Colorectal = tg1$logFC
names(DE_Colorectal) <- as.vector(tg1$ENTREZ)
ALL_Colorectal <- top$ENTREZ

prepareSPIA(biocarta[1:2], "biocartaEx")
res <- runSPIA(de=DE_Colorectal, all=ALL_Colorectal, "biocartaEx")
res


###################################################
### code chunk number 13: tg1
###################################################
library(topologyGSA)
data(examples)

kegg <- pathways("hsapiens", "kegg")
p <- convertIdentifiers(kegg[["Fc epsilon RI signaling pathway"]], "SYMBOL")
runTopologyGSA(p, "var", y1, y2, 0.05)


###################################################
### code chunk number 14: tg2
###################################################
library(clipper)
library(ALL)

kegg <- pathways("hsapiens", "kegg")
paths <-convertIdentifiers(kegg[1:5], "ENTREZID")
genes <- unlist(lapply(paths, nodes))
data(ALL)
all <- as.matrix(exprs(ALL[1:length(genes),1:20]))
classes <- c(rep(1,10), rep(2,10))
rownames(all) <- genes
clipped <- runClipper(paths, all, classes, "mean", pathThr=0.1)
resClip <- do.call(rbind,clipped$results)
resClip[,1:5]


