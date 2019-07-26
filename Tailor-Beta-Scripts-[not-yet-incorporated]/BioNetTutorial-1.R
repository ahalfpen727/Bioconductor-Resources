### R code from vignette source 'Tutorial.Rnw'

###################################################
### code chunk number 1: Tutorial.Rnw:50-52
###################################################
options(width=60)
ps.options(family="sans")


###################################################
### code chunk number 2: Tutorial.Rnw:55-59
###################################################
library(BioNet)
library(DLBCL)
data(dataLym)
data(interactome)


###################################################
### code chunk number 3: Tutorial.Rnw:63-66
###################################################
pvals <- cbind(t=dataLym$t.pval, s=dataLym$s.pval)
rownames(pvals) <- dataLym$label
pval <- aggrPvals(pvals, order=2, plot=FALSE)


###################################################
### code chunk number 4: Tutorial.Rnw:70-73
###################################################
subnet <- subNetwork(dataLym$label, interactome)
subnet <- rmSelfLoops(subnet)
subnet


###################################################
### code chunk number 5: Tutorial.Rnw:79-81
###################################################
fb <- fitBumModel(pval, plot=FALSE)
scores <- scoreNodes(subnet, fb, fdr=0.001)


###################################################
### code chunk number 6: Tutorial.Rnw:88-91
###################################################
module <- runFastHeinz(subnet, scores)
logFC <- dataLym$diff
names(logFC) <- dataLym$label


###################################################
### code chunk number 7: Tutorial.Rnw:98-99
###################################################
plotModule(module, scores=scores, diff.expr=logFC)


###################################################
### code chunk number 8: Tutorial.Rnw:154-158
###################################################
library(BioNet)
library(DLBCL)
data(exprLym)
data(interactome)


###################################################
### code chunk number 9: Tutorial.Rnw:166-167
###################################################
exprLym


###################################################
### code chunk number 10: Tutorial.Rnw:173-174
###################################################
interactome


###################################################
### code chunk number 11: Tutorial.Rnw:181-183
###################################################
network <- subNetwork(featureNames(exprLym), interactome)
network


###################################################
### code chunk number 12: Tutorial.Rnw:188-190
###################################################
network <- largestComp(network)
network


###################################################
### code chunk number 13: Tutorial.Rnw:200-204
###################################################
library(genefilter)
library(impute)
expressions <- impute.knn(exprs(exprLym))$data
t.test <- rowttests(expressions, fac=exprLym$Subgroup)


###################################################
### code chunk number 14: Tutorial.Rnw:207-208
###################################################
t.test[1:10, ]


###################################################
### code chunk number 15: Tutorial.Rnw:215-218
###################################################
library(xtable)
top.table <- xtable(t.test[1:10,], display=c("s", "f", "f", "f"))
print(top.table, floating=FALSE)


###################################################
### code chunk number 16: Tutorial.Rnw:228-233
###################################################
data(dataLym)
ttest.pval <- t.test[, "p.value"]
surv.pval <- dataLym$s.pval
names(surv.pval) <- dataLym$label
pvals <- cbind(ttest.pval, surv.pval)


###################################################
### code chunk number 17: Tutorial.Rnw:242-243
###################################################
pval <- aggrPvals(pvals, order=2, plot=FALSE)


###################################################
### code chunk number 18: Tutorial.Rnw:250-252
###################################################
fb <- fitBumModel(pval, plot=FALSE)
fb


###################################################
### code chunk number 19: Tutorial.Rnw:255-260
###################################################
dev.new(width=13, height=7)
par(mfrow=c(1,2))
hist(fb)
plot(fb)
dev.off()


###################################################
### code chunk number 20: Tutorial.Rnw:275-276
###################################################
plotLLSurface(pval, fb)


###################################################
### code chunk number 21: Tutorial.Rnw:285-286
###################################################
scores <- scoreNodes(network=network, fb=fb, fdr=0.001)


###################################################
### code chunk number 22: Tutorial.Rnw:293-296
###################################################
network <- rmSelfLoops(network)
writeHeinzEdges(network=network, file="lymphoma_edges_001", use.score=FALSE)
writeHeinzNodes(network=network, file="lymphoma_nodes_001", node.scores = scores)


###################################################
### code chunk number 23: Tutorial.Rnw:313-315
###################################################
datadir <- file.path(path.package("BioNet"), "extdata")
dir(datadir)


###################################################
### code chunk number 24: Tutorial.Rnw:321-324
###################################################
module <- readHeinzGraph(node.file=file.path(datadir, "lymphoma_nodes_001.txt.0.hnz"), network=network)
diff <- t.test[, "dm"]
names(diff) <- rownames(t.test)


###################################################
### code chunk number 25: Tutorial.Rnw:327-328
###################################################
plotModule(module, diff.expr=diff, scores=scores)


###################################################
### code chunk number 26: Tutorial.Rnw:350-353
###################################################
sum(scores[nodes(module)])
sum(scores[nodes(module)]>0)
sum(scores[nodes(module)]<0)


###################################################
### code chunk number 27: Tutorial.Rnw:375-380
###################################################
library(BioNet)
library(DLBCL)
library(ALL)
data(ALL)
data(interactome)


###################################################
### code chunk number 28: Tutorial.Rnw:387-388
###################################################
ALL


###################################################
### code chunk number 29: Tutorial.Rnw:396-397
###################################################
interactome


###################################################
### code chunk number 30: Tutorial.Rnw:407-409
###################################################
mapped.eset <- mapByVar(ALL, network=interactome, attr="geneID")
mapped.eset[1:5,1:5]


###################################################
### code chunk number 31: Tutorial.Rnw:414-415
###################################################
length(intersect(rownames(mapped.eset), nodes(interactome)))


###################################################
### code chunk number 32: Tutorial.Rnw:424-429
###################################################
network <- subNetwork(rownames(mapped.eset), interactome)
network
network <- largestComp(network)
network <- rmSelfLoops(network)
network


###################################################
### code chunk number 33: Tutorial.Rnw:441-449
###################################################
library(limma)
design <- model.matrix(~ -1+ factor(c(substr(unlist(ALL$BT), 0, 1))))
colnames(design)<- c("B", "T")
contrast.matrix <- makeContrasts(B-T, levels=design)
contrast.matrix
fit <- lmFit(mapped.eset, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


###################################################
### code chunk number 34: Tutorial.Rnw:454-455
###################################################
pval <- fit2$p.value[,1]


###################################################
### code chunk number 35: Tutorial.Rnw:466-468
###################################################
fb <- fitBumModel(pval, plot=FALSE)
fb


###################################################
### code chunk number 36: Tutorial.Rnw:471-475
###################################################
dev.new(width=13, height=7)
par(mfrow=c(1,2))
hist(fb)
plot(fb)


###################################################
### code chunk number 37: Tutorial.Rnw:487-488
###################################################
scores <- scoreNodes(network=network, fb=fb, fdr=1e-14)


###################################################
### code chunk number 38: Tutorial.Rnw:494-496
###################################################
writeHeinzEdges(network=network, file="ALL_edges_001", use.score=FALSE)
writeHeinzNodes(network=network, file="ALL_nodes_001", node.scores = scores)


###################################################
### code chunk number 39: Tutorial.Rnw:517-519
###################################################
datadir <- file.path(path.package("BioNet"), "extdata")
module <- readHeinzGraph(node.file=file.path(datadir, "ALL_nodes_001.txt.0.hnz"), network=network)


###################################################
### code chunk number 40: Tutorial.Rnw:524-529
###################################################
nodeDataDefaults(module, attr="diff") <- ""
nodeData(module, n=nodes(module), attr="diff") <- fit2$coefficients[nodes(module),1]
nodeDataDefaults(module, attr="score") <- ""
nodeData(module, n=nodes(module), attr="score") <- scores[nodes(module)]
nodeData(module)[1]


###################################################
### code chunk number 41: Tutorial.Rnw:535-536
###################################################
saveNetwork(module, file="ALL_module", type="XGMML")


###################################################
### code chunk number 42: Tutorial.Rnw:568-576 (eval = FALSE)
###################################################
## j.repl <- 100
## resampling.pvals <- list()
## for(i in 1:j.repl)
## {
##   resampling.result <- resamplingPvalues(exprMat=mapped.eset, groups=factor(c(substr(unlist(ALL$BT), 0, 1))), resampleMat=FALSE, alternative="two.sided")
##   resampling.pvals[[i]] <- resampling.result$p.values
##   print(i)
## }


###################################################
### code chunk number 43: Tutorial.Rnw:585-591 (eval = FALSE)
###################################################
## fb <- lapply(resampling.pvals, fitBumModel, plot=FALSE, starts=1)
## resampling.scores <- c()
## for(i in 1:j.repl)
## {
##   resampling.scores[[i]] <- scoreNodes(network=network, fb=fb[[i]], fdr=1e-14)
## }


###################################################
### code chunk number 44: Tutorial.Rnw:596-598 (eval = FALSE)
###################################################
## score.mat <- as.data.frame(resampling.scores)
## colnames(score.mat) <- paste("resample", (1:j.repl), sep="")


###################################################
### code chunk number 45: Tutorial.Rnw:603-605 (eval = FALSE)
###################################################
## writeHeinzEdges(network=network, file="ALL_e_resample", use.score=FALSE)
## writeHeinzNodes(network=network, file="ALL_n_resample", node.scores = score.mat)


###################################################
### code chunk number 46: Tutorial.Rnw:618-620
###################################################
datadir <- file.path(path.package("BioNet"), "extdata")
modules <- readHeinzGraph(node.file=file.path(datadir, "ALL_n_resample.txt.0.hnz"), network=network)


###################################################
### code chunk number 47: Tutorial.Rnw:629-631
###################################################
cons.scores <- consensusScores(modules, network)
writeHeinz(network=network, file="ALL_cons", node.scores=cons.scores$N.scores, edge.scores=cons.scores$E.scores)


###################################################
### code chunk number 48: Tutorial.Rnw:638-643
###################################################
datadir <- file.path(path.package("BioNet"), "extdata")
cons.module <- readHeinzGraph(node.file=file.path(datadir, "ALL_cons_n.txt.0.hnz"), network=network)
cons.edges <- sortedEdgeList(cons.module)
E.width <- 1+cons.scores$E.freq[cons.edges]*10
N.size <- 1+cons.scores$N.freq[nodes(cons.module)]*10


###################################################
### code chunk number 49: Tutorial.Rnw:647-648
###################################################
plotModule(cons.module, edge.width=E.width, vertex.size=N.size, edge.label=cons.scores$E.freq[cons.edges]*100, edge.label.cex=0.6)


