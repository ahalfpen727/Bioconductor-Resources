### R code from vignette source 'gage.Rnw'

###################################################
### code chunk number 1: synopsis1 (eval = FALSE)
###################################################
## library(gage)
## data(gse16873)
## hn=(1:6)*2-1
## dcis=(1:6)*2


###################################################
### code chunk number 2: synopsis2 (eval = FALSE)
###################################################
## data(kegg.gs)
## gse16873.kegg.p <- gage(gse16873, gsets = kegg.gs, ref = hn, samp = dcis)


###################################################
### code chunk number 3: synopsis3 (eval = FALSE)
###################################################
## library(gageData)
## data(go.sets.hs)
## data(go.subs.hs)
## gse16873.bp.p <- gage(gse16873, gsets = go.sets.hs[go.subs.hs$BP], ref = hn, samp = dcis)
## gse16873.mf.p <- gage(gse16873, gsets = go.sets.hs[go.subs.hs$MF], ref = hn, samp = dcis)
## gse16873.cc.p <- gage(gse16873, gsets = go.sets.hs[go.subs.hs$CC], ref = hn, samp = dcis)


###################################################
### code chunk number 4: install (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite(c("gage","gageData"))


###################################################
### code chunk number 5: install (eval = FALSE)
###################################################
## install.packages(c("/your/local/directory/gage_2.9.1.tar.gz",
##     "/your/local/directory/gageData_1.0.0.tar.gz"), 
##     repos = NULL, type = "source")


###################################################
### code chunk number 6: start
###################################################
options(width=80)


###################################################
### code chunk number 7: start
###################################################
library(gage)


###################################################
### code chunk number 8: start (eval = FALSE)
###################################################
## library(help=gage)


###################################################
### code chunk number 9: start (eval = FALSE)
###################################################
## help(gage)
## ?gage


###################################################
### code chunk number 10: dataPrep
###################################################
data(gse16873)
cn=colnames(gse16873)
hn=grep('HN',cn, ignore.case =T)
adh=grep('ADH',cn, ignore.case =T)
dcis=grep('DCIS',cn, ignore.case =T)
print(hn)
print(dcis)


###################################################
### code chunk number 11: dataPrep
###################################################
data(kegg.gs)
data(go.gs)
lapply(kegg.gs[1:3],head)
head(rownames(gse16873))


###################################################
### code chunk number 12: kegg.gsets
###################################################
kg.hsa=kegg.gsets("hsa")
kegg.gs=kg.hsa$kg.sets[kg.hsa$sigmet.idx]
#save(kegg.gs, file="kegg.hsa.sigmet.gsets.RData")


###################################################
### code chunk number 13: korg
###################################################
#kegg.gsets works with 3000 KEGG species,for examples:
data(korg)
head(korg[,1:3])


###################################################
### code chunk number 14: go.gsets
###################################################
go.hs=go.gsets(species="human")
go.bp=go.hs$go.sets[go.hs$go.subs$BP]
go.mf=go.hs$go.sets[go.hs$go.subs$MF]
go.cc=go.hs$go.sets[go.hs$go.subs$CC]
#save(go.bp, go.mf, go.cc, file="go.hs.gsets.RData")


###################################################
### code chunk number 15: bods
###################################################
#for Bioconductor species supported by go.gsets function:
data(bods)
print(bods)


###################################################
### code chunk number 16: gage.1direction
###################################################
gse16873.kegg.p <- gage(gse16873, gsets = kegg.gs, 
    ref = hn, samp = dcis)
#go.gs here only the first 1000 entries as a fast example.
gse16873.go.p <- gage(gse16873, gsets = go.gs, 
    ref = hn, samp = dcis)
#or we can do analysis on BP, MF, or CC subcategories if we've
#generated the data as above.
#gse16873.bp.p <- gage(gse16873, gsets = go.bp,
#                        ref = hn, samp = dcis)
str(gse16873.kegg.p, strict.width='wrap')
head(gse16873.kegg.p$greater[, 1:5], 4)
head(gse16873.kegg.p$less[, 1:5], 4)
head(gse16873.kegg.p$stats[, 1:5], 4)


###################################################
### code chunk number 17: gage.2direction
###################################################
gse16873.kegg.2d.p <- gage(gse16873, gsets = kegg.gs, 
    ref = hn, samp = dcis, same.dir = F)
head(gse16873.kegg.2d.p$greater[,1:5], 4)
head(gse16873.kegg.2d.p$stats[,1:5], 4)


###################################################
### code chunk number 18: page
###################################################
gse16873.kegg.page.p <- gage(gse16873, gsets = kegg.gs, 
    ref = hn, samp = dcis, saaTest = gs.zTest)
head(gse16873.kegg.page.p$greater[, 1:5], 4)


###################################################
### code chunk number 19: gage.alternative (eval = FALSE)
###################################################
## gse16873.kegg.t.p <- gage(gse16873, gsets = kegg.gs, 
##     ref = hn, samp = dcis, use.fold = F)


###################################################
### code chunk number 20: gage.alternative (eval = FALSE)
###################################################
## gse16873.kegg.rk.p <- gage(gse16873, gsets = kegg.gs, 
##     ref = hn, samp = dcis, rank.test = T)


###################################################
### code chunk number 21: gage.alternative (eval = FALSE)
###################################################
## gse16873.kegg.ks.p <- gage(gse16873, gsets = kegg.gs, 
##     ref = hn, samp = dcis, saaTest = gs.KSTest)


###################################################
### code chunk number 22: gage.alternative (eval = FALSE)
###################################################
## gse16873.kegg.unpaired.p <- gage(gse16873, gsets = kegg.gs, 
##     ref = hn, samp = dcis, compare = "unpaired")


###################################################
### code chunk number 23: gage.alternative (eval = FALSE)
###################################################
## gse16873.kegg.gamma.p <- gage(gse16873, gsets = kegg.gs, 
##     ref = hn, samp = dcis, use.stouffer=FALSE)


###################################################
### code chunk number 24: symbol.ID
###################################################
head(rownames(gse16873))
gse16873.sym<-gse16873
data(egSymb)
rownames(gse16873.sym)<-eg2sym(rownames(gse16873.sym))
head(rownames(gse16873.sym))


###################################################
### code chunk number 25: symbol.ID
###################################################
kegg.gs.sym<-lapply(kegg.gs, eg2sym)
lapply(kegg.gs.sym[1:3],head)
gse16873.kegg.p2 <- gage(gse16873.sym, gsets = kegg.gs.sym, 
    ref = hn, samp = dcis)


###################################################
### code chunk number 26: full.table
###################################################
write.table(gse16873.kegg.2d.p$greater, file = "gse16873.kegg.2d.p.txt", 
    sep = "\t")
write.table(rbind(gse16873.kegg.p$greater, gse16873.kegg.p$less), 
    file = "gse16873.kegg.p.txt", sep = "\t")


###################################################
### code chunk number 27: significant.genesets
###################################################
gse16873.kegg.sig<-sigGeneSet(gse16873.kegg.p, outname="gse16873.kegg")
str(gse16873.kegg.sig, strict.width='wrap')
gse16873.kegg.2d.sig<-sigGeneSet(gse16873.kegg.2d.p, outname="gse16873.kegg")
str(gse16873.kegg.2d.sig, strict.width='wrap')
write.table(gse16873.kegg.2d.sig$greater, file = "gse16873.kegg.2d.sig.txt", 
    sep = "\t")
write.table(rbind(gse16873.kegg.sig$greater, gse16873.kegg.sig$less), 
    file = "gse16873.kegg.sig.txt", sep = "\t")


###################################################
### code chunk number 28: nonredundant.genesets
###################################################
gse16873.kegg.esg.up <- esset.grp(gse16873.kegg.p$greater, 
    gse16873, gsets = kegg.gs, ref = hn, samp = dcis, 
    test4up = T, output = T, outname = "gse16873.kegg.up", make.plot = F)
gse16873.kegg.esg.dn <- esset.grp(gse16873.kegg.p$less, 
    gse16873, gsets = kegg.gs, ref = hn, samp = dcis, 
    test4up = F, output = T, outname = "gse16873.kegg.dn", make.plot = F)
names(gse16873.kegg.esg.up)
head(gse16873.kegg.esg.up$essentialSets, 4)
head(gse16873.kegg.esg.up$setGroups, 4)
head(gse16873.kegg.esg.up$coreGeneSets, 4)


###################################################
### code chunk number 29: sel.gene.expression
###################################################
rownames(gse16873.kegg.p$greater)[1:3]
gs=unique(unlist(kegg.gs[rownames(gse16873.kegg.p$greater)[1:3]]))
essData=essGene(gs, gse16873, ref =hn, samp =dcis)
head(essData, 4)
ref1=1:6
samp1=7:12
for (gs in rownames(gse16873.kegg.p$greater)[1:3]) {
    outname = gsub(" |:|/", "_", substr(gs, 10, 100))
    geneData(genes = kegg.gs[[gs]], exprs = essData, ref = ref1, 
        samp = samp1, outname = outname, txt = T, heatmap = T, 
        Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}


###################################################
### code chunk number 30: all.gene.expression
###################################################
for (gs in rownames(gse16873.kegg.p$greater)[1:3]) {
    outname = gsub(" |:|/", "_", substr(gs, 10, 100))
    outname = paste(outname, "all", sep=".")
    geneData(genes = kegg.gs[[gs]], exprs = gse16873, ref = hn, 
        samp = dcis, outname = outname, txt = T, heatmap = T, 
        Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}


###################################################
### code chunk number 31: pathview
###################################################
library(pathview)
gse16873.d <- gse16873[ ,dcis] - gse16873[ ,hn]
path.ids=c("hsa04110 Cell cycle", "hsa00020 Citrate cycle (TCA cycle)")
path.ids2 <- substr(path.ids, 1, 8)
#native KEGG view
pv.out.list <- sapply(path.ids2, function(pid) pathview(gene.data = gse16873.d[, 
                        1:2], pathway.id = pid, species = "hsa"))
#Graphviz view
pv.out.list <- sapply(path.ids2, function(pid) pathview(gene.data = gse16873.d[, 
                        1:2], pathway.id = pid, species = "hsa", kegg.native=F,
                        sign.pos="bottomleft")) 


###################################################
### code chunk number 32: pathview.all (eval = FALSE)
###################################################
## sel <- gse16873.kegg.p$greater[, "q.val"] < 0.1 & !is.na(gse16873.kegg.p$greater[,
## "q.val"])
## path.ids <- rownames(gse16873.kegg.p$greater)[sel]
## path.ids2 <- substr(path.ids, 1, 8)
## pv.out.list <- sapply(path.ids2, function(pid) pathview(gene.data = gse16873.d[,
##                         1:2], pathway.id = pid, species = "hsa"))


###################################################
### code chunk number 33: gagePipe
###################################################
#introduce another half dataset
library(gageData)
data(gse16873.2)
cn2=colnames(gse16873.2)
hn2=grep('HN',cn2, ignore.case =T)
dcis2=grep('DCIS',cn2, ignore.case =T)
#batch GAGE analysis on the combined dataset
gse16873=cbind(gse16873, gse16873.2)
dataname='gse16873' #output data prefix
sampnames=c('dcis.1', 'dcis.2')
refList=list(hn, hn2+12)
sampList=list(dcis, dcis2+12)
gagePipe(gse16873, gsname = "kegg.gs", dataname = "gse16873", 
    sampnames = sampnames, ref.list = refList, samp.list = sampList, 
    comp.list = "paired")


###################################################
### code chunk number 34: gageComp
###################################################
load('gse16873.gage.RData')
gageComp(sampnames, dataname, gsname = "kegg.gs", 
    do.plot = TRUE)


###################################################
### code chunk number 35: heter.gage
###################################################
boxplot(data.frame(gse16873))
gse16873.kegg.heter.p <- heter.gage(gse16873, gsets = kegg.gs, 
    ref.list = refList, samp.list = sampList)
gse16873.kegg.heter.2d.p <- heter.gage(gse16873, gsets = kegg.gs, 
    ref.list = refList, samp.list = sampList, same.dir = F)


