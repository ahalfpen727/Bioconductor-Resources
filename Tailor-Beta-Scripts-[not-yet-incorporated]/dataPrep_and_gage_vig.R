### R code from vignette source 'vignettes/gage/inst/doc/dataPrep.Rnw'
devtools::install_github("egeulgen/pathfindR")
source("http://bioconductor.org/biocLite.R")
biocLite(c("RITAN","RITANdata", ""))
install.packages("RITAN")
library(devtools)
devtools::install_github('oganm/MSigDB')
###################################################
### code chunk number 1: start
###################################################
Sys.getenv()
options(width=80)
###################################################
### code chunk number 2: start
###################################################
library(gage)
library(gageData)


data(kegg.gs)
data(go.gs)
data(carta.gs)

data(go.sets.hs)
data(go.subs.hs)
data(kegg.gs.dise)

data(korg)
data(bods)
data(egSymb)
#kegg.gsets works with 3000 KEGG species,for examples:


#check the data
head(korg[,1:3])
head(demo.data)
str(demo.data)
#convert the data.frame into a matrix as to speed up the computing
demo.data=as.matrix(demo.data)
str(demo.data)

###################################################
# set up kegg database
kg.hsa=kegg.gsets("hsa")
save(kg.hsa, file="kegg.hsa.sigmet.gsets.RData")
kegg.sigmet.gs <- kg.hsa$kg.sets[kg.hsa$sigmet.idx]
kegg.dise.gs <- kg.hsa$kg.sets[kg.hsa$dise.idx]


# set up go database
go.hs <- go.gsets(species="human",id.type = "eg")
go.hs <- go.gsets(species="human")

go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]
go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF]
go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC]

###################################################
### code chunk number 3: demo.data
###################################################
###################################################
### code chunk number 4: readList
###################################################
#an example GMT gene set data derived from MSigDB data
filename=system.file("extdata/c2.demo.gmt", package = "gage")
demo.gs=readList(filename)
demo.gs[1:3]
#to use these gene sets with gse16873, need to convert the gene symbols
#to Entrez IDs first
data(egSymb)
demo.gs.sym<-lapply(demo.gs, sym2eg)
demo.gs.sym[1:3]

data(hnrnp.cnts)
cnts=hnrnp.cnts
dim(cnts)

sel.rn=rowSums(cnts) != 0
cnts=cnts[sel.rn,]
dim(cnts)
###################################################
### code chunk number 5: gse16873.affyid
###################################################
library(gageData)
data(gse16873.affyid)
affyid=rownames(gse16873.affyid)

library(hgu133a.db)
egids2=hgu133aENTREZID[affyid]
annots=toTable(egids2)
str(annots)
gse16873.affyid=gse16873.affyid[annots$probe_id,]

#if multiple probe sets map to a gene, select the one with maximal IQR
iqrs=apply(gse16873.affyid, 1, IQR)
sel.rn=tapply(1:nrow(annots), annots$gene_id, function(x){
  x[which.max(iqrs[x])]
})
gse16873.egid=gse16873.affyid[sel.rn,]
rownames(gse16873.egid)=names(sel.rn)

cn=colnames(gse16873.egid)
hn=grep('HN',cn, ignore.case =T)
dcis=grep('DCIS',cn, ignore.case =T)
data(kegg.gs)
gse16873.kegg.p.affy <- gage(gse16873.egid, gsets = kegg.gs,
                             ref = hn, samp = dcis)
#result should be similar to that of using gse16873


###################################################
### code chunk number 6: pathview.conversion
###################################################
library(pathview)
data(bods)
print(bods)
#simulated human expression data with RefSeq ID
refseq.data <- sim.mol.data(mol.type = "gene", id.type = "REFSEQ",
                            nexp = 2, nmol = 1000)
#construct map between non-Entrez ID and Entrez Gene ID
id.map.refseq <- id2eg(ids = rownames(refseq.data), category =
                         "REFSEQ", org = "Hs")
#Map data onto Entrez Gene IDs, note different sum.method can be used
entrez.data <- mol.sum(mol.data = refseq.data, id.map = id.map.refseq,
                       sum.method = "mean")



###################################################
### code chunk number 1: install (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite(c("gage","gageData"))


cuff.res=read.delim(file="gene_exp.diff", sep="\t")
#notice the column name special character changes. The column used to be
#cuff.res$log2.fold_change. for older versions of Cufflinks.
cuff.fc=cuff.res$log2.fold_change.
gnames=cuff.res$gene
sel=gnames!="-"
gnames=as.character(gnames[sel])
cuff.fc=cuff.fc[sel]
names(cuff.fc)=gnames
gnames.eg=pathview::id2eg(gnames, category ="symbol")
sel2=gnames.eg[,2]>""
cuff.fc=cuff.fc[sel2]
names(cuff.fc)=gnames.eg[sel2,2]
range(cuff.fc)
cuff.fc[cuff.fc>10]=10
cuff.fc[cuff.fc< -10]=-10
exp.fc=cuff.fc
out.suffix="cuff"
###################################################
### code chunk number 7: dataPrep
###################################################
data(gse16873)
cn=colnames(gse16873)
hn=grep('HN',cn, ignore.case =T)
adh=grep('ADH',cn, ignore.case =T)
dcis=grep('DCIS',cn, ignore.case =T)
print(hn)
print(dcis)


###################################################
### code chunk number 8: dataPrep
###################################################
data(kegg.gs)
data(go.gs)
lapply(kegg.gs[1:3],head)
head(rownames(gse16873))


###################################################
### code chunk number 9: gage.1direction
###################################################
gse16873.kegg.p <- gage(gse16873, gsets = kegg.gs,
    ref = hn, samp = dcis)
#go.gs here only the first 1000 entries as a fast example.
gse16873.go.p <- gage(gse16873, gsets = go.gs,
    ref = hn, samp = dcis)
str(gse16873.kegg.p, strict.width='wrap')
head(gse16873.kegg.p$greater[, 1:5], 4)
head(gse16873.kegg.p$less[, 1:5], 4)
head(gse16873.kegg.p$stats[, 1:5], 4)


###################################################
### code chunk number 10: gage.2direction
###################################################
gse16873.kegg.2d.p <- gage(gse16873, gsets = kegg.gs,
    ref = hn, samp = dcis, same.dir = F)
head(gse16873.kegg.2d.p$greater[,1:5], 4)
head(gse16873.kegg.2d.p$stats[,1:5], 4)


###################################################
### code chunk number 11: page
###################################################
gse16873.kegg.page.p <- gage(gse16873, gsets = kegg.gs,
    ref = hn, samp = dcis, saaTest = gs.zTest)
head(gse16873.kegg.page.p$greater[, 1:5], 4)


###################################################
### code chunk number 12: gage.alternative (eval = FALSE)
###################################################
## gse16873.kegg.t.p <- gage(gse16873, gsets = kegg.gs,
##     ref = hn, samp = dcis, use.fold = F)


###################################################
### code chunk number 13: gage.alternative (eval = FALSE)
###################################################
## gse16873.kegg.rk.p <- gage(gse16873, gsets = kegg.gs,
##     ref = hn, samp = dcis, rank.test = T)


###################################################
### code chunk number 14: gage.alternative (eval = FALSE)
###################################################
## gse16873.kegg.ks.p <- gage(gse16873, gsets = kegg.gs,
##     ref = hn, samp = dcis, saaTest = gs.KSTest)


###################################################
### code chunk number 15: gage.alternative (eval = FALSE)
###################################################
## gse16873.kegg.unpaired.p <- gage(gse16873, gsets = kegg.gs,
##     ref = hn, samp = dcis, compare = "unpaired")


###################################################
### code chunk number 16: gage.alternative (eval = FALSE)
###################################################
## gse16873.kegg.gamma.p <- gage(gse16873, gsets = kegg.gs,
##     ref = hn, samp = dcis, use.stouffer=FALSE)


###################################################
### code chunk number 17: symbol.ID
###################################################
head(rownames(gse16873))
gse16873.sym<-gse16873
data(egSymb)
rownames(gse16873.sym)<-eg2sym(rownames(gse16873.sym))
head(rownames(gse16873.sym))


###################################################
### code chunk number 18: symbol.ID
###################################################
kegg.gs.sym<-lapply(kegg.gs, eg2sym)
lapply(kegg.gs.sym[1:3],head)
gse16873.kegg.p2 <- gage(gse16873.sym, gsets = kegg.gs.sym,
    ref = hn, samp = dcis)


###################################################
### code chunk number 19: full.table
###################################################
write.table(gse16873.kegg.2d.p$greater, file = "gse16873.kegg.2d.p.txt",
    sep = "\t")
write.table(rbind(gse16873.kegg.p$greater, gse16873.kegg.p$less),
    file = "gse16873.kegg.p.txt", sep = "\t")


###################################################
### code chunk number 20: significant.genesets
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
### code chunk number 21: nonredundant.genesets
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
### code chunk number 22: sel.gene.expression
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
### code chunk number 23: all.gene.expression
###################################################
for (gs in rownames(gse16873.kegg.p$greater)[1:3]) {
    outname = gsub(" |:|/", "_", substr(gs, 10, 100))
    outname = paste(outname, "all", sep=".")
    geneData(genes = kegg.gs[[gs]], exprs = gse16873, ref = hn,
        samp = dcis, outname = outname, txt = T, heatmap = T,
        Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}


###################################################
### code chunk number 24: pathview
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
### code chunk number 25: pathview.all (eval = FALSE)
###################################################
## sel <- gse16873.kegg.p$greater[, "q.val"] < 0.1 & !is.na(gse16873.kegg.p$greater[,
## "q.val"])
## path.ids <- rownames(gse16873.kegg.p$greater)[sel]
## path.ids2 <- substr(path.ids, 1, 8)
## pv.out.list <- sapply(path.ids2, function(pid) pathview(gene.data = gse16873.d[,
##                         1:2], pathway.id = pid, species = "hsa"))


###################################################
### code chunk number 26: gagePipe
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
### code chunk number 27: gageComp
###################################################
load('gse16873.gage.RData')
gageComp(sampnames, dataname, gsname = "kegg.gs",
    do.plot = TRUE)


###################################################
### code chunk number 28: heter.gage
###################################################
boxplot(data.frame(gse16873))
gse16873.kegg.heter.p <- heter.gage(gse16873, gsets = kegg.gs,
    ref.list = refList, samp.list = sampList)
gse16873.kegg.heter.2d.p <- heter.gage(gse16873, gsets = kegg.gs,
    ref.list = refList, samp.list = sampList, same.dir = F)


