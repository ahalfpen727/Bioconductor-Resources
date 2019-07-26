### R code from vignette source 'dataPrep.Rnw'

###################################################
### code chunk number 1: start
###################################################
options(width=80)


###################################################
### code chunk number 2: start
###################################################
library(gage)


###################################################
### code chunk number 3: demo.data
###################################################
filename=system.file("extdata/gse16873.demo", package = "gage")
demo.data=readExpData(filename, row.names=1)
#check the data
head(demo.data)
str(demo.data)
#convert the data.frame into a matrix as to speed up the computing
demo.data=as.matrix(demo.data)
str(demo.data)


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


