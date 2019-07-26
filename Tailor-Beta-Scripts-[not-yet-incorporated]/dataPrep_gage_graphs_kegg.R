#dataPrep.Rnw
### R code from vignette source 'dataPrep.Rnw'

library(gage); filename=system.file("extdata/gse16873.demo", package = "gage")
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
library(gage)
filename=system.file("extdata/gse16873.demo", package = "gage")
demo.data=readExpData(filename, row.names=1)
head(demo.data);str(demo.data)
#convert the data.frame into a matrix as to speed up the computing
demo.data=as.matrix(demo.data)
head(demo.data);str(demo.data)

#an example GMT gene set data derived from MSigDB data

filename=system.file("extdata/c2.demo.gmt", package = "gage")
demo.gs=readList(filename)
demo.gs

#to use these gene sets with gse16873, need to convert the gene symbols
#to Entrez IDs first
#genes in gse16873 was label by Entrez IDs
data(gse16873)
head(rownames(gse16873))
#may convert the gene IDs to official symbols
gse16873.sym<-gse16873
data(egSymb)
rownames(gse16873.sym)<-eg2sym(rownames(gse16873.sym))
head(rownames(gse16873.sym))

#convert kegg.gs correspondingly
data(kegg.gs)
kegg.gs.sym<-lapply(kegg.gs, eg2sym)
lapply(kegg.gs.sym[1:3],head)
#GAGE analysis with the converted data
cn=colnames(gse16873)
hn=grep('HN',cn, ignore.case =TRUE)
dcis=grep('DCIS',cn, ignore.case =TRUE)
gse16873.kegg.p2 <- gage(gse16873.sym, gsets = kegg.gs.sym,
    ref = hn, samp = dcis)
data(eg)
demo.gs.sym<-lapply(demo.gs, sym2eg)
(egSymb)

###################################################
### code chunk number 5: gse16873.affyid
###################################################
library(gageData)
data(gse16873.affyid)
affyid=rownames(gse16873.affyid)
# gse16873 -->	GSE16873: a breast cancer microarray dataset
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

library( gage )
kg.hsa <- kegg.gsets( "hsa" )
kegg.gs2 <- kg.hsa$kg.sets[ kg.hsa$sigmet.idx ]
res <- gage( E, gsets= kegg.gs2, 
       ref= which( group == "Control" ), 
       samp= which( group == "Treatment" ),
       compare= "unpaired", same.dir= FALSE )
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
library(gage)
###################################################
filename=system.file("extdata/gse16873.demo", package = "gage")
demo.data=readExpData(filename, row.names=1)
head(demo.data)
str(demo.data)
#convert the data.frame into a matrix (speeds up computing)
demo.data=as.matrix(demo.data)
str(demo.data)

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
### gse16873.affyid ex
###################################################
#source("http://bioconductor.org/biocLite.R")
#biocLite("gageData")
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
## Data Prep for graphs
###################################################
source("http://bioconductor.org/biocLite.R")
biocLite("RbcBook1")
#biocLite("SpikeInSubset")
library("RbcBook1")
library("Biobase")
stopifnot(package.version("RbcBook1") >= package_version("0.3.0"))
###################################################
### chunk number 2: affydist1
###################################################
library("affy")
library("SpikeInSubset")
library("RColorBrewer")
data(spikein95)
hist(spikein95,lwd=1,xlab=expression(log[2](intensity)),col=brewer.pal(8,"Dark2"))

boxplot(spikein95,names=1:6,col="#B2DF8A",
 ylab=expression(log[2](intensity)),xlab="array")
## library("affy")
## library("SpikeInSubset")
## data("spikein95")
## hist(spikein95)
## boxplot(spikein95)


###################################################
### chunk number 5: biasCalc
###################################################
data(spikein133)
Index <- which(probeNames(spikein133)%in%colnames(pData(spikein133)))
pms <- pm(spikein133)[Index,]
genenames <- probeNames(spikein133)[Index]
nominal <- t(sapply(probeNames(spikein133)[Index],function(i) pData(spikein133)[,i]))
x <- as.vector(log2(nominal))
y <- as.vector(log2(pms))
avg <- tapply(y,x,mean,na.rm=TRUE)


###################################################
### chunk number 6: bias
###################################################
plot(jitter(x),y,las=1,pch=".",
     ylab="Log (Base 2) PM Intensity",
     xlab="Nominal Log (Base 2) Concentration")
lines(as.numeric(names(avg)),avg,lwd=3,col="red")

plot(x,y,las=1,pch=".",
     ylab="Log (Base 2) PM Intensity",
     xlab="Nominal Log (Base 2) Concentration")
lines(as.numeric(names(avg)),avg,lwd=3,col="red")

###################################################
### chunk number 7: lognormalerrors
###################################################
qqnorm(y[x==0],pch=".")
qqline(y[x==0],col="red")





