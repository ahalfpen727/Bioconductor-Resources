## source('openGraphSaveGraph.r');

### R code from vignette source 'vignettes/cummeRbund/inst/doc/cummeRbund-manual.Rnw'
## options(echo=TRUE)
## args = commandArgs(trailingOnly = T)
args = commandArgs(TRUE)

for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
}

###################################################
### code chunk number 1: init
###################################################
options(width=65)


###################################################
### code chunk number 2: loadLib
###################################################
library(cummeRbund)


###################################################
### code chunk number 3: read
###################################################
## myDir<-system.file("extdata", package="cummeRbund") #You can leave blank if cwd or replace with your own directory path.
## gtfFile<-system.file("extdata/chr1_snippet.gtf",package="cummeRbund") #path to .gtf file used in cuffdiff analysis.
## cuff <- readCufflinks(dir=myDir,gtfFile=gtfFile,genome="hg19",rebuild=T)


## inDir = args[1]
## outDir = args[2]
cuff <- readCufflinks(dir=inDir)

###################################################
### code chunk number 4: read2 (eval = FALSE)
###################################################
cuff<-readCufflinks()
cuff<-readCufflinks(rebuild=T)
cuff

###################################################
### set working directory to output directory
###################################################
setwd(outDir)

###################################################
### code chunk number 5: read3
###################################################
cuff

###################################################
### code chunk number 30: data_access_0
###################################################
runInfo(cuff)
# this code allows you to see the total map mass per sample lane
# the disparate FPKM masses as a result of different concentrations
# can be grouped by sample and does not represent actual noise
# but it will be interpreted as noise and is there for a confounding factor
# suggest for future analysis consistent concentration throughout except for
# a spike in for normalization purposes.

replicates(cuff)


###################################################
### code chunk number 31: data_access_1
###################################################
gene.features<-annotation(genes(cuff))
head(gene.features)
dim(gene.features)

gene.fpkm<-fpkm(genes(cuff))
head(gene.fpkm)
dim(gene.fpkm)

gene.repFpkm<-repFpkm(genes(cuff))
head(gene.repFpkm)

gene.counts<-count(genes(cuff))
head(gene.counts)

isoform.fpkm<-fpkm(isoforms(cuff))
head(isoform.fpkm)

gene.diff<-diffData(genes(cuff))
head(gene.diff)

sample.names<-samples(genes(cuff))
head(sample.names)

gene.featurenames<-featureNames(genes(cuff))
head(gene.featurenames)

fpkm(genera)

fpkmMatrix(genera)
###################################################
### code chunk number 33: data_access_3
###################################################
gene.matrix<-fpkmMatrix(genes(cuff))
head(gene.matrix)
dim(gene.matrix)

###################################################
### code chunk number 34: data_access_4
###################################################
gene.rep.matrix<-repFpkmMatrix(genes(cuff))
head(gene.rep.matrix)

isoform.rep.matrix<-repFpkmMatrix(isoforms(cuff))
head(isoform.rep.matrix)

tss.rep.matrix<-repFpkmMatrix(TSS(cuff))
head(tss.rep.matrix)

cds.rep.matrix<-repFpkmMatrix(CDS(cuff))
head(cds.rep.matrix)

gene1.rep.matrix<-repFpkmMatrix(genera)
head(gene1.rep.matrix)

isoform1.rep.matrix<-repFpkmMatrix(genera)
head(isoform1.rep.matrix)

tss1.rep.matrix<-repFpkmMatrix(genera)
head(tss1.rep.matrix)

cds1.rep.matrix<-repFpkmMatrix(genera)
head(cds1.rep.matrix)


###################################################
### code chunk number 35: data_access_5
###################################################
gene.count.matrix<-countMatrix(genes(cuff))
head(gene.count.matrix)

gene1.count.matrix<-countMatrix(genera)
head(gene1.count.matrix)
print(gene1.count.matrix)

###################################################
### code chunk number 70: sig_mat_plot_1
###################################################
mySigMat<-sigMatrix(cuff,level='genes',alpha=0.05)

print(mySigMat)


###################################################
### code chunk number 71: get_sig_1
###################################################
mySigGeneIds<-getSig(cuff,alpha=0.05,level='genes')
head(mySigGeneIds)
length(mySigGeneIds)

siggenes<-getSig(cuff,x="LUTS",y="CTRL",alpha=0.05,level='genes')
siggenes
genera<-getGenes(cuff, siggenes)
print(CDS(genera))
###################################################
### code chunk number 73: get_sig_3
###################################################
CXCL12<-getGene(cuff,"CXCL12")
CXCL12

###################################################
### code chunk number 74: get_sig_4
###################################################
mySigTable<-getSigTable(cuff,alpha=0.05,level='genes')
head(mySigTable)
length(mySigTable)
mySigTable
###################################################
### code chunk number 6: add_features
###################################################


fet<-getFeatures(cuff, featureIdList = "gene_id", sampleIdList = NULL)
fet

wh<-addFeatures(cuff,gene.features,featureIdList="locus")
wh
###################################################
### code chunk number 7: global_dispersion
###################################################
disp<-dispersionPlot(genes(cuff))
disp

disp1<-dispersionPlot(isoforms(cuff))
disp1

disp2<-dispersionPlot(CDS(cuff))
disp2
print(disp)


###################################################
### code chunk number 9: SCV_visualization
###################################################
genes.scv<-fpkmSCVPlot(genes(cuff))
genes.scv
isoforms.scv<-fpkmSCVPlot(isoforms(cuff))
isoforms.scv
cvcds<-fpkmSCVPlot(TSS(cuff))
cvcds
cvcds1<-fpkmSCVPlot(CDS(cuff))
cvcds1

###################################################
### code chunk number 10: global_plots_1
###################################################
dens<-csDensity(genes(cuff))
dens
densRep<-csDensity(genes(cuff),replicates=T)
densRep


###################################################
### code chunk number 11: global_plots_dens
###################################################
dens<-csDensity(genes(cuff))
dens
densRep<-csDensity(genes(cuff),replicates=T)
densRep
print(dens)
print(densRep)


###################################################
### code chunk number 13: global_plots_2
###################################################
b<-csBoxplot(genes(cuff))
b
brep<-csBoxplot(genes(cuff),replicates=T)
brep


b1<-csBoxplot(isoforms(cuff))
b1
brep1<-csBoxplot(isoforms(cuff),replicates=T)
brep1

b2<-csBoxplot(CDS(cuff))
b2
brep2<-csBoxplot(CDS(cuff),replicates=T)
brep2
brep3<-csBoxplot(TSS(cuff),replicates=T)
brep3

###################################################
### code chunk number 16: global_plots_3.1
###################################################
s<-csScatterMatrix(genes(cuff))
s


###################################################
### code chunk number 17: global_plots_scatter_1
###################################################
s<-csScatterMatrix(genes(cuff))

print(s)


###################################################
### code chunk number 18: global_plots_3.2
###################################################
s<-csScatter(genes(cuff),"LUTS","CTRL",smooth=T)
s


###################################################
### code chunk number 19: global_plots_scatter_2
###################################################
s<-csScatter(genera,"LUTS","CTRL",smooth=T)
s
print(s)


###################################################
### code chunk number 20: global_plots_6
###################################################
dend<-csDendro(genes(cuff))
dend.rep<-csDendro(genes(cuff),replicates=T)


###################################################
### code chunk number 21: global_plots_dendro
###################################################
dend<-csDendro(genes(cuff))
dend.rep<-csDendro(genes(cuff),replicates=T)
plot(dend)
plot(dend.rep)


###################################################
### code chunk number 23: global_plots_4
###################################################
m<-MAplot(genes(cuff),"LUTS","CTRL")
m
mCount<-MAplot(genes(cuff),"LUTS","CTRL",useCount=T)
mCount



###################################################
### code chunk number 26: global_plots_5_1
###################################################
v<-csVolcanoMatrix(genes(cuff))
v


###################################################
### code chunk number 27: global_plots_volcano_1
###################################################
v<-csVolcanoMatrix(genes(cuff))
v
print(v)


###################################################
### code chunk number 28: global_plots_5_2
###################################################
v<-csVolcano(genes(cuff),"LUTS","CTRL")
v
print(v)


###################################################
### code chunk number 42: geneset_plots_1.5
###################################################
b<-expressionBarplot(genera)
b


###################################################
### code chunk number 43: geneset_plots_barplot
###################################################
b<-expressionBarplot(mySigTable)
b
print(b)


###################################################
### code chunk number 44: geneset_plots_2
###################################################
s<-csScatter(siggenes,"LUTS","CTRL",smooth=T)
s


###################################################
### code chunk number 46: geneset_plots_3
###################################################
v<-csVolcano(genera,"LUTS","CTRL")
v

###################################################
### code chunk number 48: geneset_plots_4
###################################################
ih<-csHeatmap(isoforms(genera),cluster='both',labRow=T)
ih
th<-csHeatmap(TSS(genera),cluster='both',labRow=T)
th

ih2<-csHeatmap(CDS(genera),cluster='both',labRow=T)
ih2
ih2<-csHeatmap(genera,cluster='both',labRow=T)
ih2

###################################################
### code chunk number 51: geneset_plots_5
###################################################
den<-csDendro(genera)
den


###################################################
### code chunk number 52: geneset_plots_dendro
###################################################
den1<-csDendro(genera, replicates=T)
den1


###################################################
### code chunk number 53: gene_level_1
###################################################
myGeneId<-"CXCL12"
myGene<-getGene(cuff,myGeneId)
myGene
myGeneId
fpkm(myGene)
fpkm(isoforms(myGene))


myGeneId1<-"CXCR4"
myGene1<-getGene(cuff,myGeneId1)
myGene1
myGeneId1
goo<-fpkm(myGene1)
woo<-fpkm(isoforms(myGene1))
goo
woo
dude<-repFpkmMatrix(myGene)
dude

dude1<-repFpkmMatrix(myGene1)
dude1

dud<-repFpkmMatrix(TSS(myGene1))
dud

duder<-repFpkmMatrix(CDS(myGene1))
duder

duders<-repFpkmMatrix(isoforms(myGene1))
duders

dde1<-repFpkmMatrix(myGene)
dde1

de1<-repFpkmMatrix(myGene1)
de1

dd<-repFpkmMatrix(TSS(myGene))
dd

dder<-repFpkmMatrix(CDS(myGene))
dder

dders<-repFpkmMatrix(isoforms(myGene))
dders
###################################################
### code chunk number 54: gene_plots_1
###################################################
gl<-expressionPlot(myGene)
gl

gl.rep<-expressionPlot(myGene,replicates=TRUE)
gl.rep

gl.iso.rep<-expressionPlot(isoforms(myGene),replicates=T)
gl.iso.rep

gl.cds.rep<-expressionPlot(myGene,level="isoforms",replicates=T)
gl.cds.rep

z<-expressionPlot(myGene)
z

z1<-expressionPlot(isoforms(myGene))
z1

z2<-expressionPlot(CDS(myGene))
z2

z3<-expressionPlot(TSS(myGene))
z3

g2<-expressionPlot(myGene1)
g2

g2.rep<-expressionPlot(myGene1,replicates=TRUE)
g2.rep

g2.iso.rep<-expressionPlot(isoforms(myGene1),replicates=T)
g2.iso.rep

g2.cds.rep<-expressionPlot(CDS(myGene1),replicates=T)
g2.cds.rep

g3.cds.rep<-expressionPlot(TSS(myGene1),replicates=T)
g3.cds.rep

###################################################
### code chunk number 55: gene_plots_line
###################################################
gl<-expressionPlot(mySigGenes)
gl

gl.rep<-expressionPlot(myGene,replicates=TRUE)
gl.rep

gl.iso.rep<-expressionPlot(isoforms(myGene),replicates=T)
gl.iso.rep

gl.cds.rep<-expressionPlot(CDS(myGene),replicates=T)
gl.cds.rep
print(gl)

###################################################
### code chunk number 62: gene_plots_3
###################################################
igb<-expressionBarplot(isoforms(myGene),replicates=T)
igb

igb1<-expressionBarplot((myGene),replicates=T)
igb1

igb2<-expressionBarplot(CDS(myGene),replicates=T)
igb2

igb3<-expressionBarplot(TSS(myGene),replicates=T)
igb3
###################################################
### code chunk number 63: gene_plots_bar_isoforms
###################################################
igb<-expressionBarplot(mySigGenes,replicates=T)
igb
print(igb)
igb<-expressionBarplot(isoforms(myGene1),replicates=T)
igb

igb1<-expressionBarplot((myGene1),replicates=T)
igb1

igb2<-expressionBarplot(CDS(myGene1),replicates=T)
igb2

igb3<-expressionBarplot(TSS(myGene1),replicates=T)
igb3

###################################################
### code chunk number 64: gene_plots_4
###################################################
gp<-csPie(myGene,level="isoforms")
gp
gp1<-csPie(myGene,fpkm,level="TSS")
gp1
gp2<-csPie(myGene,fpkm,level="CDS",replicates=T)
gp2

zp<-csPie(myGene1,level="isoforms")
zp
zp1<-csPie(myGene1,fpkm,level="TSS")
zp1
zp2<-csPie(myGene1,fpkm,level="CDS",replicates=T)
zp2
####################################
###################################################
### code chunk number 75: dist_heat_1
###################################################
myDistHeat<-csDistHeat(mySigGenes)



###################################################
### code chunk number 76: dist_heat_plot_1
###################################################
myDistHeat<-csDistHeat(genes(cuff))

print(myDistHeat)


###################################################
### code chunk number 77: dist_heat_2
###################################################
myRepDistHeat<-csDistHeat(genes(cuff),replicates=T)



###################################################
### code chunk number 78: dist_heat_plot_2
###################################################
myRepDistHeat<-csDistHeat(genes(cuff),replicates=T)

print(myRepDistHeat)


###################################################
### code chunk number 79: dim_reduction_1
###################################################
genes.PCA<-PCAplot(genes(cuff),"PC1","PC2")
genes.MDS<-MDSplot(genes(cuff))

genes.PCA.rep<-PCAplot(genes(cuff),"PC1","PC2",replicates=T)
genes.MDS.rep<-MDSplot(genes(cuff),replicates=T)


###################################################
### code chunk number 80: gene_PCA
###################################################
genes.PCA<-PCAplot(genes(cuff),"PC1","PC2")
genes.MDS<-MDSplot(genes(cuff))
print(genes.PCA)
genes.PCA.rep<-PCAplot(genes(cuff),"PC1","PC2",replicates=T)
genes.MDS.rep<-MDSplot(genes(cuff),replicates=T)
print(genes.PCA.rep)
print(genes.MDS.rep)


###################################################
### code chunk number 81: gene_MDS
###################################################

print(genes.MDS)


###################################################
### code chunk number 82: gene_PCA_rep
###################################################

print(genes.PCA.rep)


###################################################
### code chunk number 83: gene_MDS_rep
###################################################

print(genes.MDS.rep)


###################################################
### code chunk number 84: geneset_cluster_1
###################################################
ic<-csCluster(mySigGenes,k=10)
head(ic$cluster)
icp<-csClusterPlot(ic)
icp


###################################################
### code chunk number 85: geneset_plots_cluster
###################################################
print(icp)


###################################################
### code chunk number 86: specificity_1
###################################################
myGenes.spec<-csSpecificity(mySigGenes)
head(myGenes.spec)


###################################################
### code chunk number 87: similar_1
###################################################
mySimilar<-findSimilar(cuff,"CXCL12",n=20)
mySimilar.expression<-expressionPlot(mySimilar,logMode=T,showErrorbars=F)
print(mySimilar.expression)

mySimilar1<-findSimilar(cuff,"CXCR4",n=10)
mySimilar1.expression<-expressionPlot(mySimilar1,logMode=T,showErrorbars=F)
print(mySimilar1.expression)

###################################################
### code chunk number 89: similar_2
###################################################
myProfile<-c(500,500000)
mySimilar2<-findSimilar(cuff,myProfile,n=50)
mySimilar2.expression<-expressionPlot(mySimilar2,logMode=T,showErrorbars=F)
mySimilar2.expression

###################################################
### code chunk number 90: similar_plots_2
###################################################
myProfile<-c(500,0,400)
mySimilar2<-findSimilar(cuff,myProfile,n=100)
mySimilar2.expression<-expressionPlot(mySimilar2,logMode=T,showErrorbars=F)
print(mySimilar2.expression)


###################################################
### code chunk number 91: close_connection
###################################################
end<-dbDisconnect(cuff@DB)


###################################################
### code chunk number 92: session
###################################################
sessionInfo()


