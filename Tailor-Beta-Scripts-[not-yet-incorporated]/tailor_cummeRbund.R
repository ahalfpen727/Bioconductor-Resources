## source('openGraphSaveGraph.r');
options(error=traceback) # causes a traceback to appear if there is an error, and the traceback has the line number
# prefixed by #

### R code from vignette source 'vignettes/cummeRbund/inst/doc/cummeRbund-manual.Rnw'
## options(echo=TRUE)
## args = commandArgs(trailingOnly = T)
args = commandArgs(TRUE)

for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
}



###################################################
### code chunk number 1: init, 2: loadLib, 3: Read
###################################################
options(width=65)

library(cummeRbund)

##################################################
#### Get environment variable
################################################
refgtf = Sys.getenv("CUMMERBUND_INPUT_GTF")
genome_path = Sys.getenv("GENOME_PATH")
gen_v = Sys.getenv("R_GENOME")

## myDir<-system.file("extdata", package="cummeRbund") #leave blank if cwd or replace with own directory path.
## gtfFile<-system.file("extdata/chr1_snippet.gtf",package="cummeRbund") #path to .gtf file used

cuff <- readCufflinks(dir=inDir,gtfFile=refgtf,genome=genome_path,rebuild=T)

cuff
###################################################
### 5: set working directory to output directory
###################################################
setwd(outDir)

mySigGeneIds<-getSig(cuff,x=over,y=under,alpha=.05,level='genes')
head(mySigGeneIds)
length(mySigGeneIds)

sigGenes<-getGenes(cuff, mySigGeneIds)
head(sigGenes)
#length(sigGenes)


#gene.diff<-diffData(genes(cuff))
#head(gene.diff)
#length(gene.diff)

mySigGeneIds<-getSig(cuff,x=over,y=under,alpha=.05,level='genes')
#mySigGeneIds
#mySigGeneIds2<-getSig(cuff,x=over,y=under,alpha=.05,level='isoforms')
#mySigGeneIds2
mySigGeneIds3<-getSig(cuff,x=over,y=under,alpha=.05,level='CDS')
#mySigGeneIds3
#mySigGeneIds4<-getSig(cuff,x=over,y=under,alpha=.05,level='TSS')
#mySigGeneIds4
#mySigGeneIds5<-getSig(cuff,x=over,y=under,alpha=.05,level='relCDS')
#mySigGeneIds5
#mySigGeneIds1<-getSig(cuff,x=over,y=under,alpha=.05,level='splicing')
#mySigGeneIds1

runInfo(cuff)
replicates(cuff)


#gene.fpkm<-fpkm(genes(cuff))
#head(gene.fpkm)

#gene.repFpkm<-repFpkm(genes(cuff))
#head(gene.repFpkm)

#gene.counts<-count(genes(cuff))
#head(gene.counts)

#isoform.fpkm<-fpkm(isoforms(cuff))
#head(isoform.fpkm)

#head(fpkmMatrix(sigGenes))
sample.names<-samples(genes(cuff))
#head(sample.names)
gene.featurenames<-featureNames(genes(cuff))
#head(gene.featurenames)


###################################################
### code chunk number7: fpkmMatrix
###################################################
#gene.matrix<-fpkmMatrix(genes(cuff))
#head(gene.matrix)
#iso.matrix<-fpkmMatrix(isoforms(cuff))
#head(iso.matrix)
#cds.matrix<-fpkmMatrix(CDS(cuff))
#head(cds.matrix)
#tss.matrix<-fpkmMatrix(TSS(cuff))
#head(tss.matrix)

###################################################
### code chunk number 8: repFPKMMatrix
###################################################
#gene.rep.matrix<-repFpkmMatrix(genes(cuff))
#head(gene.rep.matrix)
#isoform.rep.matrix<-repFpkmMatrix(isoforms(cuff))
#head(isoform.rep.matrix)
#tss.rep.matrix<-repFpkmMatrix(TSS(cuff))
#head(tss.rep.matrix)
#cds.rep.matrix<-repFpkmMatrix(CDS(cuff))
#head(cds.rep.matrix)
#gene1.rep.matrix<-repFpkmMatrix(sigGenes)
#head(gene1.rep.matrix)
#isoform1.rep.matrix<-repFpkmMatrix(CDS(sigGenes))
#head(isoform1.rep.matrix)
#tss1.rep.matrix<-repFpkmMatrix(isoforms(sigGenes))
#head(tss1.rep.matrix)
#cds1.rep.matrix<-repFpkmMatrix(TSS(sigGenes))
#head(cds1.rep.matrix)


###################################################
### code chunk number 9: countMatrix
###################################################
#gene.count.matrix<-countMatrix(genes(cuff))
#head(gene.count.matrix)
#gene1.count.matrix<-countMatrix(sigGenes)
#head(gene1.count.matrix)
#tss.count.matrix<-countMatrix(TSS(sigGenes))
#head(tss.count.matrix)
#iso.count.matrix<-countMatrix(isoforms(sigGenes))
#head(iso.count.matrix)
#cds.count.matrix<-countMatrix(CDS(sigGenes))
#head(cds.count.matrix)

###################################################
### code chunk number 36: create_geneset_1
###################################################
geneListString = Sys.getenv("CUMMERBUND_GENE_LIST")
myGenelist = unlist(strsplit(geneListString, " "))
myGenelist
myGeneIds=myGenelist

myGenes<-getGenes(cuff,myGeneIds)
myGenes

myGenes[is.na(myGenes)] <- 0

###################################################
### code chunk number 37: create_geneset_2
###################################################
#FPKM values for genes in gene set
#head(fpkm(myGenes))

#Isoform-level FPKMs for gene set
#head(fpkm(isoforms(myGenes)))

#Replicate FPKMs for TSS groups within gene set
#head(repFpkm(TSS(myGenes)))

##################################################
### code chunk number 7: global_dispersion
###################################################
disp<-dispersionPlot(genes(cuff))
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
dens2<-csDensity(genes(cuff),replicates=T)
print(dens2)
dens3<-csDensity(isoforms(cuff), replicates=T)
print(dens3)
dens4<-csDensity(CDS(cuff), replicates=T)
print(dens4)
dens6<-csDensity(TSS(cuff), replicates=T)
print(dens6)
###################################################
### code chunk number 13: global_plots_2
###################################################

b<-csBoxplot(genes(cuff)) + ggtitle("genes")
b
brep<-csBoxplot(genes(cuff),replicates=T) + ggtitle("genes")
brep
brep1<-csBoxplot(isoforms(cuff),replicates=T) + ggtitle("isoforms")
brep1
brep2<-csBoxplot(CDS(cuff),replicates=T) + ggtitle("CDS")
brep2
brep3<-csBoxplot(TSS(cuff),replicates=T) + ggtitle("TSS")
brep3

###################################################
### code chunk number 17: global_plots_scatter_1
###################################################
s<-csScatterMatrix(genes(cuff)) + ggtitle("genes")
print(s)

s<-csScatter(genes(cuff),over,under,smooth=T) + ggtitle("genes")
s
s1<-csScatter(isoforms(cuff), over, under,smooth=T) + ggtitle("isoforms")
s1

S1<-csScatter(TSS(cuff), over, under,smooth=T) + ggtitle("TSS")
S1
s2<-csScatter(CDS(cuff), over, under,smooth=T) + ggtitle("CDS")
s2

###################################################
### code chunk number 19: global_plots_scatter_2
###################################################
s22<-csScatter(sigGenes,over,under,smooth=T) + ggtitle("sigGenes")
s22
#s2212<-csScatter(isoforms(sigGenes),over, under,smooth=T) + ggtitle("isoforms")
#s2212
#s2213<-csScatter(TSS(sigGenes),over, under,smooth=F) + ggtitle("TSS")
#s2213


###################################################
### code chunk number 20: global_plots_6
###################################################

###################################################
### code chunk number 21: global_plots_dendro
###################################################
#dend<-csDendro(genes(cuff),replicates=F)
## plot(dend) + ggtitle("genes, replicates=F")

dend.rep<-csDendro(genes(cuff),replicates=T)
## plot(dend.rep) + ggtitle("genes, replicates=T")

dend.rep<-csDendro(TSS(cuff),replicates=T)
## plot(dend.rep) + ggtitle("TSS")

dend.rep<-csDendro(isoforms(cuff),replicates=T)
## plot(dend.rep) + ggtitle("isoforms")

dend.rep<-csDendro(CDS(cuff),replicates=T)
## plot(dend.rep) + ggtitle("CDS")


###################################################
### code chunk number 23: global_plots_4
###################################################
m<-MAplot(genes(cuff),over,under) + ggtitle("genes")
m

mCount1<-MAplot(genes(cuff),over, under,useCount=T) + ggtitle("genes")
mCount1

mCount2<-MAplot(isoforms(cuff),over,under,useCount=T) + ggtitle("isoforms")
mCount2

m9<-MAplot(CDS(cuff),over,under,useCount=T) + ggtitle("CDS")
m9

mCount8<-MAplot(TSS(cuff),over, under,useCount=T) + ggtitle("TSS")
mCount8

###################################################
### code chunk number 26: global_plots_5_1
###################################################
 v12<-csVolcanoMatrix(genes(cuff)) + ggtitle("genes")
 v12
 v11<-csVolcanoMatrix(isoforms(cuff)) + ggtitle("isoforms")
 v11
 v10<-csVolcanoMatrix(CDS(cuff)) + ggtitle("CDS")
 v10
 v9<-csVolcanoMatrix(TSS(cuff)) + ggtitle("TSS")
 v9


###################################################
### code chunk number 29: global_plots_volcano_2
###################################################
vv<-csVolcano(genes(cuff),over, under) + ggtitle("genes")
vv
vy<-csVolcano(TSS(cuff),over, under) + ggtitle("TSS")
vy
vu<-csVolcano(CDS(cuff),over, under) + ggtitle("CDS")
vu
vu<-csVolcano(isoforms(cuff),over, under) + ggtitle("isoforms")
vu

###################################################
### code chunk number 38: create_geneset_3
###################################################
# this was modified to try and fix the loop
# myGenesetis<-getGenes(cuff,myGeneIds,sampleIdList=c(over,under))
# myGenesetis

###################################################
### code chunk number 40: geneset_plots_heatmap
###################################################
## h<-csHeatmap(myGenes,cluster='both')
## h
## h.rep<-csHeatmap(myGenes,cluster='both',replicates=T)
## print(h.rep)
## print(h)


##################################################
### code chunk number 42: geneset_plots_1.5
###################################################
b<-expressionBarplot(myGenes) + ggtitle("myGenes")
b
#b11<-expressionBarplot(CDS(sigGenes)) + ggtitle("sigGenes")
#b11
#b12<-expressionBarplot(isoforms(sigGenes)) + ggtitle("sigGenes")
#b12
#b13<-expressionBarplot(TSS(sigGenes)) + ggtitle("sigGenes")
#b13


###################################################
### code chunk number 44: geneset_plots_2
###################################################
sc<-csScatter(myGenes,over,under,smooth=T) + ggtitle("myGenes")
sc
print(sc)
sc1<-csScatter(sigGenes,over,under,smooth=T) + ggtitle("sigGenes")
sc1
print(sc1)


###################################################
### code chunk number 46: geneset_plots_3
###################################################
vs<-csVolcano(myGenes,under,over) + ggtitle("myGenes")
vs
print(vs)
vsa<-csVolcano(sigGenes,over,under) + ggtitle("sigGenes")
vsa
print(vsa)
vol<-csVolcano(TSS(sigGenes),over,under) + ggtitle("TSS(sigGenes)")
vol
print(vol)
volvo<-csVolcano(CDS(sigGenes),over,under) + ggtitle("CDS(sigGenes)")
volvo
print(volvo)
vol1<-csVolcano(isoforms(sigGenes),over,under) + ggtitle("isoforms(sigGenes)")
vol1
print(vol1)


###################################################
### code chunk number 48: geneset_plots_4
###################################################
h11<-csHeatmap(myGenes,cluster='sigGenes', labRow=F) + ggtitle("myGenes")
h11
h.rep123<-csHeatmap(myGenes,cluster='both',replicates=T, labRow=F) + ggtitle("myGenes")
h.rep123
h22<-csHeatmap(myGenes,cluster='both',labRow=F) + ggtitle("myGenes")
h22

h.cs.rep<-csHeatmap(sigGenes, cluster="row",labRow=F, replicates=F) + ggtitle("sigGenes")
h.cs.rep
#h.cs.rep1<-csHeatmap(isoforms(sigGenes), cluster="row",labRow=F, replicates=F) + ggtitle("isoforms(sigGenes)")
#h.cs.rep1
#h.cs.rep2<-csHeatmap(TSS(sigGenes), cluster="row",labRow=F, replicates=F) + ggtitle("TSS(sigGenes)")
#h.cs.rep2
#h.cs.rep3<-csHeatmap(CDS(sigGenes), cluster="row",labRow=F, replicates=F) + ggtitle("CDS(sigGenes)")
#h.cs.rep3

h.rep<-csHeatmap(sigGenes, cluster="row", replicates=T, labRow=F) + ggtitle("sigGenes")
h.rep

h.rep<-csHeatmap(sigGenes, cluster="column", replicates=T, labRow=F) + ggtitle("sigGenes")
h.rep
h.cs<-csHeatmap(sigGenes, cluster=over, replicates=T, labRow=F) + ggtitle("sigGenes")
h.cs
#h.cs1<-csHeatmap(isoforms(sigGenes), cluster=over, labRow=F, replicates=T) + ggtitle("isoforms(sigGenes)")
#h.cs1
#h.cs2<-csHeatmap(TSS(sigGenes), cluster=over, replicates=T) + ggtitle("TSS(sigGenes)")
#h.cs2
#h.cs3<-csHeatmap(CDS(sigGenes), cluster=over,labRow=F replicates=T) + ggtitle("CDS(sigGenes)")
#h.cs3

## ih<-csHeatmap(isoforms(myGenes),cluster='both',labRow=F)
## ih
## th<-csHeatmap(TSS(myGenes),cluster='both',labRow=F)
## th

## ih2<-csHeatmap(CDS(myGenes),cluster='both',labRow=F)
## ih2

ih4<-csHeatmap(sigGenes,labRow=F) + ggtitle("sigGenes")
ih4

#ih5<-csHeatmap(isoforms(sigGenes),labRow=F) + ggtitle("isoforms(sigGenes)")
#ih5

#ih6<-csHeatmap(TSS(sigGenes),labRow=F) + ggtitle("TSS(sigGenes)")
#ih6
#ih3<-csHeatmap(CDS(sigGenes),labRow=F) + ggtitle("CDS(sigGenes)")
#ih3

###################################################
### code chunk number 49: geneset_plots_isoform_heatmap
###################################################
## ih7<-csHeatmap(isoforms(myGenes),cluster='both',labRow=T)
## ih7
## th8<-csHeatmap(TSS(myGenes),cluster='both',labRow=T)
## th8
## ih2<-csHeatmap(CDS(myGenes),cluster='both',labRow=T)
## ih2

ih9<-csHeatmap(myGenes,cluster='both',labRow=F) + ggtitle("myGenes")
ih9

## ih<-csHeatmap(isoforms(sigGenes),cluster='both',labRow=T)
## ih
## th<-csHeatmap(TSS(sigGenes),cluster='both',labRow=T)
## th
## ih2<-csHeatmap(CDS(sigGenes),cluster='both',labRow=T)
## ih2

ih15<-csHeatmap(sigGenes,cluster='both',labRow=F) + ggtitle("sigGenes")
ih15


###################################################
### code chunk number 36: create_geneset_1
###################################################

geneListString = Sys.getenv("CUMMERBUND_GENE_LIST")
myGeneIds = unlist(strsplit(geneListString, " "))
myGeneIds
myGenes<-getGenes(cuff,myGeneIds)
myGenes
myGenes[is.na(myGenes)]<-c(0)

###################################################
### code chunk number 37: create_geneset_2
###################################################
#FPKM values for genes in gene set
head(fpkm(myGenes))
#Isoform-level FPKMs for gene set
#head(fpkm(isoforms(myGenes)))
#Replicate FPKMs for TSS groups within gene set
#head(repFpkm(TSS(myGenes)))

##############################################
### Gene Loop ################################
##############################################

# this was modified for the loop
 genesArray = c(myGeneIds)
 features(myGeneIds)
for (i in 1:length(genesArray)) {

    myGeneId = genesArray[i]
    myGeneId
    features(myGeneId)
    print(myGeneId)
###################################################
### code chunk number 53: gene_level_1
###################################################

    myGene<-getGene(cuff,myGeneId)
    myGene
    features(myGene)

###################################################
### code chunk number 53: gene_level_1
###################################################

    repFpkm<-repFpkmMatrix(myGene)
    head(repFpkm)
    repFpkm<-repFpkmMatrix(TSS(myGene))
 #   head(repFpkm)
    repFpkm<-repFpkmMatrix(CDS(myGene))
 #   head(repFpkm)
    repFpkm<-repFpkmMatrix(isoforms(myGene))
 #   head(repFpkm)

###################################################
### code chunk number 54:
###################################################
    gl.rep<-expressionPlot(myGene,replicates=TRUE)
   print(gl.rep)

    gl.cds.rep<-expressionPlot(TSS(myGene),replicates=TRUE)
   print(gl.cds.rep)

    gl.cds.rep<-expressionPlot(CDS(myGene),replicates=TRUE)
   print(gl.cds.rep)

    gl.iso.rep<-expressionPlot(isoforms(myGene),replicates=TRUE)
   print(gl.iso.rep)


###################################################
### code chunk number 55: gene_plots_line
###################################################
    gl<-expressionPlot(myGene, replicates=FALSE)
   print(gl)

    z3<-expressionPlot(TSS(myGene), replicates=FALSE)
    z3

    z2<-expressionPlot(CDS(myGene), replicates=FALSE)
    z2

    z1<-expressionPlot(isoforms(myGene), replicates=FALSE)
    z1


###################################################
### code chunk number 59: gene_plots_2
###################################################
    igb1<-expressionBarplot(myGene,replicates=TRUE)
    igb1

    igb3<-expressionBarplot(TSS(myGene),replicates=TRUE)
    igb3

    igb2<-expressionBarplot(CDS(myGene),replicates=TRUE)
    igb2

    igb<-expressionBarplot(isoforms(myGene),replicates=TRUE)
    igb

###################################################
### code chunk number 60: gene_plots_bar
###################################################

    igb1<-expressionBarplot(myGene,replicates=FALSE)
    igb1

    igb3<-expressionBarplot(TSS(myGene),replicates=FALSE)
    igb3

    igb2<-expressionBarplot(CDS(myGene),replicates=FALSE)
    igb2

###################################################
### code chunk number 64: gene_plots__pie
###################################################
    ## gp1<-csPie(myGene,fpkm,level="genes",replicates=TRUE)
    ## gp1

    gp1<-csPie(myGene,fpkm,level="TSS",replicates=TRUE)
   print(gp1)

    gp2<-csPie(myGene,fpkm,level="CDS",replicates=TRUE)
   print(gp2)

    gp3<-csPie(myGene,level="isoforms",replicates=TRUE)
   print(gp3)

  
###################################################
### code chunk number 65: features_1
###################################################
    ## head(features(myGene))

###################################################
### code chunk number 66: features_2
###################################################
     genetrack<-makeGeneRegionTrack(myGene)
     plotTracks(genetrack)



###################################################
### code chunk number 67: features_3
###################################################
 trackList<-list()
 myStart<-min(features(myGene)$start)
 myEnd<-max(features(myGene)$end)
 myChr<-unique(features(myGene)$seqnames)

 print(gen_v)

 ideoTrack <- IdeogramTrack(genome = gen_v, chromosome = myChr)
 trackList<-c(trackList,ideoTrack)

 axtrack<-GenomeAxisTrack()
 trackList<-c(trackList,axtrack)

 genetrack<-makeGeneRegionTrack(myGene)
 genetrack

 trackList<-c(trackList,genetrack)

 biomTrack<-BiomartGeneRegionTrack(genome=gen_v,chromosome=as.character(myChr),
                                   start=myStart,end=myEnd,name="ENSEMBL",showId=T)

 trackList<-c(trackList,biomTrack)


 conservation <- UcscTrack(genome = gen_v, chromosome = myChr,
                           track = "Conservation", table = "phyloP100wayAll",
                           from = myStart-2000, to = myEnd+2000, trackType = "DataTrack",
                           start = "start", end = "end", data = "score",
                           type = "hist", window = "auto", col.histogram = "darkblue",
                           fill.histogram = "darkblue", ylim = c(-3.7, 4),
                           name = "Conservation")

 trackList<-c(trackList,conservation)

 plotTracks(trackList,from=myStart-2000,to=myEnd+2000)

}


###################################################
### code chunk number 69: sig_mat_plot_1
###################################################
mySigMat<-sigMatrix(cuff,level='genes',alpha=0.05)
print(mySigMat)
sigIsoformIds1<-sigMatrix(cuff, level='isoforms',alpha=0.05)
print(sigIsoformIds1)
sigIsoformIds2<-sigMatrix(cuff, level='TSS',alpha=0.05)
print(sigIsoformIds2)
sigIsoformIds3<-sigMatrix(cuff,level='CDS', alpha=0.05)
print(sigIsoformIds3)

###################################################
### code chunk number 73: get_sig_3
###################################################
mySigGenes<-getGenes(cuff,mySigGeneIds)
mySigGenes



###################################################
### code chunk number 74: get_sig_4
###################################################
mySigTable<-getSigTable(cuff,alpha=0.05,level='genes')
head(mySigTable,20)


###################################################
### code chunk number 75: dist_heat_1
###################################################
myDistHeata<-csDistHeat(genes(cuff)) + ggtitle("genes")
myDistHeata
myDistHeatb<-csDistHeat(isoforms(cuff)) + ggtitle("isoforms")
myDistHeatb
myDistHeatc<-csDistHeat(TSS(cuff)) + ggtitle("TSS")
myDistHeatc
myDistHeatd<-csDistHeat(CDS(cuff)) + ggtitle("CDS")
myDistHeatd

#myDistHeate<-csDistHeat(TSS(sigGenes)) + ggtitle("TSS(sigGenes)")
#myDistHeate
#myDistHeatf<-csDistHeat(CDS(sigGenes)) + ggtitle("CDS(sigGenes)")
#myDistHeatf

#myDistHeatg<-csDistHeat(isoforms(sigGenes)) + ggtitle("isoforms(sigGenes)")
#myDistHeatg

###################################################
### code chunk number 77: dist_heat_2
###################################################
myRepDistHeat1<-csDistHeat(genes(cuff),replicates=TRUE) + ggtitle("genes, replicates=T")
myRepDistHeat1
myRepDistHeat2<-csDistHeat(genes(cuff),replicates=FALSE) + ggtitle("genes, replicates=F")
myRepDistHeat2
myRepDistHeat3<-csDistHeat(isoforms(cuff),replicates=FALSE) + ggtitle("genes, replicates=F")
myRepDistHeat3
myRepDistHeat4<-csDistHeat(TSS(cuff),replicates=FALSE) + ggtitle("TSS, replicates=F")
myRepDistHeat4
myRepDistHeat5<-csDistHeat(CDS(cuff),replicates=FALSE) + ggtitle("CDS, replicates=F")
myRepDistHeat5

#myRepDistHeat6<-csDistHeat(genes(sigGenes), replicates=TRUE) + ggtitle("genes(sigGenes), replicates=T")
#myRepDistHeat6
#myRepDistHeat7<-csDistHeat(isoforms(sigGenes),replicates=TRUE) + ggtitle("isoforms(sigGenes), replicates=T")
#myRepDistHeat7
#myRepDistHeat8<-csDistHeat(TSS(sigGenes),replicates=TRUE) + ggtitle("TSS(sigGenes), replicates=T")
#myRepDistHeat8
#myRepDistHeat9<-csDistHeat(CDS(sigGenes),replicates=TRUE) + ggtitle("CDS(sigGenes), replicates=T")
#myRepDistHeat9

###################################################
### code chunk number 79: dim_reduction_1
###################################################
genes.PCA<-PCAplot(genes(cuff),"PC1","PC2", replicates=FALSE) + ggtitle("genes, replicates=F")
genes.MDS<-MDSplot(genes(cuff), replicates=FALSE, k=2) + ggtitle("genes, replicates=F")

print(genes.PCA)
print(genes.MDS)

genes.PCA.rep<-PCAplot(genes(cuff),"PC1","PC2", replicates=TRUE) + ggtitle("genes, replicates=T")
genes.MDS.rep<-MDSplot(genes(cuff), replicates=TRUE) + ggtitle("genes, replicates=T")

print(genes.PCA.rep)
print(genes.MDS.rep)


genes.PCA1<-PCAplot(isoforms(cuff),"PC1","PC2",replicates=TRUE) + ggtitle("ISOFORMS, replicates=T")
genes.PCA1
genes.PCA2<-PCAplot(CDS(cuff),"PC1","PC2",replicates=TRUE) + ggtitle("CDS, replicates=T")
genes.PCA2
genes.PCA3<-PCAplot(TSS(cuff),"PC1","PC2",replicates=TRUE) + ggtitle("TSS, replicates=T")
genes.PCA3
genes.PCA4<-PCAplot(isoforms(cuff),"PC1","PC2",replicates=F) + ggtitle("ISOFORMS, replicates=F")
genes.PCA4
genes.PCA5<-PCAplot(CDS(cuff),"PC1","PC2",replicates=F) + ggtitle("CDS, replicates=F")
genes.PCA5
genes.PCA6<-PCAplot(TSS(cuff),"PC1","PC2",replicates=F) + ggtitle("TSS, replicates=F")
genes.PCA6

genes.MDS<-MDSplot(genes(cuff),replicates=FALSE)

geners.PCA1.rep<-MDSplot(genes(cuff),"PC1","PC2",replicates=F) + ggtitle("GENES, replicates=F")
geners.PCA1.rep
genes.MDS1.rep<-MDSplot(sigGenes,replicates=TRUE) + ggtitle("SIGGENES, replicates=T")
genes.MDS1.rep
enes.MDS2.rep<-MDSplot(isoforms(sigGenes), replicates=TRUE) + ggtitle("SIG_ISOFORMS, replicates=T")
enes.MDS2.rep
enes.MDS3.rep<-MDSplot(CDS(sigGenes), replicates=TRUE) + ggtitle("SIG_CDS, replicates=T")
enes.MDS3.rep
enes.MDS4.rep<-MDSplot(TSS(sigGenes), replicates=TRUE) + ggtitle("SIG_TSS, replicates=T")
enes.MDS4.rep

################################################################################################
###  Create data frames of sig DE genes, CDS, TSS, isoform....
################################################################################################
mySigGeneIds<-getSig(cuff,x=over,y=under,alpha=.05,level='genes')
mySigGeneIds
#mySigGeneIds2<-getSig(cuff,x=over,y=under,alpha=.05,level='isoforms')
#mySigGeneIds2
mySigGeneIds3<-getSig(cuff,x=over,y=under,alpha=.05,level='CDS')
#mySigGeneIds3
#mySigGeneIds4<-getSig(cuff,x=over,y=under,alpha=.05,level='TSS')
#mySigGeneIds4
#mySigGeneIds5<-getSig(cuff,x=over,y=under,alpha=.05,level='relCDS')
#mySigGeneIds5
#mySigGeneIds1<-getSig(cuff,x=over,y=under,alpha=.05,level='splicing')
#mySigGeneIds1



###################################################
### code chunk number 84: geneset_cluster_1
###################################################
# ic<-csCluster(myGenes,k=4)
# head(ic$cluster)
# icp<-csClusterPlot(ic)
# icp

###################################################
### code chunk number 86: specificity_1
###################################################
# myGenes.spec<-csSpecificity(myGenes)
# head(myGenes.spec)


## ###################################################
## ### code chunk number 87: similar_1
## ###################################################
# mySimilar<-findSimilar(cuff,myGene,n=20)
# mySimilar.expression<-expressionPlot(mySimilar,logMode=T,showErrorbars=F)


## ###################################################
## ### code chunk number 88: similar_plots_1
## ###################################################
## mySimilar<-findSimilar(cuff,"PINK1",n=20)
## mySimilar.expression<-expressionPlot(mySimilar,logMode=T,showErrorbars=F)
## print(mySimilar.expression)


###################################################
### code chunk number 89: similar_2
###################################################
# myProfile<-c(500,0,400)
# mySimilar2<-findSimilar(cuff,myProfile,n=10)
# mySimilar2.expression<-expressionPlot(mySimilar2,logMode=T,showErrorbars=F)
# mySimilar2.expression
# print(mySimilar2.expression)

###################################################
### code chunk number 90: similar_plots_2
###################################################
# mySimilar2.expression<-expressionPlot(mySimilar2,logMode=F,showErrorbars=F)
# mySimilar2.expression
# print(mySimilar2.expression)


###################################################
### code chunk 91: close_connection  92: session
###################################################
end<-dbDisconnect(cuff@DB)
sessionInfo()


