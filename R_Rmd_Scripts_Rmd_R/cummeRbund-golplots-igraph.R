### CummeRbund-GO-Graph-Enrichment
###############################################################################
# GOplots for Significantly differentially expressed features: method 1
############################################################################
## source('openGraphSaveGraph.r');
options(error=traceback) # causes a traceback to appear if there is an error, and the traceback has the line number
# prefixed by #

### original  R code from Drew'
## options(echo=TRUE)
## args = commandArgs(trailingOnly = T)
args = commandArgs(TRUE)

for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
}

###################################################
### code chunk number 30: data_access_0
###################################################

# These arguments are passed in on the command line via launchR.eachDiffDir.sh
# args:
# diffDir=\"${diffDir}\" inDir=\"${INPUT}/${diffDir}\" outDir=\"${OUTPUT}/${diffDir}\"
# under=\"${under}\" over=\"${over}\" Rplots=\"$\{OUTPUT}/${diffDir}/${R_PLOTS}\"
# FPKMmatrix=\"${OUTPUT}/${diffDir}/${FPKM_MATRIX}\" DiffTable=\"${OUTPUT}/${diffDir}/${DIFF_TABLE}\"

###################################################
# The parameters below  are for use on the cluster
# outDir=\"${OUTPUT}/${diffDir}\"
#outDir="/media/drew/67a0cf75-b688-4ee0-ac64-b8bbbb27d0c6/home/drew/umb_triley/urine1/cuffdiff_results_hg19_default/LUTS-over-CTRL"
#myDir<-file.path(outDir)
#myDir
###################################################
# The parameters below  are for local running
myDir<-file.path("/media/drew/67a0cf75-b688-4ee0-ac64-b8bbbb27d0c6/home/drew/umb_triley/urine1/cuffdiff_results_hg19_default/LUTS-over-CTRL")
setwd(myDir)
dir()

###################################################
### code chunk number 1: loadLib
###################################################
library(cummeRbund);library(edgeR)
library(keggorthology);library(limma)
library(gage);library(gageData)
library(igraph);library(KEGGgraph)
library(org.Hs.eg.db);library(biomaRt)
library(reactome.db);library(ReactomePA)
library(clusterProfiler);library(DOSE)
library(GOstats);library(topGO)
library(KEGG.db);library(GO.db)
#library(GSEABase);library(pathview)
#library(STRINGdb); source("http://igraph.sf.net")
###################################################
### code chunk number 2: read
###################################################
#gtfFile<-file.path("/media/drew/67a0cf75-b688-4ee0-ac64-b8bbbb27d0c6/home/drew/umb_triley/urine1/cuffd")
#cuff <- readCufflinks(dir=myDir,gtfFile=gtfFile,genome="hg19",rebuild=T)
cuff <- readCufflinks(dir='.',genome="hg38",rebuild=T)
cuff

rna.seq.sum<-runInfo(cuff)
RNAseq.run.info<-rna.seq.sum[1,2]
RNAseq.run.info
replicates.info<-replicates(cuff)
replicates.info
groups<-factor(replicates.info$sample_name)
samples<-replicates.info$rep_name
groups
samples
over="LUTS"
under="CTRL"

###################################################
### code chunk number 31: data_access_1
###################################################

gene_exp.diff<-diffData(genes(cuff))
head(gene_exp.diff)
dim(gene_exp.diff)

g.rep.matrix<-repFpkmMatrix(isoforms(cuff))
head(g.rep.matrix)

gene.rep.matrix<-as.data.frame(repFpkmMatrix(isoforms(cuff)))
row.names(gene.rep.matrix)->tracking_id

gene.rep.ma<-cbind(tracking_id, gene.rep.matrix)
head(gene.rep.ma)

gene.xloc.matrix<-featureNames(isoforms(cuff))
head(gene.xloc.matrix)

gene.list<-getGenes(cuff,geneId = gene.xloc.matrix)
gene.list

gene_annotation_data<-featureNames(gene.list@isoforms)
head(gene_annotation_data)

gene.features<-annotation(isoforms(cuff))
head(gene.features)

iso_exp.diff<-diffData(isoforms(cuff))
head(iso_exp.diff)
dim(iso_exp.diff)

mySigGenes<-getSig(cuff,x=over,y=under,alpha=.05,level='genes')
head(mySigGenes)
length(mySigGenes)
sigGenes<-getGenes(cuff, mySigGenes)
length(sigGenes)

sig_genes_annot<-sigGenes@annotation
head(sig_genes_annot)

sig_genes_exp.diff<-diffData(sigGenes)
sig.genes_exp.diff<-as.data.frame(cbind(logFC=sig_genes_exp.diff$log2_fold_change,q_value=sig_genes_exp.diff$q_value),
											 row.names=sig_genes_exp.diff$gene_id)

head(sig.genes_exp.diff)
gene.features<-gene.features[,c(1,2,4)]
head(gene.features)

g.exp.df<-as.data.frame(merge(x=gene.features, y=gene.rep.ma, by.x="isoform_id", by.y="tracking_id"))#,row.names = "gene_short_name")
head(g.exp.df)

G.exp.df<-as.data.frame(g.exp.df[,-c(1:2)],row.names = g.exp.df[,3])
g.rep.matrix<-G.exp.df[,-1]
head(g.rep.matrix)
dim(g.rep.matrix)

geneList<-mapIds(x = org.Hs.eg.db,
                 keys =  row.names(g.rep.matrix),
                 column = "ENTREZID",
                 keytype = "SYMBOL",
                 multiVals="first")

head(g.rep.matrix[unique(row.names(g.rep.matrix)),])
g.rep.ma<-g.rep.matrix[unique(row.names(g.rep.matrix)),]
head(g.rep.ma)
g.rep.ma<-as.matrix(g.rep.ma)

over.group<-grep(pattern = over, x = colnames(g.rep.ma),ignore.case = T)
under.group<-grep(pattern = under, x = colnames(g.rep.ma),ignore.case = T)

g.under.matrix<-g.rep.ma[,under.group]
head(g.under.matrix)
dim(g.under.matrix)

g.over.matrix<-g.rep.ma[,over.group]
head(g.over.matrix)
dim(g.over.matrix)

####################################################################################################
### code chunk number : Step 1, Read in, normalise, and identify genes with significant effects
####################################################################################################

sig_over_gene_data<-subset(gene_diff_data, (significant =='yes') & (log2_fold_change > 0))
head(sig_over_gene_data)
nrow(sig_over_gene_data)

sig_under_gene_data<-subset(gene_diff_data, (significant =='yes') & (log2_fold_change < 0))
head(sig_under_gene_data)
nrow(sig_under_gene_data)

####################################################################
### code chunk number 2: Write repFPKMMatrix
####################################################################

g.rep.matrix<-repFpkmMatrix(genes(cuff))
head(g.rep.matrix)
rep.fpkm.file = file.path("repFPKMmatrix")
write.table(g.rep.matrix,file = rep.fpkm.file, sep = "  ", row.names = T, col.names = T,quote = F)

o.group<-grep(pattern = over,colnames(g.rep.matrix),ignore.case = T)
u.group<-grep(pattern = under,colnames(g.rep.matrix),ignore.case = T)

sig.g.rep.matrix<-g.rep.matrix[which(row.names(g.rep.matrix) %in% sig_gene_data$gene_id),]
head(sig.g.rep.matrix)
dim(sig.g.rep.matrix)
sig.rep.fpkm.file = file.path("sig.repFPKMmatrix")
write.table(sig.g.rep.matrix,file = sig.rep.fpkm.file, sep = "  ", row.names = T, col.names = T,quote = F)

over.group<-grep(pattern = over,colnames(sig.g.rep.matrix),ignore.case = T)
under.group<-grep(pattern = under,colnames(sig.g.rep.matrix),ignore.case = T)

####################################################################################
### code chunk number 3: Write sig.rep.overFPKMMatrix and sig.rep.underFPKMMatrix
###################################################################################

sig.g.o.rep.matrix<-g.rep.matrix[which(row.names(g.rep.matrix) %in% sig_over_gene_data$gene_id),]
head(sig.g.o.rep.matrix)
dim(sig.g.o.rep.matrix)
sig.o.rep.fpkm.file = file.path("sig.over.repFPKMmatrix")
write.table(sig.g.o.rep.matrix,file = sig.o.rep.fpkm.file, sep = "  ", row.names = T, col.names = T,quote = F)
over.g<-grep(pattern = over,colnames(sig.g.o.rep.matrix),ignore.case = T)
sig.g.o.rep.ma<-sig.g.o.rep.matrix[over.g]
head(sig.g.o.rep.ma)
dim(sig.g.o.rep.ma)

sig.g.u.rep.matrix<-g.rep.matrix[which(row.names(g.rep.matrix) %in% sig_under_gene_data$gene_id),]
head(sig.g.u.rep.matrix)
dim(sig.g.u.rep.matrix)
sig.u.rep.fpkm.file = file.path("sig.under.repFPKMmatrix")
write.table(sig.g.u.rep.matrix,file = sig.u.rep.fpkm.file, sep = "  ", row.names = T, col.names = T,quote = F)
under.g<-grep(pattern = over,colnames(sig.g.u.rep.matrix),ignore.case = T)
sig.g.u.rep.ma<-sig.g.u.rep.matrix[under.g]
head(sig.g.u.rep.ma)
dim(sig.g.u.rep.ma)

######################################################################
### code chunk number 3: Write repCOUNTMatrix and sig.repCOUNTmatrix
######################################################################

rep.gene.counts<-repCountMatrix(isoforms(cuff))
head(rep.gene.counts)
dim(rep.gene.counts)

sig.rep.gene.counts<-rep.gene.counts[which(row.names(rep.gene.counts) %in% sig_gene_data$gene_id),]
head(sig.rep.gene.counts)
dim(sig.rep.gene.counts)

rep.count.file = file.path("repCOUNTmatrix")
write.table(rep.gene.counts,file = rep.count.file, sep = "  ", row.names = T, col.names = T,quote = F)

sig.rep.count.file = file.path("sig.repCOUNTmatrix")
write.table(sig.rep.gene.counts,file = sig.rep.count.file, sep = "  ", row.names = T, col.names = T,quote = F)

######################################################################
### code chunk number 3: Write repCOUNTMatrix and sig.repCOUNTmatrix
######################################################################
lot(gse116969.sd,main="Standard Deviation of Rh8_Rh5 and Rh8_Rh6 ")
plot(Rh8mean.expr,main="Mean gene expression of Rh8_Rh5 and Rh8_Rh6 ")
heatmap(sig.genes.gse116969)

sig.genes.gse116969.df<-as.data.frame(sig.genes.gse116969)
sig.gse.cor<-cor(sig.genes.gse116969.df)
#plotcorr(sig.gse.cor, mar = c(0.1, 0.1, 0.1, 0.1))
# Do the same, but with colors corresponding to value
colorfun <- colorRamp(c("#CC0000","white","#3366CC"), space="Lab")
plotcorr(sig.gse.cor, col=rgb(colorfun((sig.gse.cor+1)/2), maxColorValue=255),
         mar = c(0.1, 0.1, 0.1, 0.1))


sig.dist<-as.matrix(dist(t(sig.genes.gse116969)))
loc <- cmdscale(sig.dist) # Performs MDS analysis on the geographic distances between European cities.
plot(loc[,1], -loc[,2], type="n", xlab="", ylab="", main="Euclidian MDS analysis")
text(loc[,1], -loc[,2], rownames(loc), cex=1,col=c("blue","darkgreen"))
# Plots the MDS results in 2D plot. The minus is required in this example to flip the plotting orientation.

mydatascale <- t(scale(t(sig.genes.gse116969))) # Centers and scales data.
hr <- hclust(as.dist(1-cor(t(mydatascale), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
hc <- hclust(as.dist(1-cor(mydatascale, method="spearman")), method="complete")
# Clusters columns by Spearman correlation.
heatmap(sig.genes.gse116969, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), scale="row") # col=my.colorFct(),
# Plot the data table as heatmap and the cluster results as dendrograms.
mycl <- cutree(hr, h=max(hr$height)/1.5)
mycolhc <- sample(rainbow(256))
mycolhc <- mycolhc[as.vector(mycl)]
heatmap(sig.genes.gse116969, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), scale="row", RowSideColors=mycolhc) #col=my.colorFct(),
# Cut the tree at specific height and color the corresponding clusters in the heatmap color bar.

#sig.dist.t<-as.matrix(cor(t(sig.genes.gse116969)))
#fit <- cmdscale(sig.dist.t,eig=TRUE, k=2) # k is the number of dim
#fit # view results
# plot solution
#x <- fit$points[,1];y <- fit$points[,2]
#plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
#     main="Metric MDS", type="n")
#text(x, y, labels = row.names(sig.dist.t),cex=.7,col=c("blue","darkgreen"))

gse.scaled<-as.matrix(scale(t(sig.genes.gse116969)))
heatmap(gse.scaled, Colv=F, scale='none')

gsecor<-cor(sig.genes.gse116969)
round(gsecor,2)
# Do the same, but with colors corresponding to value
plot(hclust(dist(sig.genes.gse116969,method="euclidian"),method="single"))
cl <- kmeans(sig.genes.gse116969, 2)
plot(sig.genes.gse116969, col = cl$cluster)
points(cl$centers, col = 1:2, pch = 8, cex=2)

gsecor<-cor(t(sig.genes.gse116969))
round(gsecor,2)
# Do the same, but with colors corresponding to value
plot(hclust(dist(sig.genes.gse116969,method="euclidian"),method="single"))
cl <- kmeans(sig.genes.gse116969, 2)
plot(sig.genes.gse116969, col = cl$cluster)
points(cl$centers, col = 1:2, pch = 8, cex=2)


c <- cor(t(sig.genes.gse116969), method="spearman")
d <- as.dist(1-c)
hr <- hclust(d, method = "complete", members=NULL)
# This is a short example of clustering the rows of the generated sample matrix 'y' with hclust. Seven different cluster joining
# methods can be selected with the 'method' argument: ward, single, complete, average, mcquitty, median and centroid. The
# object returned by hclust is a list of class hclust which describes the tree generated by the clustering process with the
# following components: merge, height, order, labels, method, call and dist.method.
par(mfrow = c(2, 2))
plot(hr, hang = 0.1)
plot(hr, hang = -1)
# The generated tree can be plotted with the plot() function. When the hang argument is set to '-1' then all leaves end on
# one line and their labels hang down from 0. More details on the plotting behavior is provided in the hclust help document (?hclust).
plot(as.dendrogram(hr), edgePar=list(col=3, lwd=4), horiz=T,)
# To plot trees horizontally, the hclust tree has to be transformed into a dendrogram object.
par(mfrow=c(2,1))
unclass(hr) # Prints the full content of the hclust object.
str(as.dendrogram(hr)) # Prints dendrogram structure as text.
hr$labels[hr$order] # Prints the row labels in the order they appear in the tree.
hrd1 <- as.dendrogram(hr); plot(hrd1)
hrd2 <- reorder(hrd1, sample(1:10))
plot(hrd2)
labels(hrd1)
labels(hrd2)
##################################################
### code chunk number : Step 6: Statistical Tests
##################################################
t.test(x = Rh5mean.expr,
       y = Rh6mean.expr,
       alternative = "two.sided")
wilcox.test(Rh5mean.expr,Rh6mean.expr, alternative = "g")  # g for greater
var.test(RH8expr[,"Rh5mean.expr"],RH8expr[,"Rh6mean.expr"])
chisq.test(Rh5mean.expr,Rh6mean.expr)
dev.off()
####################################################################################################
### code chunk number : Step 3: #Hierarchical clustering
####################################################################################################
library(pvclust); library(gplots) # Loads the required packages.
#Hierarchical clustering
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R")
# Import an alternative color scheme for the heatmap function.
mydatascale <- t(scale(t(g.over.matrix))) # Centers and scales data.
hr <- hclust(as.dist(1-cor(t(g.over.matrix), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
hc <- hclust(as.dist(1-cor(g.over.matrix, method="spearman")), method="complete")
# Clusters columns by Spearman correlation.
heatmap(mydatascale, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=my.colorFct(), scale="row")
# Plot the data table as heatmap and the cluster results as dendrograms.
mycl <- cutree(hr, h=max(hr$height)/1.5)
mycolhc <- sample(rainbow(256))
mycolhc <- mycolhc[as.vector(mycl)]
heatmap(g.over.matrix, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),
        col=my.colorFct(), scale="row", RowSideColors=mycolhc)
# Cut the tree at specific height and color the corresponding clusters in the heatmap color bar.

pv <- pvclust(scale(t(g.over.matrix)), method.dist="correlation", method.hclust="complete", nboot=10)
# Perform the hierarchical cluster analysis.
#Due to time resrictions, we are using here only 10 bootstrap repetitions.
# Usually, one should use at least 1000 repetitions.
plot(pv, hang=-1); pvrect(pv, alpha=0.95)
# Plots result as a dendrogram where the significant clusters are highlighted with red rectangles.
clsig <- unlist(pvpick(pv, alpha=0.95, pv="au", type="geq", max.only=TRUE)$clusters) # Retrieve members of significant clusters.

source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/dendroCol.R") # Import tree coloring function.
dend_colored <- dendrapply(as.dendrogram(pv$hclust), dendroCol, keys=clsig, xPar="edgePar", bgr="black", fgr="red", pch=20)
# Create dendrogram object where the significant clusters are labeled in red.
heatmap(sig.genes.gse116969, Rowv=dend_colored, Colv=as.dendrogram(hc), col=my.colorFct(), scale="row", RowSideColors=mycolhc)
# Plot the heatmap from above, but with the significant clusters in red and the cluster bins from the tree cutting step in
# the color bar.
x11(height=12); heatmap.2(mydata, Rowv=dend_colored, Colv=as.dendrogram(hc), col=my.colorFct(), scale="row", trace="none", RowSideColors=mycolhc) # Plot heatmap with heatmap.2() function which scales better for many entries.
mydatasort <- mydata[pv$hclust$labels[pv$hclust$order], hc$labels[hc$order]] # Sort rows in data table by 'dend_colored' and its colums by 'hc'.
x11(height=16, width=12); par(mfrow=c(1,2)); plot(dend_colored, horiz=T, yaxt="n"); image(scale(t(mydatasort)), col=my.colorFct(), xaxt="n",yaxt="n") # Plot heatmap with bootstrap tree in larger format using instead of heatmap the image function.
pdf("pvclust.pdf", height=21, width=6); plot(dend_colored, horiz=T, yaxt="n"); dev.off(); pdf("heatmap.pdf", height=20, width=6); image(scale(t(mydatasort)), col=my.colorFct(), xaxt="n",yaxt="n"); dev.off()



####################################################################################################
### code chunk number : Step 3: #Hierarchical clustering
####################################################################################################

library(pvclust); library(gplots) # Loads the required packages.

#Hierarchical clustering

source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R")
# Import an alternative color scheme for the heatmap function.
mydatascale <- t(scale(t(g.over.matrix))) # Centers and scales data.
hr <- hclust(as.dist(1-cor(t(mydatascale), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
hc <- hclust(as.dist(1-cor(mydatascale, method="spearman")), method="complete")
# Clusters columns by Spearman correlation.
heatmap(g.rep.ma, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=my.colorFct(), scale="row")
# Plot the data table as heatmap and the cluster results as dendrograms.
mycl <- cutree(hr, h=max(hr$height)/1.5)
mycolhc <- sample(rainbow(256))
mycolhc <- mycolhc[as.vector(mycl)]
heatmap(sig.g.rep.matrix, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),
        col=my.colorFct(), scale="row", RowSideColors=mycolhc)
# Cut the tree at specific height and color the corresponding clusters in the heatmap color bar.

pv <- pvclust(scale(t(sig.g.rep.matrix)), method.dist="correlation", method.hclust="complete", nboot=10)
# Perform the hierarchical cluster analysis.
#Due to time resrictions, we are using here only 10 bootstrap repetitions.
# Usually, one should use at least 1000 repetitions.
plot(pv, hang=-1); pvrect(pv, alpha=0.95)
# Plots result as a dendrogram where the significant clusters are highlighted with red rectangles.
clsig <- unlist(pvpick(pv, alpha=0.95, pv="au", type="geq", max.only=TRUE)$clusters) # Retrieve members of significant clusters.

source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/dendroCol.R") # Import tree coloring function.
dend_colored <- dendrapply(as.dendrogram(pv$hclust), dendroCol, keys=clsig, xPar="edgePar", bgr="black", fgr="red", pch=20)
# Create dendrogram object where the significant clusters are labeled in red.
heatmap(sig.g.rep.matrix, Rowv=dend_colored, Colv=as.dendrogram(hc), col=my.colorFct(), scale="row", RowSideColors=mycolhc)
# Plot the heatmap from above, but with the significant clusters in red and the cluster bins from the tree cutting step in
# the color bar.
x11(height=12); heatmap.2(mydata, Rowv=dend_colored, Colv=as.dendrogram(hc), col=my.colorFct(), scale="row", trace="none", RowSideColors=mycolhc) # Plot heatmap with heatmap.2() function which scales better for many entries.
mydatasort <- mydata[pv$hclust$labels[pv$hclust$order], hc$labels[hc$order]] # Sort rows in data table by 'dend_colored' and its colums by 'hc'.
x11(height=16, width=12); par(mfrow=c(1,2)); plot(dend_colored, horiz=T, yaxt="n"); image(scale(t(mydatasort)), col=my.colorFct(), xaxt="n",yaxt="n") # Plot heatmap with bootstrap tree in larger format using instead of heatmap the image function.
pdf("pvclust.pdf", height=21, width=6); plot(dend_colored, horiz=T, yaxt="n"); dev.off(); pdf("heatmap.pdf", height=20, width=6); image(scale(t(mydatasort)), col=my.colorFct(), xaxt="n",yaxt="n"); dev.off()


####################################################################################################
### code chunk number : Step 2, Read in, normalise, and identify genes with significant effects
####################################################################################################

gene_diff_data<-diffData(isoforms(cuff))
head(gene_diff_data)
dim(gene_diff_data)

mySigGenes<-getSig(cuff,x=over,y=under,alpha=.05,level='genes')
head(mySigGenes)
length(mySigGenes)
sigGenes<-getGenes(cuff, mySigGenes)
length(sigGenes)
sig_genes_annot<-sigGenes@annotation
head(sig_genes_annot)

sig_genes_exp.diff<-diffData(sigGenes)
sig.genes_exp.diff<-as.data.frame(cbind(logFC=sig_genes_exp.diff$log2_fold_change,q_value=sig_genes_exp.diff$q_value),
											 row.names=sig_genes_exp.diff$gene_id)

head(sig.genes_exp.diff)
sig_gene_data<-subset(sig.genes_exp.diff, q_value < 0.005 & abs(logFC) > 2)
head(sig_gene_data)
dim(sig_gene_data)
genes<-row.names(sig_gene_data)
fold.change<-sig_gene_data[,"logFC"]
sig.gene.List<-as.data.frame(fold.change, row.names = genes)
head(sig.gene.List)

sig.gene.order<-order(fold.change, decreasing = T)
fold<-sig.gene.List[sig.gene.order,]
genes<-rownames(sig.gene.List)

sig_genes<-which(sig.genes_exp.diff[,"q_value"] < 0.005 & abs(sig.genes_exp.diff[,"logFC"]) > 2)
length(sig_genes)
head(sig.genes_exp.diff[sig_genes,])
logfold=sig.genes_exp.diff[sig_genes,c("logFC")]
sig.gene.List<-as.data.frame(logfold, row.names = genes)
head(sig.gene.List)

sorted.fold.vals<-order(logfold, decreasing = T)
gene.eg.fc<-genes[sorted.fold.vals]
sorted.fc<-logfold[sorted.fold.vals]
sig.gene.fold.list<-as.data.frame(sorted.fc,row.names = gene.eg.fc)
head(sig.gene.fold.list)
dim(sig.gene.fold.list)

de <- row.names(sig.gene.fold.list)[abs(sig.gene.fold.list) > 1.5]
head(de)
length(de)

eg = bitr(de, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
sym.eg<-as.factor(eg$SYMBOL)
dim(eg)

symbols<-row.names(sig.gene.fold.list)
head(symbols)
length(symbols)

entrez.syms<-which(row.names(sig.gene.fold.list) %in% sym.eg)
length(entrez.syms)
head(sig.gene.fold.list)
foldchange<-sig.gene.fold.list[entrez.syms,]
SYMBOL<-symbols[entrez.syms]
siggene.fc.eg<-cbind(SYMBOL, foldchange)
head(siggene.fc.eg)
siggene.eg.fc<-merge(eg, siggene.fc.eg, by.x="SYMBOL", by.y="SYMBOL")
SigGeneList<-as.data.frame(siggene.eg.fc[,"foldchange"], row.names=siggene.eg.fc$ENTREZID)
colnames(SigGeneList)<-"log2.fold.change"
head(SigGeneList)
head(row.names(SigGeneList))
de <- SigGeneList
head(de)
x <- enrichPathway(gene=SigGeneList,pvalueCutoff=0.05, readable=T)
head(as.data.frame(x))

## ----fig.height=6, fig.width=12------------------------------------------
barplot(x, showCategory=8)
dotplot(x, showCategory=15)
cnetplot(x, categorySize="pvalue", foldChange=fold.change)
SigGeneList<-SigGeneList[order(SigGeneList$log2.fold.change, decreasing = T),]
## ------------------------------------------------------------------------
y <- gsePathway(SigGeneList, nPerm=10000,
                pvalueCutoff=0.2,
                pAdjustMethod="BH", verbose=FALSE)
res <- as.data.frame(y); head(res)

## ----fig.height=8, fig.width=8-------------------------------------------
#emapplot(y, color="pvalue")
#gseaplot(y, geneSetID = "R-HSA-69242")
#viewPathway("E2F mediated regulation of DNA replication", readable=TRUE, foldChange=geneList)
gene.df <- bitr(names(geneList), fromType = "SYMBOL",toType = c("ENTREZID", "SYMBOL"),OrgDb = org.Hs.eg.db)

gene.df <- bitr(gene_diff_data$gene_id, fromType = "SYMBOL",toType = c("ENTREZID", "SYMBOL"),OrgDb = org.Hs.eg.db)
head(gene.df)
ggoCC <- groupGO(gene     = de,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

ggoBP <- groupGO(gene     = de,
                 OrgDb    = org.Hs.eg.db,
                 ont      = "BP",
                 level    = 3,
                 readable = TRUE)

ggoMF <- groupGO(gene     = de,
                 OrgDb    = org.Hs.eg.db,
                 ont      = "MF",
                 level    = 3,
                 readable = TRUE)
head(ggo)
## ----fig.height=5, fig.width=9-------------------------------------------
barplot(ggoCC, drop=TRUE, showCategory=12)
barplot(ggoMF, drop=TRUE, showCategory=12)
barplot(ggoBP, drop=TRUE, showCategory=12)

barplot(ego, showCategory=8)
dotplot(ego)

## ----fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8----
gseaplot(kk2, geneSetID = "hsa04145")

## ------------------------------------------------------------------------
ego <- enrichGO(gene          = de,
                universe      = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
cnetplot(ego, categorySize="pvalue")

## categorySize can be scaled by 'pvalue' or 'geneNum'
cnetplot(ego, categorySize="pvalue", foldChange=geneList)

## ----fig.height=12, fig.width=8, eval=FALSE------------------------------
goplot(ego)

## ----eval=FALSE----------------------------------------------------------
#  ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
#                  OrgDb         = org.Hs.eg.db,
#                  keyType       = 'ENSEMBL',
#                  ont           = "CC",
#                  pAdjustMethod = "BH",
#                  pvalueCutoff  = 0.01,
#                  qvalueCutoff  = 0.05)

## ----eval=FALSE----------------------------------------------------------
#  ego2 <- setReadable(ego2, OrgDb = org.Hs.eg.db)

## ----eval=FALSE----------------------------------------------------------
  egoCC <- gseGO(geneList     = de,
                OrgDb        = org.Hs.eg.db,
                ont          = "CC",
                nPerm        = 1000,
                minGSSize    = 100,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                verbose      = FALSE)
egoBP <- gseGO(geneList     = de,
               OrgDb        = org.Hs.eg.db,
               ont          = "BP",
               nPerm        = 1000,
               minGSSize    = 100,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
egoMF <- gseGO(geneList     = de,
               OrgDb        = org.Hs.eg.db,
               ont          = "MF",
               nPerm        = 1000,
               minGSSize    = 100,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
## ------------------------------------------------------------------------
kk <- enrichKEGG(gene         = de,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)

## ------------------------------------------------------------------------
kk2 <- gseKEGG(geneList     = de,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)

## ----eval = FALSE--------------------------------------------------------
#  mkk <- enrichMKEGG(gene = gene,
#                     organism = 'hsa')

## ----eval=FALSE----------------------------------------------------------
#  mkk2 <- gseMKEGG(geneList = geneList,
#                   species = 'hsa')

## ----eval=FALSE----------------------------------------------------------
#  david <- enrichDAVID(gene = gene,
#                       idType = "ENTREZ_GENE_ID",
#                       listType = "Gene",
#                       annotation = "KEGG_PATHWAY",
#                       david.user = "clusterProfiler@hku.hk")

## ------------------------------------------------------------------------
search_kegg_organism('ece', by='kegg_code')
ecoli <- search_kegg_organism('Escherichia coli', by='scientific_name')
dim(ecoli)
head(ecoli)

gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)

egmt <- enricher(gene, TERM2GENE=c5)
head(egmt)

egmt2 <- GSEA(geneList, TERM2GENE=c5, verbose=FALSE)
head(egmt2)

## ----fig.height=5, fig.width=9-------------------------------------------
barplot(ggo, drop=TRUE, showCategory=12)

## ----fig.height=5, fig.width=8-------------------------------------------
barplot(ego, showCategory=8)

## ------------------------------------------------------------------------
dotplot(ego)

## ----fig.cap="enrichment map of enrichment result", fig.align="center", fig.height=8, fig.width=8----
emapplot(ego)

## ----fig.height=8, fig.width=8-------------------------------------------
## categorySize can be scaled by 'pvalue' or 'geneNum'
cnetplot(ego, categorySize="pvalue", foldChange=geneList)

## ----fig.height=12, fig.width=8, eval=FALSE------------------------------
#  goplot(ego)

## ----fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8----
gseaplot(kk, geneSetID = "hsa04145")

## ----eval=FALSE----------------------------------------------------------
#  browseKEGG(kk, 'hsa04110')

## ----eval=FALSE----------------------------------------------------------
#  library("pathview")
#  hsa04110 <- pathview(gene.data  = geneList,
#                       pathway.id = "hsa04110",
#                       species    = "hsa",
#                       limit      = list(gene=max(abs(geneList)), cpd=1))

## ------------------------------------------------------------------------
data(gcSample)
lapply(gcSample, head)

## ------------------------------------------------------------------------
ck <- compareCluster(geneCluster = gcSample, fun = "enrichKEGG")
head(as.data.frame(ck))

## ------------------------------------------------------------------------
mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"

formula_res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichKEGG")
head(as.data.frame(formula_res))
## ----fig.height=7, fig.width=9-------------------------------------------
dotplot(ck)
dotplot(formula_res)
dotplot(formula_res, x=~group) + ggplot2::facet_grid(~othergroup)


########################################################################
### code chunk number 4: Write adj mtrix and plot graph for over genes
########################################################################
g.under.matrix<-g.rep.ma[,under.group]
head(g.under.matrix)
dim(g.under.matrix)

g.over.matrix<-g.rep.ma[,over.group]
head(g.over.matrix)
dim(g.over.matrix)

#Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
over.graph <- graph.adjacency(as.matrix(as.dist(cor(t(g.over.matrix), method="pearson"))),
                              mode="undirected", weighted=TRUE, diag=FALSE)

#NB - Euclidean distances can also be used instead of correlation, but this changes the tutorial slightly:
#g <- graph.adjacency(as.matrix(dist(estrogenMainEffects, method="euclidean")),
                      #mode="undirected", weighted=TRUE, diag=FALSE)
#Simplfy the adjacency object
over.graph <- simplify(over.graph, remove.multiple=TRUE, remove.loops=TRUE)
#Colour negative correlation edges as blue
E(over.graph)[which(E(over.graph)$weight<0)]$color <- "darkblue"
#Colour positive correlation edges as red
E(over.graph)[which(E(over.graph)$weight>0)]$color <- "darkred"
#Convert edge weights to absolute values
E(over.graph)$weight <- abs(E(over.graph)$weight)
#Change arrow size
#For directed graphs only
#E(over.graph)$arrow.size <- 1.0

#Remove edges below absolute Pearson correlation 0.8
over.graph <- delete_edges(over.graph, E(over.graph)[which(E(over.graph)$weight<0.9)])
#Assign names to the graph vertices (optional)
V(over.graph)$name <- V(over.graph)$name
#Change shape of graph vertices
V(over.graph)$shape <- "sphere"
#Change colour of graph vertices
V(over.graph)$color <- "skyblue"
#Change colour of vertex frames
V(over.graph)$vertex.frame.color <- "white"

#Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
#Multiple scaled vales by a factor of 10
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
vSizes <- (scale01(apply(sig.g.rep.matrix, 1, mean)) + 1.0) * 10

#Amplify or decrease the width of the edges
edgeweights <- E(over.graph)$weight * 2.0
#Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst <- mst(over.graph, algorithm="prim")

#Plot the tree object
plot(mst, layout=layout.fruchterman.reingold, edge.curved=TRUE,
     vertex.size=vSizes, vertex.label.dist=-0.5, vertex.label.color="black",
     asp=FALSE, vertex.label.cex=0.6, edge.width=edgeweights, edge.arrow.mode=0,
     main="Network of significant genes")

########################################################################
### code chunk number 5: Write adj mtrix and plot graph for under genes
########################################################################

#Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
under.graph <- graph.adjacency(as.matrix(as.dist(cor(t(sig.g.u.rep.ma), method="pearson"))),
                              mode="undirected", weighted=TRUE, diag=FALSE)

#NB - Euclidean distances can also be used instead of correlation, but this changes the tutorial slightly:
#g <- graph.adjacency(as.matrix(dist(estrogenMainEffects, method="euclidean")),
#mode="undirected", weighted=TRUE, diag=FALSE)
#Simplfy the adjacency object
under.graph <- simplify(under.graph, remove.multiple=TRUE, remove.loops=TRUE)
#Colour negative correlation edges as blue
E(under.graph)[which(E(under.graph)$weight<0)]$color <- "darkblue"
#Colour positive correlation edges as red
E(under.graph)[which(E(under.graph)$weight>0)]$color <- "darkred"
#Convert edge weights to absolute values
E(under.graph)$weight <- abs(E(under.graph)$weight)
#Change arrow size
#For directed graphs only
#E(over.graph)$arrow.size <- 1.0

#Remove edges below absolute Pearson correlation 0.8
under.graph <- delete_edges(under.graph, E(under.graph)[which(E(under.graph)$weight<0.8)])
#Assign names to the graph vertices (optional)
V(under.graph)$name <- V(under.graph)$name
#Change shape of graph vertices
V(under.graph)$shape <- "sphere"
#Change colour of graph vertices
V(under.graph)$color <- "skyblue"
#Change colour of vertex frames
V(under.graph)$vertex.frame.color <- "white"

#Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
#Multiple scaled vales by a factor of 10
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
vSizes <- (scale01(apply(sig.g.u.rep.ma, 1, mean)) + 1.0) * 10

#Amplify or decrease the width of the edges
edgeweights <- E(under.graph)$weight * 2.0
#Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst <- mst(under.graph, algorithm="prim")

#Plot the tree object
plot(mst, layout=layout.fruchterman.reingold, edge.curved=TRUE,
     vertex.size=vSizes, vertex.label.dist=-0.5, vertex.label.color="black",
     asp=FALSE, vertex.label.cex=0.6, edge.width=edgeweights, edge.arrow.mode=0,
     main=under "patient significant gene network")

#################################################################################################
### code chunk number 4: Identify communities in the tree object based on 'edge betweenness'
###############################################################################################

mst.communities <- edge.betweenness.community(mst, directed=FALSE)
mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
V(mst)$color <- mst.communities$membership + 1

par(mfrow=c(1,2))
plot(
  mst.clustering, mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="LUTS patient significant gene network"
)

plot(
  mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="LUTS patient significant gene network"
)

#################################################################################################
### code chunk number 5: Step 4, Further analyses
###############################################################################################

#Check the vertex degree, i.e., number of connections to each vertex
degree(mst)

#Output information for each community, including vertex-to-community assignments and modularity
commSummary <- data.frame(
  mst.communities$names,
  mst.communities$membership,
  mst.communities$modularity
)
colnames(commSummary) <- c("Gene", "Community", "Modularity")
options(scipen=999)
commSummary

#Compare community structures using the variance of information (vi) metric (not relevant here)
#Community structures that are identical will have a vi=0
compare(mst.communities, mst.communities, method="vi")


dim(sig.genes_exp.diff)
head(sig.genes_exp.diff)

genes<-row.names(sig.genes_exp.diff)
foldchange<-sig.genes_exp.diff[,"logFC"]
qval<-sig.genes_exp.diff[,"q_value"]

under.inf = which(foldchange == "-Inf" & qval < 0.05 | foldchange < -1 & foldchange != "-Inf" & qval < 0.05)

over.inf = which(foldchange == "Inf" & qval < 0.05 | foldchange > 1 &  foldchange != "Inf" & qval < 0.05)


HIexp.inOVER<-as.data.frame(sig.genes_exp.diff[over.inf,])
head(HIexp.inOVER)
HI.exp.OVER.genes<-row.names(HIexp.inOVER)
length(HI.exp.OVER.genes)
Qval.hi_exprOVER<-HIexp.inOVER[,"q_value"]
names(Qval.hi_exprOVER)<-HI.exp.OVER.genes
#HI.exp.OVER.genes<-row.names(Qval.hi_exprOVER)
logfc.hi_exprOVER<-HIexp.inOVER[,"logFC"]
names(logfc.hi_exprOVER)<-HI.exp.OVER.genes

HIexp.inUNDER<-as.data.frame(sig.genes_exp.diff[under.inf,])
head(HIexp.inUNDER)
HI.exp.UNDER.genes<-row.names(HIexp.inUNDER)
length(HI.exp.UNDER.genes)
Qval.hi_exprUNDER<-HIexp.inUNDER[,"q_value"]
names(Qval.hi_exprUNDER)<-HI.exp.UNDER.genes
#HI.exp.UNDER.genes<-row.names(Qval.hi_exprUNDER)
logfc.hi_exprUNDER<-HIexp.inUNDER[,"logFC"]
names(logfc.hi_exprUNDER)<-HI.exp.UNDER.genes


#QvalnoOVER #QvalnoUNDER
ENTREZQval.hi_exprOVER<-mapIds(x = org.Hs.eg.db,
                         keys =  names(Qval.hi_exprOVER),
                         column = "ENTREZID",
                         keytype = "SYMBOL",
                         multiVals="first")

ENTREZQval.hi_exprUNDER<-mapIds(x = org.Hs.eg.db,
                          keys =  names(Qval.hi_exprUNDER),
                          column = "ENTREZID",
                          keytype = "SYMBOL",
                          multiVals="first")

ENTREZsiggenes<-mapIds(x = org.Hs.eg.db,
                       keys = sig_genes_exp.diff$gene_id,
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals="first")
ENTREZsiggenes

#-------------------------------------------------------------------------------------------------------

xx <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
topDiffGenes <- function(allScore) {
  return(allScore < 0.01)}

#-------------------------------------------------------------------------------------------------------

# Significantly Differentially Expressed by the Under group
# i.e Luts over Ctrl --> Ctrl = Under group, Under group highly
# expresses and No expression in the Over group

HIunder.BP.GOdata <- new("topGOdata",ontology = "BP",allGenes = Qval.hi_exprUNDER, nodeSize = 5,
	       	  			     annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
					     	       ID = "symbol")
HIunder.MF.GOdata <- new("topGOdata",ontology = "MF",allGenes = Qval.hi_exprUNDER, nodeSize = 5,
                        annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                        ID = "symbol")
HIunder.CC.GOdata <- new("topGOdata",ontology = "CC",allGenes = Qval.hi_exprUNDER, nodeSize = 5,
                        annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                        ID = "symbol")

HIover.BP.GOdata <- new("topGOdata",ontology = "BP",allGenes = Qval.hi_exprOVER, nodeSize = 5,
		    			       annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
					       	         ID = "symbol")
HIover.MF.GOdata <- new("topGOdata",ontology = "MF",allGenes = Qval.hi_exprOVER, nodeSize = 5,
	       	  			     annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
					     	       ID = "symbol")

HIover.CC.GOdata <- new("topGOdata",ontology = "CC",allGenes = Qval.hi_exprOVER, nodeSize = 5,
                      annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                      ID = "symbol")

HIunder.BPtKS <- runTest(HIunder.BP.GOdata, algorithm = "classic", statistic = "ks")
HIunder.BPtKS
HIunder.BPFisher <- runTest(HIunder.BP.GOdata, algorithm = "classic", statistic = "fisher")
HIunder.BPFisher
HIunder.BPtKS.elim <- runTest(HIunder.BP.GOdata, algorithm = "elim", statistic = "ks")
HIunder.BPtKS.elim

HIunder.MFtKS <- runTest(HIunder.MF.GOdata, algorithm = "classic", statistic = "ks")
HIunder.MFtKS
HIunder.MFFisher <- runTest(HIunder.MF.GOdata, algorithm = "classic", statistic = "fisher")
HIunder.MFFisher
HIunder.MFFtKS.elim <- runTest(HIunder.MF.GOdata, algorithm = "elim", statistic = "ks")
HIunder.MFFtKS.elim

HIunder.CCtKS <- runTest(HIunder.CC.GOdata, algorithm = "classic", statistic = "ks")
HIunder.CCtKS
HIunder.CCFisher <- runTest(HIunder.CC.GOdata, algorithm = "classic", statistic = "fisher")
HIunder.CCFisher
HIunder.CCtKS.elim <- runTest(HIunder.CC.GOdata, algorithm = "elim", statistic = "ks")
HIunder.CCtKS.elim

pdf("GOplots_grch38.CTRL_highly_expressed.pdf")

showSigOfNodes(HIunder.BP.GOdata, score(HIunder.BPFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(HIunder.BP.GOdata, score(HIunder.BPtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(HIunder.BP.GOdata, score(HIunder.BPtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(HIunder.BP.GOdata, HIunder.BPFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(HIunder.BP.GOdata, HIunder.BPtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(HIunder.BP.GOdata, HIunder.BPtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(HIunder.MF.GOdata, score(HIunder.MFFisher), firstSigNodes = 5, useInfo = "all" ,)
showSigOfNodes(HIunder.MF.GOdata, score(HIunder.MFtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(HIunder.MF.GOdata, score(HIunder.MFFtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(HIunder.MF.GOdata, HIunder.MFFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(HIunder.MF.GOdata, HIunder.MFFtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(HIunder.MF.GOdata, HIunder.MFtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(HIunder.CC.GOdata, score(HIunder.CCFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(HIunder.CC.GOdata, score(HIunder.CCtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(HIunder.CC.GOdata, score(HIunder.CCtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(HIunder.CC.GOdata, HIunder.CCFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(HIunder.CC.GOdata, HIunder.CCtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(HIunder.CC.GOdata, HIunder.CCtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

o.group<-grep(pattern = over, x = colnames(g.count.matrix),ignore.case = T)
u.group<-grep(pattern = under, x = colnames(g.count.matrix),ignore.case = T)

over.group<-grep(pattern = over, x = colnames(g.rep.matrix),ignore.case = T)
under.group<-grep(pattern = under, x = colnames(g.rep.matrix),ignore.case = T)

isexpr <- rowSums(cpm(g.count.matrix) > 1) >= 2
gse116969<-gse116969[isexpr,]
dim(gse116969)

sig.g.reps<-which(row.names(g.rep.matrix) %in% sigGenes@ids)
sig.g.rep.matrix<-as.matrix(g.rep.matrix[sig.g.reps,])
head(sig.g.rep.matrix)
dim(sig.g.rep.matrix)

min.samps<-(length(samps)/2)*.75
min.expr<-round(min.samps)

isexpr <- rowSums(cpm(g.count.matrix[,o.group]) > 1) >= min.expr | rowSums(cpm(g.count.matrix[,u.group]) > 1) >= min.expr
is.expr <- rowSums(cpm(g.rep.matrix[,over.group]) > 1) >= min.expr | rowSums(cpm(g.rep.matrix[,under.group]) > 1) >= min.expr
sig.g.rep.ma<-g.rep.matrix[is.expr,]
dim(sig.g.rep.ma)
head(sig.g.rep.ma)


RNA.Design = data.frame(row.names = replicates.info$rep_name,
                       condition = replicates.info$sample_name)

conditions<-factor(RNA.Design$condition)
RNA.Design
conditions

OverGroup = RNA.Design$condition == over
UnderGroup = RNA.Design$condition == under

Rh5samples.gene_exp = gse116969[,Rh5group]
Rh6samples.gene_exp = gse116969[,Rh6group]
Rh5.1<-Rh5samples.gene_exp[,"R8_Rh5_d1_rep1"]
Rh5.2<-Rh5samples.gene_exp[,"R8_Rh5_d1_rep2"]

###############################################################################
# Map gene symbols to entrezids in order to perform GSEA and GO analysis
################################################################################
if(require(org.Hs.eg.db)) {
	tnodes <- nodes(toprSub)
	tgeneids <- translateKEGGID2GeneID(tnodes)
	tgenesymbols <- sapply(mget(tgeneids, org.Hs.egSYMBOL, ifnotfound=NA), "[[",1)
	toprSubSymbol <- toprSub
	nodes(toprSubSymbol) <- tgenesymbols
	plot(toprSubSymbol, "neato",attrs=list(node=list(font=5, fillcolor="lightblue")))


library(biomaRt)
hsapiens <- useMart("ensembl","hsapiens_gene_ensembl" )
filters <- listFilters(hsapiens)
getBM(attributes=c("entrezgene","hgnc_symbol"),
      filters="entrezgene",
      values=toprbccGeneID, mart=hsapiens)

columns(org.Hs.eg.db)

geneList<-mapIds(x = org.Hs.eg.db,
                 keys =  rownames(iso.rep.matrix),
                 column = "ENTREZID",
                 keytype = "REFSEQ",
                 multiVals="first")

gl<-as.data.frame(geneList)
head(gl)

xx <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
topDiffGenes <- function(allScore) {
  return(allScore < 0.01)}

# Significantly Differentially Expressed by the Under group
# i.e Luts over Ctrl --> Ctrl = Under group, Under group highly
# expresses and No expression in the Over group
loOVERbpGOdata <- new("topGOdata",ontology = "BP",allGenes = genes_exp.diff, nodeSize = 5,
                      annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                      ID = "symbol")

loOVERbpGOdata <- new("topGOdata",ontology = "BP",allGenes = Qval.hi_exprUNDER, nodeSize = 5,
                      annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                      ID = "symbol")
HIunder.MFGOdata <- new("topGOdata",ontology = "MF",allGenes = Qval.hi_exprUNDER, nodeSize = 5,
                        annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                        ID = "symbol")
HIunder.CCGOdata <- new("topGOdata",ontology = "CC",allGenes = Qval.hi_exprUNDER, nodeSize = 5,
                        annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                        ID = "symbol")

LOunder.BPGOdata <- new("topGOdata",ontology = "BP",allGenes = Qval.hi_exprOVER, nodeSize = 5,
                        annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                        ID = "symbol")
hiOVERmfGOdata <- new("topGOdata",ontology = "MF",allGenes = Qval.hi_exprOVER, nodeSize = 5,
                      annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                      ID = "symbol")

hiOVERccGOdata <- new("topGOdata",ontology = "CC",allGenes = Qval.hi_exprOVER, nodeSize = 5,
                      annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                      ID = "symbol")

LOover.BPtKS <- runTest(loOVERbpGOdata, algorithm = "classic", statistic = "ks")
LOover.BPFisher <- runTest(loOVERbpGOdata, algorithm = "classic", statistic = "fisher")
LOover.BPtKS.elim <- runTest(loOVERbpGOdata, algorithm = "elim", statistic = "ks")

LOunder.BPtKS <- runTest(LOunder.BPGOdata, algorithm = "classic", statistic = "ks")
LOunder.BPFisher <- runTest(LOunder.BPGOdata, algorithm = "classic", statistic = "fisher")
LOunder.BPtKS.elim <- runTest(LOunder.BPGOdata, algorithm = "elim", statistic = "ks")

HIover.MFtKS <- runTest(hiOVERmfGOdata, algorithm = "classic", statistic = "ks")
HIover.MFFisher <- runTest(hiOVERmfGOdata, algorithm = "classic", statistic = "fisher")
HIoverMFtKS.elim <- runTest(hiOVERmfGOdata, algorithm = "elim", statistic = "ks")

HIunder.MFtKS <- runTest(HIunder.MFGOdata, algorithm = "classic", statistic = "ks")
HIunder.MFFisher <- runTest(HIunder.MFGOdata, algorithm = "classic", statistic = "fisher")
HIunder.MFtKS.elim <- runTest(HIunder.MFGOdata, algorithm = "elim", statistic = "ks")

HIover.CCtKS <- runTest(hiOVERccGOdata, algorithm = "classic", statistic = "ks")
HIover.CCFisher <- runTest(hiOVERccGOdata, algorithm = "classic", statistic = "fisher")
HIover.CCtKS.elim <- runTest(hiOVERccGOdata, algorithm = "elim", statistic = "ks")

HIunder.CCtKS <- runTest(HIunder.CCGOdata, algorithm = "classic", statistic = "ks")
HIunder.CCFisher <- runTest(HIunder.CCGOdata, algorithm = "classic", statistic = "fisher")
HIunder.CCtKS.elim <- runTest(HIunder.CCGOdata, algorithm = "elim", statistic = "ks")

pdf("GOplots_grch38.CTRL_highly_expressed.pdf")

showSigOfNodes(loOVERbpGOdata, score(LOover.BPFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(loOVERbpGOdata, score(LOover.BPtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(loOVERbpGOdata, score(LOover.BPtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(loOVERbpGOdata, LOover.BPFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(loOVERbpGOdata, LOover.BPtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(loOVERbpGOdata, LOover.BPtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(HIunder.MFGOdata, score(HIunder.MFFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(HIunder.MFGOdata, score(HIunder.MFtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(HIunder.MFGOdata, score(HIunder.MFtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(HIunder.MFGOdata, HIunder.MFFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(HIunder.MFGOdata, HIunder.MFtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(HIunder.MFGOdata, HIunder.MFtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(HIunder.CCGOdata, score(HIunder.CCFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(HIunder.CCGOdata, score(HIunder.CCtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(HIunder.CCGOdata, score(HIunder.CCtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(HIunder.CCGOdata, HIunder.CCFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(HIunder.CCGOdata, HIunder.CCtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(HIunder.CCGOdata, HIunder.CCtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

dev.off()
pdf("GOplots_grch38.LUTS_highly_expressed.pdf")

showSigOfNodes(LOunder.BPGOdata, score(LOunder.BPFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(LOunder.BPGOdata, score(LOunder.BPtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(LOunder.BPGOdata, score(LOunder.BPtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(LOunder.BPGOdata, LOunder.BPFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(LOunder.BPGOdata, LOunder.BPtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(LOunder.BPGOdata, LOunder.BPtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(hiOVERmfGOdata, score(HIover.MFFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(hiOVERmfGOdata, score(HIover.MFtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(hiOVERmfGOdata, score(HIoverMFtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(hiOVERmfGOdata, HIover.MFFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(hiOVERmfGOdata, HIoverMFtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(hiOVERmfGOdata, HIover.MFtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(hiOVERccGOdata, score(HIover.CCFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(hiOVERccGOdata, score(HIover.CCtKS), firstSigNodes = 5, useInfo = "def" )
showSigOfNodes(hiOVERccGOdata, score(HIover.CCtKS.elim), firstSigNodes = 5, useInfo = "def" )
printGraph(hiOVERccGOdata, HIover.CCFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(hiOVERccGOdata, HIover.CCtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(hiOVERccGOdata, HIover.CCtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

###################################################
### code chunk number 38: create_geneset_3
###################################################
My.Sig.Genes<-getSig(cuff,x=over,y=under,alpha=.001,level='genes')
MY.sig.Genes<-getGenes(cuff, My.Sig.Genes)

sigIsos<-getSig(cuff,x=over,y=under,alpha=.01,level='isoforms')
mySigIsos<-getGenes(cuff, sigIsos)

sigTSS<-getSig(cuff,x=over,y=under,alpha=.01,level='TSS')
mySigTSS<-getGenes(cuff, sigTSS)

sigCDS<-getSig(cuff,x=over,y=under,alpha=.01,level='CDS')
mySigCDS<-getGenes(cuff, sigCDS)

mySigrelCDS<-getSig(cuff,x=over,y=under,alpha=.01,level='relCDS')
sigrelCDS<-getGenes(cuff, mySigrelCDS)

mySigPromoter<-getSig(cuff,x=over,y=under,alpha=.01,level='promoters')
head(mySigPromoter)
length(mySigPromoter)
sigPromoter<-getGenes(cuff, mySigPromoter)
length(sigPromoter)

mySigSplice<-getSig(cuff,x=over,y=under,alpha=.01,level='splicing')
head(mySigSplice)
length(mySigSplice)
sigSplice<-getGenes(cuff, mySigSplice)
length(sigSplice)

plot(csDendro(MY.sig.Genes, replicates=T,logMode = T),
     main="Dendrogram built from the Significantly Diff-Expr Genes at alpha=0.01")
plot(csDendro(mySigIsos, replicates=T,logMode = T),
     main="Dendrogram built from the Significantly Diff-Expr Isoforms at alpha=0.01")
plot(csDendro(mySigTSS, replicates=T,logMode = T),
     main="Dendrogram built from the Significantly Diff-Expr TSS at alpha=0.01")
plot(csDendro(mySigCDS, replicates=T,logMode = T),
     main="Dendrogram built from the Significantly Diff-Expr CDS at alpha=0.01")
plot(csDendro(sigrelCDS, replicates=T,logMode = T),
     main="Dendrogram built from the Significantly Differentially loaded CDS features at alpha=0.01")
plot(csDendro(sigPromoter, replicates=T,logMode = T),
     main="Dendrogram built from the Significantly Differentially loaded Promoters at alpha=0.01")
plot(csDendro(sigSplice, replicates=T,logMode = T),
     main="Dendrogram built from the Significantly Differentially loaded Splicing at alpha=0.01")
plot(dispersionPlot(cuff) + ggtitle("Gene Level Dispersion"))


plot(dispersionPlot(isoforms(cuff)) + ggtitle("ISoform Level Dispersion"))
plot(dispersionPlot(CDS(cuff)) + ggtitle("CDS Level Dispersion"))
plot(dispersionPlot(TSS(cuff)) + ggtitle("TSS Level Dispersion"))

###################################################
### code chunk number 9: SCV_visualization
###################################################

plot(fpkmSCVPlot(isoforms(cuff)))
plot(fpkmSCVPlot(TSS(cuff)))
plot(fpkmSCVPlot(CDS(cuff)))

plot(csDensity(isoforms(cuff),replicates=T))
plot(csDensity(CDS(cuff),replicates=T))
plot(csDensity(TSS(cuff),replicates=T))


###################################################
### code chunk number 13: global_plots_2
###################################################
plot(csBoxplot(isoforms(cuff),logMode = T, replicates=T))
plot(csBoxplot(TSS(cuff),logMode = T,replicates=T))
plot(csBoxplot(CDS(cuff),logMode = T, replicates=T))

csSM<-csScatterMatrix(isoforms(cuff),logMode = T,replicates = T)
plot(csScatterMatrix(TSS(cuff),logMode = T,replicates = T))
plot(csScatterMatrix(CDS(cuff),logMode = T,replicates = T))

plot(csScatter(isoforms(cuff),,logMode = T,replicates = T))
plot(csScatter(TSS(cuff),logMode = T,replicates = T))
plot(csScatte(CDS(cuff),logMode = T,replicates = T))


###################################################
### code chunk number 17: global_plots_scatter_1
###################################################

sigGOcc <- groupGO(gene=ENTREZQval.hi_exprUNDER, OrgDb=org.Hs.eg.db, ont="CC",
                   level    = 3,readable = TRUE)
sigGOmf <- groupGO(gene=ENTREZQval.hi_exprUNDER, OrgDb=org.Hs.eg.db, ont="MF",
                   level    = 3,readable = TRUE)
sigGObp <- groupGO(gene=ENTREZQval.hi_exprUNDER, OrgDb=org.Hs.eg.db, ont="BP",
                   level    = 3,readable = TRUE)
barplot(sigGOcc, drop=TRUE, showCategory=12,args.legend="Sig Diff Expr Genes Grouping by GO-CC")
barplot(sigGOmf, drop=TRUE, showCategory=12,args.legend="Sig Diff Expr Genes Grouping by GO-MF")
barplot(sigGObp, drop=TRUE, showCategory=12,args.legend="Sig Diff Expr Genes Grouping by GO-BP")

###################################################
### code chunk number 19: global_plots_scatter_2
###################################################

sigGOcc <- groupGO(gene=ENTREZQval.hi_exprOVER, OrgDb=org.Hs.eg.db, ont="CC",
                   level    = 3,readable = TRUE)
sigGOmf <- groupGO(gene=ENTREZQval.hi_exprOVER, OrgDb=org.Hs.eg.db, ont="MF",
                   level    = 3,readable = TRUE)
sigGObp <- groupGO(gene=ENTREZQval.hi_exprOVER, OrgDb=org.Hs.eg.db, ont="BP",
                   level    = 3,readable = TRUE)
barplot(sigGOcc, drop=TRUE, showCategory=12,args.legend="Sig Diff Expr Genes Grouping by GO-CC")
barplot(sigGOmf, drop=TRUE, showCategory=12,args.legend="Sig Diff Expr Genes Grouping by GO-MF")
barplot(sigGObp, drop=TRUE, showCategory=12,args.legend="Sig Diff Expr Genes Grouping by GO-BP")

###################################################
### code chunk number 27: global_plots_volcano_1
###################################################
v<-csVolcanoMatrix(isoforms(cuff),replicates = T,logMode = T)
v
v<-csVolcanoMatrix(CDS(cuff),replicates = T,logMode = T)
v
v<-csVolcanoMatrix(TSS(cuff),replicates = T,logMode = T)
v

###################################################
### code chunk number 28: global_plots_5_2
###################################################
v<-csVolcano(isoforms(cuff),over, under, alpha=0.01, showSignificant=TRUE)
v
v<-csVolcano(CDS(cuff),over, under, alpha=0.01, showSignificant=TRUE)
v
v<-csVolcano(TSS(cuff),over, under, alpha=0.01, showSignificant=TRUE)
v

###################################################
### code chunk number 48: geneset_plots_4
###################################################
ih<-csHeatmap(isoforms(myGenes),cluster='both',labRow=F)
ih
th<-csHeatmap(TSS(myGenes),cluster='both',labRow=F)
th
th<-csHeatmap(CDS(myGenes),cluster='both',labRow=F)
th

###################################################
### code chunk number 73: get_sig_3
###################################################
mySigGenes<-getGenes(cuff,mySigGeneIds)
mySigGenes

###################################################
### code chunk number 75: dist_heat_1
###################################################
myDistHeat<-csDistHeat(isoforms(cuff))
print(myDistHeat)

myDistHeat<-csDistHeat(CDS(cuff))
print(myDistHeat)

myDistHeat<-csDistHeat(TSS(cuff))
print(myDistHeat)

###################################################
### code chunk number 76: dist_heat_plot_1
###################################################
myRepDistHeat<-csDistHeat(isoforms(cuff),replicates=T)
print(myRepDistHeat)

myRepDistHeat<-csDistHeat(TSS(cuff),replicates=T)
print(myRepDistHeat)

myRepDistHeat<-csDistHeat(CDS(cuff),replicates=T)
print(myRepDistHeat)


###################################################
### code chunk number 79: dim_reduction_1
###################################################
genes.PCA<-PCAplot(isoforms(cuff), replicates = T)
genes.PCA

genes.PCA<-PCAplot(CDS(cuff), replicates = T)
genes.PCA

genes.PCA<-PCAplot(TSS(cuff), replicates = T)
genes.PCA

genes.MDS.rep<-MDSplot(isoforms(cuff),replicates=T)
genes.MDS.rep

genes.MDS.rep<-MDSplot(isoforms(cuff),replicates=T)
genes.MDS.rep

genes.MDS.rep<-MDSplot(isoforms(cuff),replicates=T)
genes.MDS.rep

genes.PCA.rep<-PCAplot(genes(cuff),"PC1","PC2",replicates=T)
genes.MDS.rep<-MDSplot(genes(cuff),replicates=T)


###################################################
### code chunk number 80: gene_PCA
###################################################
genes.PCA<-PCAplot(genes(cuff),"PC1","PC2")
genes.MDS<-MDSplot(genes(cuff))

genes.PCA.rep<-PCAplot(genes(cuff),"PC1","PC2",replicates=T)
genes.MDS.rep<-MDSplot(genes(cuff),replicates=T)
	print(genes.PCA)

###################################################
### code chunk number 84: geneset_cluster_1
###################################################
ic<-csCluster(myGenes,k=4)
head(ic$cluster)
icp<-csClusterPlot(ic)
icp
print(icp)
myGenes.spec<-csSpecificity(myGenes)
head(myGenes.spec)

###################################################
### code chunk number 88: similar_plots_1
###################################################

mySimilar<-findSimilar(cuff,"CXCL12",n=20)
mySimilar
mySimilar.expression<-expressionPlot(mySimilar,logMode=T,showErrorbars=T)
print(mySimilar.expression)

mySimilar<-findSimilar(cuff,"CXCR4",n=20)
mySimilar.expression<-expressionPlot(mySimilar,logMode=T,showErrorbars=T)
print(mySimilar.expression)

myProfile<-c(800,0,100)
mySimilar2<-findSimilar(cuff,myProfile,n=10)
mySimilar2.expression<-expressionPlot(mySimilar2,logMode=T,showErrorbars=T)
print(mySimilar2.expression)

myProfile<-c(600,0,200)
mySimilar2<-findSimilar(cuff,myProfile,n=10)
mySimilar2.expression<-expressionPlot(mySimilar2,logMode=T,showErrorbars=T)
print(mySimilar2.expression)


###################################################
### code chunk number 91: close_connection
###################################################
end<-dbDisconnect(cuff@DB)
sessionInfo()

