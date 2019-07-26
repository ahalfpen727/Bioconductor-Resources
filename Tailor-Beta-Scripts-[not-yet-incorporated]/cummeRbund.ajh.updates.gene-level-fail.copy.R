##############################################################################################
### CummeRbund Tailor Script gene level fail
##############################################################################################
##### To install old version of packages:
#source("https://bioconductor.org/biocLite.R")
#biocLite(); install.packages("devtools")
#install_version("RSQLite", version = "1.1-2", repos = "http://cran.us.r-project.org")
#require(devtools);library(spliceR)
##############################################################################################
### Library+Variable --> 1) Load 2) import and print env variables 3) initialize cuff-obj
##############################################################################################
options(width=65)
library(cummeRbund)
#source('openGraphSaveGraph.R');

# causes a traceback to appear if there is an error, and the traceback has the line number prefixed by {#}
options(error=traceback)

#options(echo=TRUE)
#args = commandArgs(trailingOnly = T)

# Gets arguments that were passed in via command line
args = commandArgs(TRUE)

for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
}
print(args)
pdf(file.path(Rplots))

cuff<-readCufflinks(dir=inDir,gtfFile=refgtf,genome=genomeR,out=outDir,rebuild=T)
cuff

##############################################################################################
# code chunk for debugging purposes, must have all cuffdiff output (from *-over-* directory)
# define gtf and genome files also
#############################################################################################
#save.image(file.path("cuffData.hg19_default.RData"))
# These arguments are passed in on the command line via launchR.eachDiffDir.sh
# args:
# diffDir=\"${diffDir}\" inDir=\"${INPUT}/${diffDir}\" outDir=\"${OUTPUT}/${diffDir}\"  ,
# under=\"${under}\" over=\"${over}\" Rplots=\"$\{OUTPUT}/${diffDir}/${R_PLOTS}\",
# FPKMmatrix=\"${OUTPUT}/${diffDir}/${FPKM_MATRIX}\" DiffTable=\"${OUTPUT}/${diffDir}/${DIFF_TABLE}\"
diffDir<-file.path("./urine1/cuffdiff_results_hg38_default/LUTS-over-CTRL")
setwd(diffDir);dir()
over="LUTS";under="CTRL"
cuff=readCufflinks(dir=".",rebuild=F,verbose = T)
cuff<-cuff.hg19
cuff.hg19

runInfo(cuff.hg19)
run.info<-replicates(cuff)
conditions<-factor(run.info$sample_name)
samples<-factor(run.info$rep_name)

###########################################################################
### Chunk_2: Significantly differentially expressed features: method 1
############################################################################

genes_exp.diff<-diffData(genes(cuff))
head(genes_exp.diff);dim(genes_exp.diff)
g.rep.matrix<-repFpkmMatrix(genes(cuff))
head(g.rep.matrix);dim(g.rep.matrix)
g.cnt.matrix<-repCountMatrix(genes(cuff))
head(g.cnt.matrix);dim(g.cnt.matrix)

fpkm.file = file.path("FPKMmatrix")
write.table(g.rep.matrix,file = fpkm.file, sep = "  ", row.names = F, col.names = T,quote = F)
fpkm.file = file.path("COUNTmatrix")
write.table(g.cnt.matrix,file = fpkm.file, sep = "  ", row.names = F, col.names = T,quote = F)

over.group<-grep(pattern=over, colnames(g.rep.matrix))
under.group<-grep(pattern=under, colnames(g.rep.matrix))
g.o.rep.matrix<-g.rep.matrix[,over.group]
g.u.rep.matrix<-g.rep.matrix[,under.group]
g.o.cnt.matrix<-g.cnt.matrix[,over.group]
g.u.cnt.matrix<-g.cnt.matrix[,under.group]

gene.xloc.matrix<-featureNames(genes(cuff))
head(gene.xloc.matrix);length(gene.xloc.matrix)

gene.annotation<-getGenes(cuff,geneId = gene.xloc.matrix)
gene.annotation;length(gene.annotation)

mySigGenes<-getSig(cuff,x=over,y=under,alpha=.05,level='genes')
length(mySigGenes); head(mySigGenes)

sigGenes<-getGenes(cuff, mySigGenes)
sigGenes; length(sigGenes)

sig_genes_exp.diff<-diffData(sigGenes)
head(sig_genes_exp.diff);dim(sig_genes_exp.diff)
sig_genes_exp.diff<-subset(sig_genes_exp.diff, abs(sig_genes_exp.diff$log2_fold_change) > 2 & sig_genes_exp.diff$q_value < 0.05)
head(sig_genes_exp.diff);dim(sig_genes_exp.diff)

sig.g.rep.matrix<-g.rep.matrix[which(row.names(g.rep.matrix) %in% sig_genes_exp.diff$gene_id),]
head(sig.g.rep.matrix);dim(sig.g.rep.matrix)
sig.g.o.rep.matrix<-sig.g.rep.matrix[,over.group]
sig.g.u.rep.matrix<-sig.g.rep.matrix[,under.group]
head(sig.g.o.rep.matrix);dim(sig.g.o.rep.matrix)
head(sig.g.u.rep.matrix);dim(sig.g.u.rep.matrix)

sig.g.cnt.matrix<-g.cnt.matrix[which(row.names(g.cnt.matrix) %in% sig_genes_exp.diff$gene_id),]
head(sig.g.cnt.matrix);dim(sig.g.cnt.matrix)
sig.g.o.cnt.matrix<-sig.g.cnt.matrix[,over.group]
sig.g.u.cnt.matrix<-sig.g.cnt.matrix[,under.group]


sig.o.genes_exp.diff<-subset(sig_genes_exp.diff, sig_genes_exp.diff$log2_fold_change > 2 & sig_genes_exp.diff$q_value < 0.01)
head(sig.o.genes_exp.diff);dim(sig.o.genes_exp.diff)

sig.u.genes_exp.diff<-subset(sig_genes_exp.diff, sig_genes_exp.diff$log2_fold_change < 2 & sig_genes_exp.diff$q_value < 0.01)
head(sig.u.genes_exp.diff);dim(sig.u.genes_exp.diff)

s.o.g.rep.matrix<-g.rep.matrix[which(row.names(g.rep.matrix) %in% sig.o.genes_exp.diff$gene_id),]
head(s.o.g.rep.matrix);dim(s.o.g.rep.matrix)

s.u.g.rep.matrix<-g.rep.matrix[which(row.names(g.rep.matrix) %in% sig.u.genes_exp.diff$gene_id),]
head(s.u.g.rep.matrix);dim(s.u.g.rep.matrix)

over.grp.s.o.g.rep.ma<-s.o.g.rep.matrix[,over.group]
under.grp.s.o.g.rep.ma<-s.o.g.rep.matrix[,under.group]

over.grp.s.u.g.rep.ma<-s.u.g.rep.matrix[,over.group]
under.grp.s.u.g.rep.ma<-s.u.g.rep.matrix[,under.group]

s.o.g.cnt.matrix<-g.cnt.matrix[which(row.names(g.cnt.matrix) %in% sig.o.genes_exp.diff$gene_id),]
head(s.o.g.cnt.matrix);dim(s.o.g.cnt.matrix)

s.u.g.cnt.matrix<-g.cnt.matrix[which(row.names(g.cnt.matrix) %in% sig.u.genes_exp.diff$gene_id),]
head(s.u.g.cnt.matrix);dim(s.u.g.cnt.matrix)

over.grp.s.o.g.rep.ma<-s.o.g.cnt.matrix[,over.group]
under.grp.s.o.g.rep.ma<-s.o.g.cnt.matrix[,under.group]

over.grp.s.u.g.rep.ma<-s.u.g.cnt.matrix[,over.group]
under.grp.s.u.g.rep.ma<-s.u.g.cnt.matrix[,under.group]

###########################################################################
### Chunk_2.b: Significantly differentially expressed ISOFORMS
############################################################################
sig.genes.exp.df<-as.data.frame(cbind(log2FC=sig_genes_exp.diff$log2_fold_change,  q_value=sig_genes_exp.diff$q_value),row.names=sig_genes_exp.diff$gene_id)
head(sig.genes.exp.df)

s.o.g.rep.ma<-as.matrix(s.o.g.rep.matrix)
s.u.g.rep.ma<-as.matrix(s.u.g.rep.matrix)
over.grp.s.o.g.rep.ma<-as.matrix(s.o.g.rep.matrix[,over.group])
under.grp.s.o.g.rep.ma<-as.matrix(s.o.g.rep.matrix[,under.group])
over.grp.s.u.g.rep.ma<-as.matrix(s.u.g.rep.matrix[,over.group])
under.grp.s.u.g.rep.ma<-as.matrix(s.u.g.rep.matrix[,under.group])

head(s.o.g.rep.ma)
head(s.u.g.rep.ma)
head(over.grp.s.o.g.rep.ma)
head(under.grp.s.o.g.rep.ma)
head(over.grp.s.u.g.rep.ma)
head(over.grp.s.o.g.rep.ma)
head(under.grp.s.u.g.rep.ma)


heatmap(s.o.g.rep.ma)
heatmap(s.u.g.rep.ma)
heatmap(over.grp.s.o.g.rep.ma)
heatmap(under.grp.s.o.g.rep.ma)
heatmap(over.grp.s.u.g.rep.ma)
heatmap(under.grp.s.u.g.rep.ma)

cor(s.o.g.rep.ma)
s.o.g.rep.df<-as.data.frame(s.o.g.rep.ma)
sig.gse.cor<-cor(s.o.g.rep.df)
#plotcorr(sig.gse.cor,mar = c(0.1, 0.1, 0.1, 0.1))
# Do the same, but with colors corresponding to value
colorfun <- colorRamp(c("#CC0000","white","#3366CC"), space="Lab")
plotcorr(sig.gse.cor, col=rgb(colorfun((sig.gse.cor+1)/2), maxColorValue=255),
         mar = c(0.1, 0.1, 0.1, 0.1))

cor(s.u.g.rep.ma)
s.u.g.rep.df<-as.data.frame(s.u.g.rep.ma)
sig.gse.cor<-cor(s.u.g.rep.df)
#plotcorr(sig.gse.cor, mar = c(0.1, 0.1, 0.1, 0.1))
# Do the same, but with colors corresponding to value
colorfun <- colorRamp(c("#CC0000","white","#3366CC"), space="Lab")
plotcorr(sig.gse.cor, col=rgb(colorfun((sig.gse.cor+1)/2), maxColorValue=255),
         mar = c(0.1, 0.1, 0.1, 0.1))

cor(s.o.g.rep.ma)
s.o.g.rep.df<-as.data.frame(t(s.o.g.rep.ma))
sig.gse.cor<-cor(s.o.g.rep.df)
#plotcorr(sig.gse.cor,mar = c(0.1, 0.1, 0.1, 0.1))
# Do the same, but with colors corresponding to value
colorfun <- colorRamp(c("#CC0000","white","#3366CC"), space="Lab")
plotcorr(sig.gse.cor, col=rgb(colorfun((sig.gse.cor+1)/2), maxColorValue=255),
         mar = c(0.1, 0.1, 0.1, 0.1))


head(over.grp.s.o.g.rep.ma)
head(under.grp.s.u.g.rep.ma)

mean.cor<-function(x){Xsum<-apply(x,1,sum); x<-Xsum/length(x[1,])}

cor(s.u.g.rep.ma)
s.u.g.rep.df<-as.data.frame(t(s.u.g.rep.ma))
sig.gse.cor<-as.dist(1-s.u.g.rep.df)
sig.u.cor.mean<-mean.cor(sig.gse.cor)
head(sig.u.cor.mean)
head(sig.gse.cor)
sig.g.rep.matrix<-as.matrix(s.u.g.rep.ma)

mydatascale <- t(over.grp.s.o.g.rep.ma)# Centers and scales data.
hr <- cor(mydatascale, method="pearson") # Cluster rows by Pearson correlation.
dim(hr)
sig.ov.grp.o.rep.ma<-mean.cor(hr)
length(sig.ov.grp.o.rep.ma)
ns.ov.grp.o.rep.ma<- sig.ov.grp.o.rep.ma[which(sig.ov.grp.o.rep.ma > 0.8)]
names(ns.ov.grp.o.rep.ma)


sig.g.rep.matrix<-as.matrix(t(s.u.g.rep.ma))
head(sig.g.rep.matrix)
mydatascale <- sig.g.rep.matrix # Centers and scales data.
hr <- as.dist(1-cor(mydatascale, method="pearson"), method="complete") # Cluster rows by Pearson correlation.
hc <- hclust(as.dist(1-cor(mydatascale, method="spearman")), method="complete")
# Clusters columns by Spearman correlation.
heatmap(sig.g.rep.matrix, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), scale="row") # col=my.colorFct(),
# Plot the data table as heatmap and the cluster results as dendrograms.
mycl <- cutree(hr, h=max(hr$height)/1.5)
mycolhc <- sample(rainbow(256))
mycolhc <- mycolhc[as.vector(mycl)]
heatmap(sig.g.rep.matrix, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), scale="row", RowSideColors=mycolhc) #col=my.colorFct(),


apply(sig.gse.cor, 1, sum/length())
#plotcorr(sig.gse.cor, mar = c(0.1, 0.1, 0.1, 0.1))
# Do the same, but with colors corresponding to value
colorfun <- colorRamp(c("#CC0000","white","#3366CC"), space="Lab")
plotcorr(sig.gse.cor, col=rgb(colorfun((sig.gse.cor+1)/2), maxColorValue=255),
         mar = c(0.1, 0.1, 0.1, 0.1))


cor(over.grp.s.o.g.rep.ma)
over.grp.s.o.g.rep.df<-as.data.frame(over.grp.s.o.g.rep.ma)
sig.gse.cor<-cor(over.grp.s.o.g.rep.df)
#plotcorr(sig.gse.cor, mar = c(0.1, 0.1, 0.1, 0.1))
# Do the same, but with colors corresponding to value
colorfun <- colorRamp(c("#CC0000","white","#3366CC"), space="Lab")
plotcorr(sig.gse.cor, col=rgb(colorfun((sig.gse.cor+1)/2), maxColorValue=255),
         mar = c(0.1, 0.1, 0.1, 0.1))

dim(under.grp.s.u.g.rep.ma)
under.grp.s.u.g.rep.df<-as.data.frame(t(under.grp.s.u.g.rep.ma))
sig.gse.cor<-cor(under.grp.s.u.g.rep.df)
#plotcorr(sig.gse.cor, mar = c(0.1, 0.1, 0.1, 0.1))
# Do the same, but with colors corresponding to value
colorfun <- colorRamp(c("#CC0000","white","#3366CC"), space="Lab")
plotcorr(sig.gse.cor, col=rgb(colorfun((sig.gse.cor+1)/2), maxColorValue=255),
         mar = c(0.1, 0.1, 0.1, 0.1))

head(sig.g.o.rep.matrix);dim(sig.g.o.rep.matrix)
head(sig.g.u.rep.matrix);dim(sig.g.u.rep.matrix)

sig.dist<-as.matrix(dist(t(sig.g.rep.matrix)))
loc <- cmdscale(sig.dist) # Performs MDS analysis on the geographic distances between European cities.
plot(loc[,1], -loc[,2], type="n", xlab="", ylab="", main="Euclidian MDS analysis")
text(loc[,1], -loc[,2], rownames(loc), cex=1,col=c("blue","darkgreen"))
# Plots the MDS results in 2D plot. The minus is required in this example to flip the plotting orientation.
sig.g.rep.matrix<-as.matrix(sig.g.rep.matrix)
mydatascale <- t(scale(t(sig.g.rep.matrix))) # Centers and scales data.
hr <- hclust(as.dist(1-cor(t(mydatascale), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
hc <- hclust(as.dist(1-cor(mydatascale, method="spearman")), method="complete")
# Clusters columns by Spearman correlation.
heatmap(sig.g.rep.matrix, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), scale="row") # col=my.colorFct(),
# Plot the data table as heatmap and the cluster results as dendrograms.
mycl <- cutree(hr, h=max(hr$height)/1.5)
mycolhc <- sample(rainbow(256))
mycolhc <- mycolhc[as.vector(mycl)]
heatmap(sig.g.rep.matrix, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), scale="row", RowSideColors=mycolhc) #col=my.colorFct(),
# Cut the tree at specific height and color the corresponding clusters in the heatmap color bar.

#sig.dist.t<-as.matrix(cor(t(sig.genes.gse116969)))
#fit <- cmdscale(sig.dist.t,eig=TRUE, k=2) # k is the number of dim
#fit # view results
# plot solution
#x <- fit$points[,1];y <- fit$points[,2]
#plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
#     main="Metric MDS", type="n")
#text(x, y, labels = row.names(sig.dist.t),cex=.7,col=c("blue","darkgreen"))

gse.scaled<-as.matrix(scale(t(sig.g.rep.matrix)))
heatmap(gse.scaled, Colv=F, scale='none')

gsecor<-cor(sig.g.rep.matrix)
round(gsecor,2)
# Do the same, but with colors corresponding to value
plot(hclust(dist(sig.g.rep.matrix,method="euclidian"),method="single"))
cl <- kmeans(t(sig.g.rep.matrix),2 )
plot(sig.g.rep.matrix, col = cl$cluster)
points(cl$centers, col = 1:18, pch = 8, cex=2)

gsecor<-cor(t(sig.g.rep.matrix))
round(gsecor,2)
# Do the same, but with colors corresponding to value
plot(hclust(dist(sig.g.rep.matrix,method="euclidian"),method="single"))
cl <- kmeans(sig.g.rep.matrix, 2)
plot(sig.g.rep.matrix, col = cl$cluster)
points(cl$centers, col = 1:2, pch = 8, cex=2)


c <- cor(t(sig.g.rep.matrix), method="spearman")
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
mydatascale <- t(scale(t(sig.genes.gse116969))) # Centers and scales data.
hr <- hclust(as.dist(1-cor(t(sig.genes.gse116969), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
hc <- hclust(as.dist(1-cor(sig.genes.gse116969, method="spearman")), method="complete")
# Clusters columns by Spearman correlation.
heatmap(mydatascale, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=my.colorFct(), scale="row")
# Plot the data table as heatmap and the cluster results as dendrograms.
mycl <- cutree(hr, h=max(hr$height)/1.5)
mycolhc <- sample(rainbow(256))
mycolhc <- mycolhc[as.vector(mycl)]
heatmap(sig.genes.gse116969, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),
        col=my.colorFct(), scale="row", RowSideColors=mycolhc)
# Cut the tree at specific height and color the corresponding clusters in the heatmap color bar.

pv <- pvclust(scale(t(sig.genes.gse116969)), method.dist="correlation", method.hclust="complete", nboot=10)
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



mean.x.minus.xbar<-function(x){X.bar<-apply(x,1,mean); z<-abs(x-X.bar); apply(z, 1, mean)}

mad.g.o.rep.ma<-mean.x.minus.xbar(g.o.rep.matrix)
summary(mad.g.o.rep.ma)

mad.g.u.rep.ma<-mean.x.minus.xbar(g.u.rep.matrix)
summary(mad.g.u.rep.ma)

t.test(mad.g.u.rep.ma, mad.g.o.rep.ma)
var.test(mad.g.u.rep.ma, mad.g.o.rep.ma)

p.vector<-apply(sig.g.u.rep.matrix,2,function(x)
  shapiro.test(residuals(lm(x~sig.g.o.rep.matrix[,2])))$p.value)

head(p.vector);length(p.vector);mean(p.vector)
sd(p.vector);IQR(p.vector); mad(p.vector)
boxplot.stats(p.vector, coef=.5, do.conf=F, do.out=F)

boxplot.default(p.vector,range=1.4,notch=T,
                col=c("blue","green","black"),
                horizontal=T,
                main="P-values from Shapiro-Wilk Normality Test" )

hist(p.vector, ylab="Frequency",
     xlab="P-Values",
     main="P-Values from Shapiro-Wilk normality test")

grubbs.test(p.vector);sum(p.vector <= 0.05);sum(p.vector > 0.06)
length(p.vector);ad.test(p.vector); summary(p.vector)


t.test(x = Rh5.1, y = Rh5.2,alternative = "two.sided")
wilcox.test(Rh5.1,Rh5.2, exact = FALSE, correct = FALSE) #large sample
var.test(Rh5.1,Rh5.2)
chisq.test(Rh5.1,Rh5.2)
ks.test(Rh5.1,Rh5.2)  # perform ks test
cor.test(Rh5.1,Rh5.2)

Rh6.1<-Rh6samples.gene_exp[,"R8_Rh6_d1_rep1"]
Rh6.2<-Rh6samples.gene_exp[,"R8_Rh6_d1_rep2"]

t.test(x = Rh6.1, y = Rh6.2,alternative = "two.sided")
wilcox.test(Rh6.1,Rh6.2, exact = FALSE, correct = FALSE) #large sample
var.test(Rh6.1,Rh6.2)
chisq.test(Rh6.1,Rh6.2)
ks.test(Rh6.1,Rh6.2)  # perform ks test
cor.test(Rh6.1,Rh6.2)
###################################################
### code chunk number 4: Descriptive stats
###################################################

Rh5var<-apply(Rh5samples.gene_exp, 1,var)
Rh6var<-apply(Rh6samples.gene_exp, 1,var)
summary(Rh5var)
summary(Rh6var)

t.test(x = Rh5var, y = Rh6var,alternative = "two.sided")
wilcox.test(Rh5var,Rh6var, exact = FALSE, correct = FALSE) #large sample
var.test(Rh5var,Rh6var)
chisq.test(Rh5var,Rh6var)
ks.test(Rh5var,Rh6var)  # perform ks test
cor.test(Rh5var,Rh6var)

gse116969.var<-cbind(Rh5var, Rh6var)
head(gse116969.var)

Rh5sd<-apply(Rh5samples.gene_exp, 1,sd)
Rh6sd<-apply(Rh6samples.gene_exp, 1,sd)
summary(Rh5sd)
summary(Rh6sd)
gse116969.sd<-cbind(Rh5sd, Rh6sd)
head(gse116969.sd)

t.test(x = Rh5sd, y = Rh6sd,alternative = "two.sided")
wilcox.test(Rh5sd,Rh6sd, exact = FALSE, correct = FALSE) #large sample
var.test(Rh5sd,Rh6sd)
chisq.test(Rh5sd,Rh6sd)
ks.test(Rh5sd,Rh6sd)  # perform ks test
cor.test(Rh5sd,Rh6sd)

gse116969.sd.df<-cbind(gse116969, gse116969.sd)
head(gse116969.sd.df)

sd.ratio<-(c(Rh5sd+.01)/c(Rh6sd+.01))
plot(sd.ratio)
summary(sd.ratio)
sd.zero<-sd.ratio[which(sd.ratio < 1)]

Rh5mean.expr<-apply(Rh5samples.gene_exp, 1,mean)
Rh6mean.expr<-apply(Rh6samples.gene_exp, 1,mean)
head(Rh5mean.expr)
head(Rh6mean.expr)

Rh8mean.expr<-cbind(Rh5mean.expr,Rh6mean.expr)
head(Rh8mean.expr)

gse116969.mean.expr.df<-cbind(gse116969, Rh8mean.expr)
head(gse116969.mean.expr.df)

gse116969.expr.df<-cbind(gse116969.mean.expr.df, sd.ratio)
head(gse116969.expr.df)

gse116969.sd.df<-cbind(gse116969, sd.ratio)
head(gse116969.sd.df)

sorted.sd<-order(sd.ratio, decreasing = F)
head(gse116969.expr.df[sorted.sd,])
gse116969.df<-gse116969.expr.df[sorted.sd,]

t.test(x = Rh5mean.expr, y = Rh6mean.expr,alternative = "two.sided")
wilcox.test(Rh5mean.expr,Rh6mean.expr, exact = FALSE, correct = FALSE) #large sample
var.test(Rh5mean.expr,Rh6mean.expr)
chisq.test(Rh5mean.expr,Rh6mean.expr)
ks.test(Rh5mean.expr,Rh6mean.expr)  # perform ks test
cor.test(Rh5mean.expr,Rh6mean.expr)


isos_exp.diff<-diffData(isoforms(cuff))
head(isos_exp.diff);dim(isos_exp.diff)
iso.xloc.matrix<-featureNames(isoforms(cuff.hg19))
head(iso.xloc.matrix); length(iso.xloc.matrix)

iso.list<-getGenes(cuff,geneId = iso.xloc.matrix)
iso.rep.ma<-iso.list@genome
iso.rep.ma<-iso.list@repFpkm
head(iso.rep.ma)

iso.repFPKM.ma<-iso.list@isoforms@annotation
head(iso.repFPKM.ma)
fpkm.per.condition<-as.factor(iso.repFPKM.ma$sample_name)
fpkm.per.sample<-as.factor(iso.repFPKM.ma$rep_name)
over.iso.repFPKM.ma<-subset(iso.repFPKM.ma, fpkm.per.condition == over & iso.repFPKM.ma$fpkm > 0 & iso.repFPKM.ma$status == "OK")
under.iso.repFPKM.ma<-subset(iso.repFPKM.ma, fpkm.per.condition == under & iso.repFPKM.ma$fpkm > 0 & iso.repFPKM.ma$status == "OK")
t.under.iso.repFPKM.ma<-t(under.iso.repFPKM.ma)
head(under.iso.repFPKM.ma)
under.cnt.var<-lm(under.iso.repFPKM.ma$fpkm ~ under.iso.var.ma$variance)
plot(under.iso.repFPKM.ma$fpkm ~ under.iso.var.ma$variance)
plot(under.iso.var.ma[,c(3:7)])
hist(under.iso.var.ma$variance,breaks=10)
iso_annotation_data<-featureNames(iso.list@isoforms)
head(iso_annotation_data); dim(iso_annotation_data)


iso.count.ma<-iso.list@count
head(iso.count.ma)
group.variation<-as.factor(iso.count.ma$sample_name)

under.iso.var.ma<-subset(iso.count.ma, group.variation == under & iso.count.ma$count > 0 & iso.count.ma$status == "OK")
head(under.iso.var.ma)
under.cnt.var<-lm(under.iso.var.ma$count ~ under.iso.var.ma$variance)
plot(under.iso.var.ma$variance ~ under.iso.var.ma$count)
plot(under.iso.var.ma[,c(3:7)])
hist(under.iso.var.ma$variance,breaks=10)
iso_annotation_data<-featureNames(iso.list@isoforms)
head(iso_annotation_data); dim(iso_annotation_data)


under.iso.var.ma<-subset(iso.count.ma, group.variation == over & iso.count.ma$count > 0 & iso.count.ma$status == "OK")
head(under.iso.var.ma)
under.cnt.var<-lm(under.iso.var.ma$count ~ under.iso.var.ma$variance)
plot(under.iso.var.ma$variance ~ under.iso.var.ma$count)
plot(under.iso.var.ma[,c(3:7)])

iso_annotation_data<-featureNames(iso.list@isoforms)
head(iso_annotation_data); dim(iso_annotation_data)
one.iso.genes<-unique(iso_annotation_data$gene_short_name)
iso_annotation_data
# gets significant isoform features
mySigIsos<-getSig(cuff,x=over,y=under,level='isoforms')
length(mySigIsos); head(mySigIsos)
# gets significant genes correspondint to significant isoform features
# sigIsos < mySigIsos
# sigIsos = genes.of.sigIsos
genes.of.sigIsos<-getGenes(cuff, mySigIsos)
gene.sig_isos_annot<-genes.of.sigIsos@promoters
gene.sig_isos_annot;length(gene.sig_isos_annot)

iso_annotation_data<-featureNames(genes.of.sigIsos)
head(iso_annotation_data); dim(iso_annotation_data)

sig_isos_exp.diff<-diffData(isoforms(sigIsos))
head(sig_isos_exp.diff);dim(sig_isos_exp.diff)

sig_cds_exp.diff["gene_id"]<-sig_isos_annot

iso.rep.matrix<-repFpkmMatrix(isoforms(cuff))
head(iso.rep.matrix);dim(iso.rep.matrix)


iso.rep.matrix<-as.data.frame(repFpkmMatrix(isoforms(cuff)))
row.names(iso.rep.matrix)->tracking_id

gene.rep.ma<-cbind(tracking_id, gene.rep.matrix)
head(gene.rep.ma)

sig_cds_exp.diff["gene_id"]<-gene_annotation_data["gene_short_name"]

gene.features<-annotation(isoforms(cuff))
head(gene.features)
sig_cds_exp.diff["gene_id"]<-sig_cds_annot["gene_short_name"]
sig_gene_data  <- subset(genes_exp.diff, (log2_fold_change > 5 & q_value <= 0.01))
head(sig_gene_data);dim(sig_gene_data)
sig_geneids <- c(sig_gene_data$gene_id)

###########################################################################
### Chunk_2.b: Significantly differentially expressed CDS
############################################################################

mySigCDS<-getSig(cuff,x=over,y=under,alpha=.05,level='CDS')
head(mySigCDS)
length(mySigCDS)
mySigCDSTable<-getSigTable(cuff,alpha=0.05,level='CDS')
head(mySigCDSTable,20)
length(mySigCDSTable)
sigCDS<-getGenes(cuff, mySigCDS)
sigCDS
length(sigCDS)
sig_cds_exp.diff<-diffData(sigCDS)
sig_cds_annot<-sigCDS@annotation
head(sig_cds_exp.diff)
head(sig_cds_annot)
sig_cds_exp.diff["gene_id"]<-sig_cds_annot["gene_short_name"]
head(sig_cds_exp.diff)

###########################################################################
### Chunk_2.b: Significantly differentially expressed TSS
############################################################################

mySigTSS<-getSig(cuff,x=over,y=under,alpha=.05,level='TSS')
head(mySigTSS)
length(mySigTSS)
mySigTSSTable<-getSigTable(cuff,alpha=0.05,level='TSS')
head(mySigTSSTable,20)
length(mySigTSSTable)
sigTSS<-getGenes(cuff, mySigTSS)
sigTSS
length(sigTSS)
sig_tss_exp.diff<-diffData(sigTSS)
sig_tss_annot<-sigTSS@annotation
head(sig_tss_exp.diff)
head(sig_tss_annot)
sig_tss_exp.diff["gene_id"]<-sig_tss_annot["gene_short_name"]
head(sig_tss_exp.diff)
dim(sig_tss_exp.diff)

###################################################
### Chunk_3: Differentially Loaded Features
###################################################

mySigrelCDS<-getSig(cuff,x=over,y=under,alpha=.05,level='relCDS')
head(mySigrelCDS)
length(mySigrelCDS)
sigrelCDS<-getGenes(cuff, mySigrelCDS)
sigrelCDS
length(sigrelCDS)
relcds.diff<-distValues(relCDS(cuff))
relcds.sigdiff<-subset(relcds.diff, significant=="yes")
head(relcds.sigdiff)
dim(relcds.sigdiff)

mySigSplices<-getSig(cuff,x=over,y=under,alpha=.05,level='splicing')
head(mySigSplices)
length(mySigSplices)
sigSplices<-getGenes(cuff, mySigSplices)
sigSplices
length(sigSplices)
splice.diff<-distValues(splicing(cuff))
splice.sigdiff<-subset(splice.diff, significant=="yes")
head(splice.sigdiff)
dim(splice.sigdiff)

mySigPromoters<-getSig(cuff,x=over,y=under,alpha=.05,level='promoters')
head(mySigPromoters)
length(mySigPromoters)
sigPromoters<-getGenes(cuff, mySigPromoters)
sigPromoters
length(sigPromoters)
promoter.diff<-distValues(promoters(cuff))
promoter.sigdiff<-subset(promoter.diff, significant=="yes")
head(promoter.sigdiff)
dim(promoter.sigdiff)

###################################################
### Chunk_4: Sig Matrix_plots
###################################################

mySigMat<-sigMatrix(cuff,level='genes',alpha=0.05)
print(mySigMat)
sigIsoformIds1<-sigMatrix(cuff, level='isoforms',alpha=0.05)
print(sigIsoformIds1)
sigIsoformIds2<-sigMatrix(cuff, level='TSS',alpha=0.05)
print(sigIsoformIds2)
sigIsoformIds3<-sigMatrix(cuff, level='CDS', alpha=0.05)
print(sigIsoformIds3)

###############################################################################
# Chunk_5: Features, Counts, and all inclusive tables
# (CDS-p_ids, TSS-TSS ids, isoforms- XLOC and TCONS ids, and Genes XLOC ids)
################################################################################
data(geneList)
de <- names(geneList)[abs(geneList) > 1.5]
head(de)
myGenes <- getGenes(cuff, sig_geneids)
myGenes

plot(expressionBarplot(myGenes,logMode=T,pseudocount=1))

## [1] "4312"  "8318"  "10874" "55143" "55388" "991"

x <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
head(as.data.frame(x))


runInfo(cuff)
reps<-replicates(cuff)
replicates<-reps$rep_name
conditions(cuff)

gene.counts<-count(cuff)
head(gene.counts)
gene.featurenames<-featureNames(genes(cuff))
head(gene.featurenames)

iso.counts<-count(isoforms(cuff))
head(iso.counts)
iso.featurenames<-featureNames(isoforms(cuff))
head(iso.featurenames)

cds.counts<-count(CDS(cuff))
head(cds.counts)
cds.featurenames<-featureNames(CDS(cuff))
head(cds.featurenames)

tss.counts<-count(TSS(cuff))
head(tss.counts)
tss.featurenames<-featureNames(TSS(cuff))
head(tss.featurenames)

###################################################################
### Chunk_7: get gene symbols for CDS TSS Isoforms and Genes
#####################################################################

gene.xloc.matrix<-featureNames(genes(cuff))
head(gene.xloc.matrix)
gene.list<-getGenes(cuff,geneId = gene.xloc.matrix)
gene.list
gene_annotation_data<-featureNames(gene.list)
head(gene_annotation_data)

isos.tcons.matrix<-featureNames(isoforms(cuff))
head(isos.tcons.matrix)
iso.list<-getGenes(cuff,geneId = isos.tcons.matrix)
iso.list
iso_annotation_data<-featureNames(iso.list)
head(iso_annotation_data)

cds.p_id.matrix<-featureNames(CDS(cuff))
head(cds.p_id.matrix)
cds.list<-getGenes(cuff,geneId = cds.p_id.matrix)
cds.list
cds_annotation_data<-featureNames(cds.list)
head(cds_annotation_data)

tss.tssgroup.matrix<-featureNames(TSS(cuff))
head(tss.tssgroup.matrix)
tss.list<-getGenes(cuff,geneId = tss.matrix)
tss.list
tss_annotation_data<-featureNames(tss.list)
head(tss_annotation_data)

####################################################################
### code chunk number 8: Write repFPKMMatrix and DiffTable to file
####################################################################

g.rep.matrix<-repFpkmMatrix(genes(cuff))
head(g.rep.matrix)
gene.xloc.matrix<-featureNames(genes(cuff))
gene.list<-getGenes(cuff,geneId = gene.xloc.matrix)
gene.list

gene_annotation_data<-featureNames(gene.list)
head(iso_annotation_data)

gene_rep_matrix<-cbind2(gene_annotation_data, g.rep.matrix)
gene.rep.matrix<-gene_rep_matrix[,-1]
fpkm.file = file.path("FPKMmatrix")
write.table(gene.rep.matrix,file = fpkm.file, sep = "  ", row.names = F, col.names = T,quote = F)

# This table are huge and chew up a lot of memory
isodiff<-diffTable(isoforms(cuff))
diff.table = file.path("DiffTable.dataframe")
write.table(isodiff, file = diff.table, sep = "  ", row.names = F , col.names = T,quote = F)


##################################################
### code chunk number 10: global_dispersion
###################################################

disp_cuff<-dispersionPlot(cuff) + ggtitle("Overloading of All Cuffdiff Output")
disp_cuff
disp_gene<-dispersionPlot(genes(cuff)) + ggtitle("CUFFDIFF Genes Overloading")
disp_gene

###################################################
### code chunk number 11: SCV_visualization
###################################################
genes.scv<-fpkmSCVPlot(genes(cuff))
genes.scv
genes.scv_pool<-fpkmSCVPlot(genes(cuff), showPool=T)
genes.scv_pool
genes.scv<-fpkmSCVPlot(genes(cuff), showPool=F)
genes.scv

isoforms.scv<-fpkmSCVPlot(isoforms(cuff))
isoforms.scv
isoforms.scv_pool<-fpkmSCVPlot(isoforms(cuff), showPool=T)
isoforms.scv_pool
isoforms.scv<-fpkmSCVPlot(isoforms(cuff), showPool=F)
isoforms.scv

cvtss<-fpkmSCVPlot(TSS(cuff))
cvtss
cvtss_pool<-fpkmSCVPlot(TSS(cuff), showPool=T)
cvtss_pool
tss.scv<-fpkmSCVPlot(TSS(cuff), showPool=F)
tss.scv


cvcds1<-fpkmSCVPlot(CDS(cuff))
cvcds1
cvcds_pool<-fpkmSCVPlot(CDS(cuff), showPool=T,)
cvcds_pool
cvcds<-fpkmSCVPlot(CDS(cuff), showPool=F)
cvcds

###################################################
### code chunk number 12: global_plots_1
###################################################
cuffgenes<-genes(cuff)
cuffgene.table<-cuff@tables
cuffgene.table$countTable

dens1<-csDensity(genes(cuff),replicates=F,size=1)
print(dens1)

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
box_g<-csBoxplot(genes(cuff),replicates=T, logMode=T) + ggtitle("genes")
box_g
box_i<-csBoxplot(isoforms(cuff),replicates=T, logMode=T) + ggtitle("isoforms")
box_i
box_c<-csBoxplot(CDS(cuff),replicates=T, logMode=T) + ggtitle("CDS")
box_c
box_t<-csBoxplot(TSS(cuff),replicates=T, logMode=T) + ggtitle("TSS")
box_t
###################################################
### code chunk number 14: global_plots_scatter_1
###################################################

gene_scatter<-csScatter(genes(cuff),over,under,smooth=T) + ggtitle("genes")
gene_scatter
iso_scatter<-csScatter(isoforms(cuff), over, under,smooth=T, logMode=T) + ggtitle("isoforms")
iso_scatter
tss_scatter<-csScatter(TSS(cuff), over, under,smooth=T) + ggtitle("TSS")
tss_scatter
cds_scatter<-csScatter(CDS(cuff), over, under,smooth=T) + ggtitle("CDS")
cds_scatter
scgene<-csScatter(sigGenes,over,under,smooth=T) + ggtitle("sig Genes")
scgene
sciso<-csScatter(sigIsos,over, under,smooth=T) + ggtitle("Sig Isoforms")
sciso
sctss<-csScatter(sigTSS,over, under,smooth=F) + ggtitle("Sig TSS Groups")
sctss
sctss<-csScatter(sigCDS,over, under,smooth=F) + ggtitle("Sig CDS")
sctss
##################################################
### code chunk number 15: global_plots_dendro
###################################################

csDendro(genes(cuff), replicates=T, logMode=F)
csDendro(isoforms(cuff),replicates=T, logMode=F)
csDendro(TSS(cuff),replicates=T, logMode=F)
csDendro(CDS(cuff),replicates=T, logMode=F)

plot(csDendro(sigGenes,replicates=T, logMode=F,), main="Sig Genes")
plot(csDendro(sigTSS,replicates=T, logMode=F,), main="Sig TSS")
plot(csDendro(sigCDS,replicates=T, logMode=F,), main="Sig CDS")
plot(csDendro(sigIsos,replicates=T, logMode=F,), main="Sig Isoforms")

###################################################
### code chunk number 16: global_plots_4
###################################################

ma_g<-MAplot(genes(cuff),over,under) + ggtitle("genes")
ma_g
ma_i<-MAplot(isoforms(cuff),over,under,useCount=T) + ggtitle("isoforms")
ma_i
ma_c<-MAplot(CDS(cuff),over,under,useCount=T) + ggtitle("CDS")
ma_c
ma_t<-MAplot(TSS(cuff),over, under,useCount=T) + ggtitle("TSS")
ma_t

###################################################
### Code Chunk 17:  volcanos
###################################################

vv<-csVolcano(genes(cuff),over, under,) + ggtitle("genes")
vv
vy<-csVolcano(TSS(cuff),over, under) + ggtitle("TSS")
vy
vu<-csVolcano(CDS(cuff),over, under) + ggtitle("CDS")
vu
vu<-csVolcano(isoforms(cuff),over, under) + ggtitle("Isoforms")
vu

vol_g<-csVolcanoMatrix(genes(cuff)) + ggtitle("genes")
vol_g
vol_i<-csVolcanoMatrix(isoforms(cuff)) + ggtitle("isoforms")
vol_i
vol_c<-csVolcanoMatrix(CDS(cuff)) + ggtitle("CDS")
vol_c
vol_t<-csVolcanoMatrix(TSS(cuff)) + ggtitle("TSS")
vol_t

vol_sigT<-csVolcano(sigIsos,over,under, features=T) + ggtitle("Sig Isoforms")
vol_sigT
vol_sigC<-csVolcano(sigCDS,over,under, features=T) + ggtitle("Sig CDS")
vol_sigC
vol_sigI<-csVolcano(sigTSS,over,under, features=T) + ggtitle("Sig TSS_group")
vol_sigI
vol_sigG<-csVolcano(sigGenes,over,under, features=T) + ggtitle("Sig Genes")
vol_sigG
vol_sigmyG<-csVolcano(myGenes,over,under, features=T) + ggtitle("Genes of Interest")
vol_sigmyG

##############################################################
### code chunk number 18: create_geneset/genes_of_interest
###############################################################

gene.counts<-count(genes(cuff))
head(gene.counts)

gene.featurenames<-featureNames(genes(cuff))
head(gene.featurenames)

iso.counts<-count(isoforms(cuff))
head(iso.counts)

iso.featurenames<-featureNames(isoforms(cuff))
head(iso.featurenames)

cds.counts<-count(CDS(cuff))
head(cds.counts)

cds.featurenames<-featureNames(CDS(cuff))
head(cds.featurenames)

tss.counts<-count(TSS(cuff))
head(tss.counts)

tss.featurenames<-featureNames(TSS(cuff))
head(tss.featurenames)
###################################################

# the 'geneListString' variable is set in the settings file
# The environment variables are initialized at the begining of this script
geneListString=c("CXCL12","CXCR4","TGFB1","TGFB","IL6","IL13","NFKB")
myGeneIds = unlist(strsplit(geneListString, " "))
myGeneIds

myGeneSet<-getGenes(cuff,myGeneIds)
myGeneSet<-getGenes(cuff,myGeneIds, sampleIdList=c(over,under))
myGeneSet

myGenes<-getGenes(cuff,geneIdList=myGeneIds)
myGenes

gene.counts<-count(genes(myGenes))
head(gene.counts)
gene.featurenames<-featureNames(genes(myGenes))
head(gene.featurenames)

GeneSet.counts<-count(genes(myGeneSet))
head(GeneSet.counts)
GeneSet.featurenames<-featureNames(genes(myGeneSet))
head(GeneSet.featurenames)

###################################################
### code chunk number 19: geneset_plots
###################################################

exP.iso<-expressionPlot(myGenes, replicates=T,logMode=T,showErrorbars =T) + ggtitle("Genes of interest")
exP.iso
exp.cds2<-expressionPlot(CDS(myGenes), replicates=T, logMode=T) + ggtitle("CDS of Genes of Interest")
exp.cds2
exp.tss2<-expressionPlot(TSS(myGenes), replicates=T,logMode=T,showErrorbars =T) + ggtitle("TSS of Genes of Interest")
exp.tss2
exp.iso2<-expressionPlot(isoforms(myGenes), replicates=T, logMode=T) + ggtitle("Isoforms of Genes of Interest")
exp.iso2

#####################################################
### code chunk number 20: genes_of_interest_heatmaps
#####################################################

hmiso.rep<-csHeatmap(sigIsos,cluster='both',replicates=T,labRow=F,labCol=T, logMode = T ) + ggtitle("Sig Isoforms")
print(hmiso.rep)
hmiso<-csHeatmap(sigIsos,cluster='both',replicates=F,labRow=F,labCol=T, logMode = T ) + ggtitle("Sig Isoforms")
print(hmiso)

hmtss.rep<-csHeatmap(sigTSS,cluster='both',replicates=T,labRow = F,labCol=T, logMode = T) + ggtitle("Sig TSS_Group")
print(hmtss.rep)
hmtss<-csHeatmap(sigTSS, cluster="row",labRow=F, replicates=F) + ggtitle("sigTSS")
hmtss

hmcds.rep<-csHeatmap(sigCDS,cluster='both',replicates=T,labRow = F,labCol=T, logMode = T) + ggtitle("Sig CDS")
print(hmcds.rep)
hmcds<-csHeatmap(sigCDS,cluster='both',replicates=F,labRow = F,labCol=T, logMode = T) + ggtitle("Sig CDS")
print(hmcds)

hmsigG<-csHeatmap(sigGenes,cluster='both', labRow=F,replicates=F) + ggtitle("Sig Genes")
hmsigG
h.repsigG<-csHeatmap(sigGenes,cluster='both',replicates=T, labRow=F) + ggtitle("Sig Genes")
h.repsigG

#h.cs1<-csHeatmap(isoforms(sigGenes), cluster='both', labRow=F, replicates=T) + ggtitle("isoforms(sigGenes)")
#h.cs1
#h.cs2<-csHeatmap(TSS(sigGenes), cluster='both', replicates=T) + ggtitle("TSS(sigGenes)")
#h.cs2
#h.cs3<-csHeatmap(CDS(sigGenes), cluster='both',labRow=F, replicates=T) + ggtitle("CDS(sigGenes)")
#h.cs3

hg<-csHeatmap(myGenes,cluster="both", replicates=T, labRow=T) + ggtitle("Genes of Interest")
hg
ihiso<-csHeatmap(isoforms(myGenes),cluster='both',labRow=T) + ggtitle("Isoforms of interest")
ihiso
thtss<-csHeatmap(TSS(myGenes),cluster='both',labRow=T) + ggtitle("TSS_Groups of Interest")
thtss
ihcds<-csHeatmap(CDS(myGenes),cluster='both',labRow=T) + ggtitle("CDS_Groups of Interest")
ihcds

hg<-csHeatmap(genes(cuff),cluster="both", replicates=T, labRow=T) + ggtitle("All Genes Features")
hg
ihiso<-csHeatmap(isoforms(cuff),cluster='both',labRow=T) + ggtitle("All Isoforms Features")
ihiso
thtss<-csHeatmap(TSS(cuff),cluster='both',labRow=T) + ggtitle("All TSS_Groups Features")
thtss
ihcds<-csHeatmap(CDS(cuff),cluster='both',labRow=T) + ggtitle("All CDS Features")
ihcds

##################################################
### code chunk number 21: geneset_plots_1.5
###################################################

expbar_i<-expressionBarplot(sigIsos, showErrorbars = T, logMode = T) + ggtitle("Significantly Differentially Expressed Isofroms")
expbar_i
expbar_t<-expressionBarplot(sigTSS, showErrorbars = T, logMode = T) + ggtitle("Significantly Differentially Expressed TSS_groups")
expbar_t
expbar_c<-expressionBarplot(sigCDS, showErrorbars = T, logMode = T) + ggtitle("Significantly Differentially Expressed CDS")
expbar_c
expbar_g3<-expressionBarplot(sigGenes, showErrorbars = T, logMode = T) + ggtitle("Significantly Differentially Expressed Genes")
expbar_g3
expbar_myg<-expressionBarplot(myGenes) + ggtitle("Genes of Interest")
expbar_myg

####################################
### Code Chunk 22-23:  Gene Loop
####################################
myGeneList=c("CXCL12","CXCR4","IL4","TNF")
myGeneIds<-unlist(myGeneList)
genesArray = c(myGeneIds)
for (i in 1:length(genesArray)) {
	print(genesArray[i])
    myGeneId = genesArray[i]
    myGeneId
    #myGene<-getGene(cuff,genesArray[i])
    myGene<-getGene(cuff,myGeneId)
    myGene
    gene.feats<-features(gene.sig_isos_annot)
    genefeats<-featureNames(gene.sig_isos_annot)
    names(myGene)}
  #  myGene[is.na(myGene)]<-c(0)
##########################################################
### Code Chunk 24:  ExpressionPlots/charts - Bar and Pie
###########################################################
exp.myg<-expressionPlot(myGene, replicates=FALSE)
print(exp.myg)
exP.tss.rep<-expressionPlot(isoforms(myGene), replicates=T,logMode=T, drawSummary = F) + ggtitle("Isoform of Genes of Interest")
print(exP.tss.rep)
exP.tss<-expressionPlot(TSS(myGene), replicates=T,logMode=T, showErrorbars=T, drawSummary=T) + ggtitle("TSS of Genes of Interest")
print(exP.tss)
exP.cds<-expressionPlot(CDS(myGene), replicates=T,logMode=T,showErrorbars =T, drawSummary = T) + ggtitle("CDS of Genes of Interest")
print(exP.cds)
##########################################################
### Code Chunk 24.5 Pie
###########################################################

gp1<-csPie(myGene,fpkm,level="CDS",replicates=TRUE)
print(gp1)
gp2<-csPie(myGene,level="isoforms",replicates=TRUE)
print(gp2)
gp3<-csPie(myGene,fpkm,level="TSS",replicates=TRUE)
print(gp3)
gp4<-csPie(myGene,level="genes",replicates=TRUE)
print(gp4)
??makeGeneRegionTrack
###################################################
### code chunk number 25: features_2
###################################################
 BiomartGeneRegionTrack(start=myGene@features)
genetrack<-makeGeneRegionTrack(object=myGene@features)
 plotTracks(genetrack)

 trackList<-list()
 myStart<-min(features(myGene)$start)
 myEnd<-max(features(myGene)$end)
 myChr<-unique(features(myGene)$seqnames)

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
### code chunk number 26: dist_heat_maps
###################################################

DistHeatG<-csDistHeat(genes(cuff),replicates=F) + ggtitle("Genes (w/o replicates)")
DistHeatG
RepDistHeatG<-csDistHeat(genes(cuff),replicates=T) + ggtitle("Genes (w/ replicates)")
RepDistHeatG

DistHeatI<-csDistHeat(isoforms(cuff),replicates=F) + ggtitle("genes (w/o replicates)")
DistHeatI
RepDistHeatI<-csDistHeat(isoforms(cuff),replicates=T) + ggtitle("Genes (w/ replicates)")
RepDistHeatI

DistHeatT<-csDistHeat(TSS(cuff),replicates=FALSE) + ggtitle("CDS, (w/o replicates)")
DistHeatT
RepDistHeatT<-csDistHeat(TSS(cuff),replicates=T) + ggtitle("TSS, (w/ replicates)")
RepDistHeatT

DistHeatC<-csDistHeat(CDS(cuff),replicates=FALSE) + ggtitle("CDS (w/o replicates)")
DistHeatC
RepDistHeatC<-csDistHeat(CDS(cuff),replicates=T) + ggtitle("CDS (w/o replicates)")
RepDistHeatC

Heat_rel<-csDistHeat(relCDS(cuff)) + ggtitle("relCDS")
Heat_rel
Heat_promoters<-csDistHeat(promoters(cuff)) + ggtitle("Promoters")
Heat_promoters
Heat_splicing<-csDistHeat(splicing(cuff)) + ggtitle("Splicing")
Heat_splicing

sigDistHeatG<-csDistHeat(sigGenes) + ggtitle("Genes-sigGenes")
sigDistHeatG
sigRepDistHeatG<-csDistHeat(sigGenes, replicates=T) + ggtitle("Genes-sigGenes")
sigRepDistHeatG
sigHeat_iso_rep<-csDistHeat(isoforms(sigGenes), replicates=T) + ggtitle("RelCDS-sigGenes")
sigHeat_iso_rep
sigHeat_tss_rep<-csDistHeat(TSS(sigGenes), replicates=T) + ggtitle("Promoters-sigGenes")
sigHeat_tss_rep
sigHeat_cds_rep<-csDistHeat(CDS(sigGenes), replicates=T) + ggtitle("Splicing-sigGenes")
sigHeat_cds_rep

sigHeat_rel<-csDistHeat(relCDS(sigGenes)) + ggtitle("RelCDS-sigGenes")
sigHeat_rel
sigHeat_promoters<-csDistHeat(promoters(sigGenes)) + ggtitle("Promoters-sigGenes")
sigHeat_promoters
sigHeat_splicing<-csDistHeat(splicing(sigGenes)) + ggtitle("Splicing-sigGenes")
sigHeat_splicing
###################################################
### code chunk number 28: dim_reduction_1 PCA
###################################################
genes.PCA<-PCAplot(genes(cuff),"PC1","PC2",scale=T, replicates=FALSE) + ggtitle("Genes")
genes.PCA
rep_genes.PCA<-PCAplot(genes(cuff),"PC1","PC2",scale=T, replicates=TRUE) + ggtitle("Genes, replicates=T")
rep_genes.PCA

cds.PCA<-PCAplot(CDS(cuff),"PC1","PC2",scale =T, replicates=F) + ggtitle("CDS")
cds.PCA
rep_cds.PCA<-PCAplot(CDS(cuff),"PC1","PC2",scale =T, replicates=TRUE) + ggtitle("CDS, replicates=T")
rep_cds.PCA

tss.PCA<-PCAplot(TSS(cuff),"PC1","PC2",replicates=F,scale =T) + ggtitle("TSS")
tss.PCA
rep_tss.PCA<-PCAplot(TSS(cuff),"PC1","PC2",replicates=TRUE,scale =T) + ggtitle("TSS, replicates=T")
rep_tss.PCA

isos.PCA<-PCAplot(isoforms(cuff),"PC1","PC2",replicates=F,scale =T) + ggtitle("Isoforms")
isos.PCA
rep_isos.PCA<-PCAplot(isoforms(cuff),"PC1","PC2", replicates=T, scale=T, pseudocount = 1) + ggtitle("Isoforms, replicates=T")
rep_isos.PCA

###################################################
### code chunk number 29: dim_reduction_2 MDS
###################################################

genes.MDS.rep<-MDSplot(genes(cuff), replicates=T, logMode=T) + ggtitle("Genes, replicates=T")
print(genes.MDS.rep)
genes.MDS.iso<-MDSplot(isoforms(cuff),logMode=T, replicates=TRUE) + ggtitle("Isoforms, replicates=T")
genes.MDS.iso
genes.MDS.cds<-MDSplot(CDS(cuff),logMode=T, replicates=TRUE) + ggtitle("CDS, replicates=T")
genes.MDS.cds
genes.MDS.tss<-MDSplot(TSS(cuff),logMode=T, replicates=TRUE) + ggtitle("TSS, replicates=T")
genes.MDS.tss

siggenes.MDS<-MDSplot(sigGenes,replicates=TRUE, logMode=T) + ggtitle("SIGGENES, replicates=T")
siggenes.MDS
sig_genes.MDS.iso<-MDSplot(isoforms(sigGenes),logMode=T, replicates=TRUE) + ggtitle("SIG_ISOFORMS, replicates=T")
sig_genes.MDS.iso
sig_genes.MDS.cds<-MDSplot(CDS(sigGenes),logMode=T, replicates=TRUE) + ggtitle("SIG_CDS, replicates=T")
sig_genes.MDS.cds
sig_genes.MDS.tss<-MDSplot(TSS(sigGenes),logMode=T, replicates=TRUE) + ggtitle("SIG_TSS, replicates=T")
sig_genes.MDS.tss


#################################################
### code chunk number 29: geneset_cluster_1
###################################################

 ic<-csCluster(genes(myGenes),k=4)
 head(ic$cluster)
 icp<-csClusterPlot(ic)
 icp

###################################################
### code chunk number 30: specificity_1
###################################################

 myGenes.spec<-csSpecificity(genes(myGenes))
 head(myGenes.spec)

#####################################################
##### code chunk number 31: similar_1
#####################################################

 mySimilar<-findSimilar(genes(cuff),myGene,n=20)
 mySimilar.expression<-expressionPlot(mySimilar,logMode=T,showErrorbars=F)

######################################################
##### code chunk number 32: similar_plots_1
#####################################################

mySimilar<-findSimilar(cuff,diffgene,n=20)
mySimilar.expression<-expressionPlot(mySimilar,logMode=T,showErrorbars=F)
print(mySimilar.expression)

###################################################
### code chunk number 33: similar_2
###################################################
 myProfile<-c(500,0,400)
 mySimilar2<-findSimilar(CDS(cuff),myProfile,n=10)
 mySimilar2.expression<-expressionPlot(mySimilar2,logMode=T,showErrorbars=F)
 mySimilar2.expression
 print(mySimilar2.expression)
###################################################
### code chunk number 34: similar_plots_2
###################################################
 mySimilar2.expression<-expressionPlot(mySimilar2,logMode=F,showErrorbars=F)
 mySimilar2.expression
 print(mySimilar2.expression)

###################################################
#gene.fpkm<-fpkm(genes(cuff))
#head(gene.fpkm)
#gene.repFpkm<-repFpkm(genes(cuff))
#head(gene.repFpkm)
#gene.matrix<-fpkmMatrix(genes(cuff))
#head(gene.matrix)
#head(fpkmMatrix(isoforms(cuff))); #head(fpkmMatrix(CDS(cuff)))
#head(fpkmMatrix(TSS(cuff))); #head(fpkmMatrix(sigGenes))
#head(fpkmMatrix(sigIsos)); #head(fpkmMatrix(sigCDS))
#head(fpkmMatrix(sigTSS)); #head(repFpkmMatrix(isoforms(cuff)))
#head(repFpkmMatrix(TSS(cuff))); #head(repFpkmMatrix(CDS(cuff)))
#head(repFpkmMatrix(sigGenes)); #head(repFpkmMatrix(CDS(sigGenes)))
#head(repFpkmMatrix(isoforms(sigGenes))); #head(repFpkmMatrix(TSS(sigGenes)))

###################################################
### code chunk number 9: countMatrix
###################################################
#head(countMatrix(genes(cuff))); #head(countMatrix(TSS(sigGenes)))
#head(countMatrix(isoforms(sigGenes))); #head(countMatrix(CDS(sigGenes)))
#head(countMatrix(sigGenes)); #head(countMatrix(sigIsos))
#head(countMatrix(sigTSS)); #head(countMatrix(sigCDS))
dev.off()

end<-dbDisconnect(cuff@DB)
dev.off(file.path(Rplots))
