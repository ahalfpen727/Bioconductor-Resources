#drosophila analysis of Rh8_Rh5 [2 samples] against Rh8_Rh6 [2 samples]
#######################################################################
### code chunk number 1: Library Load
#######################################################################
library(igraph);library(GEOquery)
library(Hmisc);library(pvclust)
library(outliers);library(gplots)
library(ellipse);library(limma)
library(readxl);library(edgeR)
library(lmtest);library(multtest)
library(ape);library(nortest)
#######################################################################
### code chunk number 2: Load the count matrix
#######################################################################
files<-dir()
files
gse116969<-grep(pattern="Rister_",files, ignore.case = T)
gse116969.xlsx<-files[gse116969]
gse116969.xlsx<-read_excel(gse116969.xlsx)
gse116969.xlsx<-as.data.frame(gse116969.xlsx)
colnames(gse116969.xlsx)

transcripts<-unique(gse116969.xlsx[,"transcript_id"])
row.names(gse116969.xlsx)=transcripts
head(gse116969.xlsx)

gse116969.annotation<-as.data.frame(gse116969.xlsx[,c(1:3)])
gse116969<-as.matrix(gse116969.xlsx[,-c(1:3)])
#gse116969<-as.matrix(gse116969.xlsx)
samps<-colnames(gse116969)
genes<-row.names(gse116969)
head(gse116969)
dim(gse116969)
cg.genes<-grep(pattern="CG",x=gse116969.annotation$gene_name, ignore.case=F)
unknown.cg.genes<-gse116969.annotation[cg.genes,]
unknowncg.genes<-gse116969.xlsx[cg.genes,]
head(unknowncg.genes);dim(unknowncg.genes)
gse116969.df<-as.data.frame(unknowncg.genes[,-c(1:3)], row.names=unknowncg.genes$gene_name)
head(gse116969.df)

isexpr <- rowSums(cpm(gse116969.df) > 1) >= 2
gse116969.df<-gse116969.df[isexpr,]
head(gse116969.df);dim(gse116969.df)
gse116969<-gse116969.df

Rh6samps<-grep(pattern="Rh6", colnames(gse116969), ignore.case = T)
Rh5samps<-grep(pattern="Rh5", colnames(gse116969), ignore.case = T)
head(gse116969[,Rh6samps])
head(gse116969[,Rh5samps])

Rh8_Rh6s<-factor(grep(pattern="Rh6", colnames(gse116969), ignore.case = T))
Rh8_Rh5s<-factor(grep(pattern="Rh5", colnames(gse116969), ignore.case = T))
Rh8_Rh6group<-gse116969[,Rh8_Rh6s]
head(Rh8_Rh6group)
Rh8_Rh5group<-gse116969[,Rh8_Rh5s]
head(Rh8_Rh5group)

##############################################################
### code chunk number 3: Design Matrix and pair like-samples
##############################################################

gseDesign = data.frame(row.names = colnames(gse116969),
                       condition = c("Rh5", "Rh5", "Rh6","Rh6"))
conditions<-factor(gseDesign$condition)

Rh5group = gseDesign$condition == "Rh5"
Rh6group = gseDesign$condition == "Rh6"

Rh5samples.gene_exp = gse116969[,Rh5group]
Rh6samples.gene_exp = gse116969[,Rh6group]

Rh5.1<-Rh5samples.gene_exp[,"R8_Rh5_d1_rep1"]
Rh5.2<-Rh5samples.gene_exp[,"R8_Rh5_d1_rep2"]

t.test(x = Rh5.1, y = Rh5.2,alternative = "two.sided")
wilcox.test(Rh5.1,Rh5.2, exact = FALSE, correct = FALSE) #large sample
var.test(Rh5.1,Rh5.2)
cor.test(Rh5.1,Rh5.2)

Rh6.1<-Rh6samples.gene_exp[,"R8_Rh6_d1_rep1"]
Rh6.2<-Rh6samples.gene_exp[,"R8_Rh6_d1_rep2"]

t.test(x = Rh6.1, y = Rh6.2,alternative = "two.sided")
wilcox.test(Rh6.1,Rh6.2, exact = FALSE, correct = FALSE) #large sample
var.test(Rh6.1,Rh6.2)
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
head(gse116969)

Rh6.1.sd<-sd(Rh6.1)
Rh6.2.sd<-sd(Rh6.2)
Rh5.1.sd<-sd(Rh5.1)
Rh5.2.sd<-sd(Rh5.2)
Rh6.1.sd;Rh6.2.sd;Rh5.1.sd;Rh5.2.sd

Rh6.1.IQR<-range(Rh6.1)
Rh6.2.IQR<-IQR(Rh6.2)
Rh5.1.IQR<-IQR(Rh5.1)
Rh5.2.IQR<-IQR(Rh5.2)
Rh6.1.IQR;Rh6.2.IQR;Rh5.1.IQR;Rh5.2.IQR

Rh6.1.mu<-mean(Rh6.1)
Rh6.2.mu<-mean(Rh6.2)
Rh5.1.mu<-mean(Rh5.1)
Rh5.2.mu<-mean(Rh5.2)
Rh6.1.mu;Rh6.2.mu;Rh5.1.mu;Rh5.2.mu

Rh6var<-apply(Rh6samples.gene_exp, 1,var)
summary(Rh6var)
Rh5var<-apply(Rh5samples.gene_exp, 1,var)
summary(Rh5var)

Rh6varsq<-apply(Rh6samples.gene_exp, 1,sd)
summary(Rh6varsq)
Rh5varsq<-apply(Rh5samples.gene_exp, 1,sd)
plot(Rh5varsq)

gse116969sd<-apply(gse116969, 2,sd)
gse116969mu<-apply(gse116969, 2,mean)
gse116969var<-apply(gse116969, 2,var)

Rh5sd<-apply(Rh5samples.gene_exp, 2,sd)
Rh6sd<-apply(Rh6samples.gene_exp, 2,sd)
summary(Rh5sd)
summary(Rh6sd)

Rh5sd<-apply(Rh5samples.gene_exp, 1,sd)
Rh6sd<-apply(Rh6samples.gene_exp, 1,sd)
summary(Rh5sd)
summary(Rh6sd)
gse116969.sd<-cbind(Rh5sd, Rh6sd)
head(gse116969.sd)

t.test(x = Rh5sd, y = Rh6sd,alternative = "two.sided")
wilcox.test(Rh5sd,Rh6sd, exact = FALSE, correct = FALSE) #large sample
var.test(Rh5sd,Rh6sd)
cor.test(Rh5sd,Rh6sd)

gse116969.sd.df<-cbind(gse116969, gse116969.sd)
head(gse116969.sd.df)

sd.ratio<-(c(Rh5sd+0.01)/c(Rh6sd+0.01))
head(sd.ratio)
plot(sd.ratio)
summary(sd.ratio)
sd.zero<-c(sd.ratio < 1.2 & sd.ratio > 0.8)
head(sd.zero); length(sd.zero)
table(sd.zero)
sd.zero<-c(sd.ratio < 1.1 & sd.ratio > 0.9)
head(sd.zero); length(sd.zero)
table(sd.zero)

Rh5mean.expr<-apply(Rh5samples.gene_exp, 1,mean)
Rh6mean.expr<-apply(Rh6samples.gene_exp, 1,mean)
head(Rh5mean.expr)
head(Rh6mean.expr)

cor.test(Rh6samples.gene_exp$R8_Rh6_d1_rep1, Rh6samples.gene_exp$R8_Rh6_d1_rep2)
cor.test(Rh5samples.gene_exp$R8_Rh5_d1_rep1, Rh5samples.gene_exp$R8_Rh5_d1_rep2)

gse116969.mean.expr.df<-cbind(gse116969, Rh8mean.expr)
head(gse116969.mean.expr.df)

logbase2fold<-log(x=c(Rh5mean.expr/Rh6mean.expr), base=2)

gse116969.mu.expr.df<-cbind(gse116969.mean.expr.df, logbase2fold)
head(gse116969.mu.expr.df)

gse116969.exp.df<-cbind(gse116969.mu.expr.df, gse116969.sd)
head(gse116969.exp.df)

gse116969.exp.sd.df<-cbind(gse116969.exp.df, sd.ratio)
head(gse116969.exp.sd.df)
dim(gse116969.exp.sd.df[sd.zero,])
gse116969.expr.df<-cbind(gse116969, logbase2fold)
gse116969.exp.sd.df<-gse116969.expr.df[sd.zero,]
head(gse116969.exp.sd.df)

gse116969.sd.df<-cbind(gse116969, sd.ratio)
head(gse116969.sd.df)

gse116969.lowsd.df<-gse116969.exp.sd.df[sd.zero,]
head(gse116969.lowsd.df); dim(gse116969.lowsd.df)

gse116969.sorted.df<-gse116969.lowsd.df[order(gse116969.lowsd.df["sd.ratio"], decreasing = F),]
gse116969.sorted.df<-gse116969.lowsd.df[order(gse116969.lowsd.df["logbase2fold"] ,decreasing=T),]
head(gse116969.sorted.df);dim(gse116969.sorted.df)

Rh5.hi.genes<-gse116969.sorted.df[,"logbase2fold"] >= 2
Rh6.hi.genes<-gse116969.sorted.df[,"logbase2fold"] <= -2

Rh6.higenes.df<-gse116969.sorted.df[Rh6.hi.genes,]
head(Rh6.higenes.df); dim(Rh6.higenes.df)

Rh5.higenes.df<-gse116969.sorted.df[Rh5.hi.genes,]
head(Rh5.higenes.df); dim(Rh5.higenes.df)


##################################################
### code chunk number : Step 6: Clustering
##################################################
pdf("Drosophila_comparison_of_Rh8_Rh5_to_Rh8_Rh6")
hist(logbase2fold)
plot(gse116969.sd,main="Standard Deviation of Rh8_Rh5 and Rh8_Rh6 ")
plot(Rh8mean.expr,main="Mean gene expression of Rh8_Rh5 and Rh8_Rh6 ")
heatmap(sig.genes.gse116969)

sig.genes.gse116969.df<-as.data.frame(sig.genes.gse116969)
sig.gse.cor<-cor(Rh5.higenes.df)
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

####################################################################################################
### code chunk number : Step 3: #Hierarchical clustering
####################################################################################################

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


