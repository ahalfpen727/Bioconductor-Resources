######################################################################################################
##Tailor script igraph netowrk analysis for cuffdiff data
## built from Tutorial: Network plot from expression data in R using igraph
## https://www.biostars.org/p/285296/ library(devtools);install_github("genomicsclass/maPooling")

######################################################################################################
### Library+Variable --> 1) Load
######################################################################################################
options(width=65)
#source("http://igraph.sf.net")
library(cummeRbund); library(DOSE);library(clusterProfiler);library(biomaRt)
library(ReactomePA);library(GSEABase); library(GO.db); library(org.Hs.eg.db)
library(KEGG.db);library(reactome.db); library(gage);library(pathview)
library(limma);library(ReactomePA); library(KEGGgraph);library(keggorthology)
library(EGSEAdata);library(gageData); library(msigdf)
library(igraph);library(Hmisc);library(ellipse)
  # These arguments are passed in on the command line via launchR.eachDiffDir.sh
  # args:
  # diffDir=\"${diffDir}\" inDir=\"${INPUT}/${diffDir}\" outDir=\"${OUTPUT}/${diffDir}\"
  # under=\"${under}\" over=\"${over}\" Rplots=\"$\{OUTPUT}/${diffDir}/${R_PLOTS}\"
  # FPKMmatrix=\"${OUTPUT}/${diffDir}/${FPKM_MATRIX}\" DiffTable=\"${OUTPUT}/${diffDir}/${DIFF_TABLE}\"

######################################################################################################
### Library+Variable --> 2) import and print env variables 3) initialize cuff-obj
##############################################################################################
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
# These arguments are passed in on the command line via launchR.eachDiffDir.sh
# args:
# diffDir=\"${diffDir}\" inDir=\"${INPUT}/${diffDir}\" outDir=\"${OUTPUT}/${diffDir}\"  ,
# under=\"${under}\" over=\"${over}\" Rplots=\"$\{OUTPUT}/${diffDir}/${R_PLOTS}\",
# FPKMmatrix=\"${OUTPUT}/${diffDir}/${FPKM_MATRIX}\" DiffTable=\"${OUTPUT}/${diffDir}/${DIFF_TABLE}\"
#pdf(file.path("cummerbund_results_hg38_gtf_guided_Rplots.pdf"))
diffDir<-file.path("./umb_triley/urine1/cuffdiff_results_hg38_default/LUTS-over-CTRL")
setwd(diffDir);dir()
over="LUTS";under="CTRL"
cuff=readCufflinks(dir=".",rebuild=T,verbose = T)
cuff.hg38<-cuff
cuff.hg19=readCufflinks(dir=diffDir,rebuild=T,verbose = T)
cuff<-cuff.hg19
cuff

###########################################################################
### Chunk_2: Get Gene Expression tables with counts and FPKM
############################################################################

gene_exp.diff<-diffData(genes(cuff))
head(gene_exp.diff); dim(gene_exp.diff)
gene.exp.List<-as.data.frame(cbind(logFC=gene_exp.diff$log2_fold_change,
													 p_value=gene_exp.diff$p_value,
													 q_value=gene_exp.diff$q_value),
											 row.names=gene_exp.diff$gene_id)
head(gene.exp.List); dim(gene.exp.List)

g.rep.matrix<-repFpkmMatrix(genes(cuff))
head(g.rep.matrix)
over.group<-grep(pattern = over, x = colnames(g.rep.matrix),ignore.case = T)
under.group<-grep(pattern = under, x = colnames(g.rep.matrix),ignore.case = T)

rep.fpkm.file = file.path("repFPKMmatrix")
write.table(g.rep.matrix,file = rep.fpkm.file, sep = "  ", row.names = T, col.names = T,quote = F)

o.g.rep.matrix<-g.rep.matrix[over.group]
u.g.rep.matrix<-g.rep.matrix[under.group]

g.count.matrix<-repCountMatrix(genes(cuff))
head(g.count.matrix);dim(g.count.matrix)
rep.count.file = file.path("repCOUNTmatrix")
write.table(g.count.matrix,file = rep.count.file, sep = "  ", row.names = T, col.names = T,quote = F)

o.g.count.matrix<-g.count.matrix[over.group]
u.g.count.matrix<-g.count.matrix[under.group]

###########################################################################
### Chunk_2: Get Significant Gene Expression tables with counts and FPKM
############################################################################

mySigGenes<-getSig(cuff,x=over,y=under,alpha=.05,level='genes')
length(mySigGenes); head(mySigGenes)

sigGenes<-getGenes(cuff, mySigGenes)
sigGenes; length(sigGenes)

sig_genes_exp.diff<-diffData(sigGenes)
head(sig_genes_exp.diff); dim(sig_genes_exp.diff)

sig.gene.exp.List<-as.data.frame(cbind(logFC=sig_genes_exp.diff$log2_fold_change,
													 p_value=sig_genes_exp.diff$p_value,
													 q_value=sig_genes_exp.diff$q_value),
											 row.names=sig_genes_exp.diff$gene_id)
head(sig.gene.exp.List); dim(sig.gene.exp.List)

sig_gene.df  <- subset(gene_exp.diff, (log2_fold_change > 2 & q_value <= 0.01))
head(sig_gene_data);dim(sig_gene_data)
sig_gene_ids <-sig_gene_data$gene_id
length(sig_gene_ids)
sigfold <-sig_gene_data$log2_fold_change
s.g.rep.matrix<-subset(g.rep.matrix, row.names(g.rep.matrix) %in% sig_gene_ids)
head(s.g.rep.matrix);dim(s.g.rep.matrix)
sig.rep.fpkm.file = file.path("sig.repFPKMmatrix")
write.table(s.g.rep.matrix,file = sig.rep.fpkm.file, sep = "  ", row.names = T, col.names = T,quote = F)

o.s.g.rep.matrix<-s.g.rep.matrix[over.group]
u.s.g.rep.matrix<-s.g.rep.matrix[under.group]

s.g.count.matrix<-subset(g.count.matrix, row.names(g.count.matrix) %in% sig_gene_ids)
head(s.g.count.matrix);dim(s.g.count.matrix)
sig.rep.count.file = file.path("sig.repCOUNTmatrix")
write.table(s.g.count.matrix,file = sig.rep.count.file, sep = "  ", row.names = T, col.names = T,quote = F)

o.s.g.count.matrix<-s.g.count.matrix[over.group]
u.s.g.count.matrix<-s.g.count.matrix[under.group]

####################################################################################################
### code chunk number : Group significant genes by the condition with higher expression
####################################################################################################

sig_over_gene_data<-subset(gene_exp.diff, (q_value <= 0.01) & (log2_fold_change > 2))
head(sig_over_gene_data)
nrow(sig_over_gene_data)

sig_under_gene_data<-subset(gene_exp.diff, (q_value <= 0.01) & (log2_fold_change < -2))
head(sig_under_gene_data)
nrow(sig_under_gene_data)

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

####################################################################################################
### code chunk number : Step 6: Clustering
##################################################

genes.rep.fpkm.df<-as.data.frame(s.g.rep.matrix)
gse.cor<-cor(genes.rep.fpkm.df)
#plotcorr(sig.gse.cor, mar = c(0.1, 0.1, 0.1, 0.1))
# Do the same, but with colors corresponding to value
colorfun <- colorRamp(c("#CC0000","white","#3366CC"), space="Lab")
plotcorr(gse.cor, col=rgb(colorfun((sig.gse.cor+1)/2), maxColorValue=255),
         mar = c(0.1, 0.1, 0.1, 0.1))

gene.dist<-as.matrix(dist(t(s.g.rep.matrix)))
loc <- cmdscale(gene.dist) # Performs MDS analysis on the geographic distances between European cities.
plot(loc[,1], -loc[,2], type="n", xlab="", ylab="", main="Euclidian MDS analysis")
text(loc[,1], -loc[,2], rownames(loc), cex=1,col=c("blue","darkgreen"))
# Plots the MDS results in 2D plot. The minus is required in this example to flip the plotting orientation.
mydatascale <- t(scale(t(s.g.rep.matrix))) # Centers and scales data.
hr <- hclust(as.dist(1-cor(t(mydatascale), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
hc <- hclust(as.dist(1-cor(mydatascale, method="spearman")), method="complete")
# Clusters columns by Spearman correlation.
heatmap(gene.dist, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), scale="row") # col=my.colorFct(),
# Plot the data table as heatmap and the cluster results as dendrograms.
mycl <- cutree(hr, h=max(hr$height)/1.5)
mycolhc <- sample(rainbow(256))
mycolhc <- mycolhc[as.vector(mycl)]
heatmap(gene.dist, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), scale="row", RowSideColors=mycolhc) #col=my.colorFct(),
# Cut the tree at specific height and color the corresponding clusters in the heatmap color bar.

sig.dist.t<-as.matrix(cor(t(s.g.rep.matrix)))
fit <- cmdscale(sig.dist.t,eig=TRUE, k=2) # k is the number of dim
fit # view results
# plot solution
x <- fit$points[,1];y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS", type="n")
text(x, y, labels = row.names(sig.dist.t),cex=.7,col=c("blue","darkgreen"))

gse.scaled<-as.matrix(scale(t(g.rep.matrix)))
heatmap(gse.scaled, Colv=F, scale='none')

gsecor<-cor(g.rep.matrix)
round(gsecor,2)
# Do the same, but with colors corresponding to value
plot(hclust(dist(g.rep.matrix,method="euclidian"),method="single"))
cl <- kmeans(g.rep.matrix, 2)
plot(g.rep.matrix, col = cl$cluster)
points(cl$centers, col = 1:2, pch = 8, cex=2)

gsecor<-cor(t(s.g.rep.matrix))
round(gsecor,2)
# Do the same, but with colors corresponding to value
plot(hclust(dist(s.g.rep.matrix,method="euclidian"),method="single"))
cl <- kmeans(s.g.rep.matrix, 2)
plot(s.g.rep.matrix, col = cl$cluster)
points(cl$centers, col = 1:2, pch = 8, cex=2)


c <- cor(t(s.g.rep.matrix), method="spearman")
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

###################################################################################################
### code chunk number : Step 6: Statistical Tests on count data
###################################################################################################
u.g.count.matrix<-as.matrix(u.g.count.matrix)
o.g.count.matrix<-as.matrix(o.g.count.matrix)

t.test(x = o.g.count.matrix,
       y = u.g.count.matrix,
       alternative = "two.sided")


cor.test(x = o.g.count.matrix,
       y = u.g.count.matrix,
       alternative = "two.sided")

wilcox.test(o.g.count.matrix,u.g.count.matrix, alternative = "g")  # g for greater

var.test(o.g.count.matrix,u.g.count.matrix)
chisq.test(o.g.count.matrix,u.g.count.matrix)

###################################################################################################
### code chunk number : Step 6: Statistical Tests on FPKM data
###################################################################################################
u.g.rep.ma<-as.matrix(u.g.rep.matrix)
o.g.rep.ma<-as.matrix(o.g.rep.matrix)

var.test(o.g.rep.ma,u.g.count.matrix)
chisq.test(o.g.rep.ma,u.g.count.matrix)

t.test(x = o.g.rep.ma,
       y = u.g.rep.ma,
       alternative = "two.sided")

cor.test(x = o.g.rep.ma,
       y = u.g.rep.ma,
       alternative = "two.sided")

wilcox.test(o.g.rep.ma,u.g.rep.ma, alternative = "g")  # g for greater

####################################################################################################
### code chunk number : Step 3: #Hierarchical clustering
####################################################################################################
library(pvclust); library(gplots) # Loads the required packages.
#Hierarchical clustering
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R")
# Import an alternative color scheme for the heatmap function.
mydatascale <- t(scale(t(u.g.rep.ma))) # Centers and scales data.
hr <- hclust(as.dist(1-cor(t(u.g.rep.ma), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
hc <- hclust(as.dist(1-cor(u.g.rep.ma, method="spearman")), method="complete")
# Clusters columns by Spearman correlation.
heatmap(mydatascale, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=my.colorFct(), scale="row")
# Plot the data table as heatmap and the cluster results as dendrograms.
mycl <- cutree(hr, h=max(hr$height)/1.5)
mycolhc <- sample(rainbow(256))
mycolhc <- mycolhc[as.vector(mycl)]
heatmap(sig.genes.gse116969, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),
        col=my.colorFct(), scale="row", RowSideColors=mycolhc)
# Cut the tree at specific height and color the corresponding clusters in the heatmap color bar.

pv <- pvclust(scale(t(o.g.rep.ma)), method.dist="correlation", method.hclust="complete", nboot=10)
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


geneList<-mapIds(x = org.Hs.eg.db,
                 keys =  gene_exp.diff$gene_id,
                 column = "ENTREZID",
                 keytype = "SYMBOL",
                 multiVals="first")
geneList

s.g.rep.matrix<-as.matrix(s.g.rep.matrix)
####################################################################################################
### code chunk number : Step 3: #Hierarchical clustering
####################################################################################################
library(pvclust); library(gplots) # Loads the required packages.

#Hierarchical clustering

source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R")
# Import an alternative color scheme for the heatmap function.
mydatascale <- t(scale(t(s.g.rep.matrix))) # Centers and scales data.
hr <- hclust(as.dist(1-cor(t(mydatascale), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
hc <- hclust(as.dist(1-cor(mydatascale, method="spearman")), method="complete")
# Clusters columns by Spearman correlation.
heatmap(s.g.rep.matrix, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=my.colorFct(), scale="row")
# Plot the data table as heatmap and the cluster results as dendrograms.
mycl <- cutree(hr, h=max(hr$height)/1.5)
mycolhc <- sample(rainbow(256))
mycolhc <- mycolhc[as.vector(mycl)]
heatmap(s.g.rep.matrix, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),
        col=my.colorFct(), scale="row", RowSideColors=mycolhc)
# Cut the tree at specific height and color the corresponding clusters in the heatmap color bar.

pv <- pvclust(scale(t(s.g.rep.matrix)), method.dist="correlation", method.hclust="complete", nboot=10)
# Perform the hierarchical cluster analysis.
#Due to time resrictions, we are using here only 10 bootstrap repetitions.
# Usually, one should use at least 1000 repetitions.
plot(pv, hang=-1); pvrect(pv, alpha=0.95)
# Plots result as a dendrogram where the significant clusters are highlighted with red rectangles.
clsig <- unlist(pvpick(pv, alpha=0.95, pv="au", type="geq", max.only=TRUE)$clusters) # Retrieve members of significant clusters.

source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/dendroCol.R") # Import tree coloring function.
dend_colored <- dendrapply(as.dendrogram(pv$hclust), dendroCol, keys=clsig, xPar="edgePar", bgr="black", fgr="red", pch=20)
# Create dendrogram object where the significant clusters are labeled in red.
heatmap(s.g.rep.matrix, Rowv=dend_colored, Colv=as.dendrogram(hc), col=my.colorFct(), scale="row", RowSideColors=mycolhc)
# Plot the heatmap from above, but with the significant clusters in red and the cluster bins from the tree cutting step in
# the color bar.


########################################################################
### code chunk number 4: Write adj mtrix and plot graph for over genes
########################################################################

#Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
over.graph <- graph.adjacency(as.matrix(as.dist(cor(t(sig.g.u.rep.matrix), method="pearson"))),
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
vSizes <- (scale01(apply(o.s.g.rep.matrix, 1, mean)) + 1.0) * 10

#Amplify or decrease the width of the edges
edgeweights <- E(over.graph)$weight * 2.0
#Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst <- mst(over.graph, algorithm="prim")

plot(mst,vertex.size=vSizes, vertex.label.dist=-0.5, vertex.label.color="black",
     asp=FALSE, vertex.label.cex=0.6, edge.width=edgeweights, edge.arrow.mode=0,
     main="Network of significant genes")

#Plot the tree object
plot(mst, layout=layout.fruchterman.reingold, edge.curved=TRUE,
     vertex.size=vSizes, vertex.label.dist=-0.5, vertex.label.color="black",
     asp=FALSE, vertex.label.cex=0.6, edge.width=edgeweights, edge.arrow.mode=0,
     main="Network of significant genes")

########################################################################
### code chunk number 5: Write adj mtrix and plot graph for under genes
########################################################################

#Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
under.graph <- graph.adjacency(as.matrix(as.dist(cor(t(sig.g.u.rep.matrix), method="pearson"))),
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
     main="Under group's significant gene network")

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


####################################################################################################
### code chunk number : Step 2, Read in, normalise, and identify genes with significant effects
####################################################################################################



gene_diff_data<-diffData(genes(cuff))
head(gene_diff_data)
dim(gene_diff_data)

sig_gene_data<-subset(gene_diff_data, (significant =='yes' & q_value < 0.005 & abs(log2_fold_change) > 2))
head(sig_gene_data)
dim(sig_gene_data)
genes<-sig_gene_data[,"gene_id"]
fold.change<-sig_gene_data[,"log2_fold_change"]
sig.gene.List<-as.data.frame(fold.change, row.names = genes)
s.gene.List<-as.matrix(sig.gene.List)
head(s.gene.List)
s.gene.order<-s.gene.List[order(s.gene.List[,"fold.change"], decreasing = T),]
head(s.gene.order)

de <-names(s.gene.order)[abs(s.gene.order) > 1.5]
head(de)
length(de)

eg = bitr(de, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
sym.eg<-as.factor(eg$SYMBOL)
head(eg)
head(sym.eg)

symbols<-names(s.gene.order)
head(symbols)
length(symbols)

entrez.syms<-which(names(s.gene.order) %in% eg$SYMBOL)
length(entrez.syms)
head(entrez.syms)
head(sig.gene.fold.list)
s.gene.order.eg<-s.gene.order[entrez.syms]
s.gene.order.syeg<-cbind(eg,s.gene.order.eg)
head(s.gene.order.syeg)
log2FC<-s.gene.order.syeg[,"s.gene.order.eg"]
SigGeneList<-as.data.frame(log2FC, row.names=s.gene.order.syeg$ENTREZID)
head(SigGeneList)

de <- row.names(SigGeneList)
x <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
head(as.data.frame(x))

## ----fig.height=6, fig.width=12------------------------------------------
barplot(x, showCategory=8)
dotplot(x, showCategory=15)
cnetplot(x, categorySize="pvalue", foldChange=fold.change)

SigGene.ma<-as.matrix(SigGeneList)
row.names(SigGene.ma)
SigGene.ma
## ------------------------------------------------------------------------
y <- gsePathway(row.names(SigGene.ma), nPerm=10000,
                pvalueCutoff=0.2,
                pAdjustMethod="BH", verbose=FALSE)
res <- as.data.frame(y); head(res)

## ----fig.height=8, fig.width=8-------------------------------------------
#emapplot(y, color="pvalue")
#gseaplot(y, geneSetID = "R-HSA-69242")
#viewPathway("E2F mediated regulation of DNA replication", readable=TRUE, foldChange=geneList)
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
head(ggoCC)
## ----fig.height=5, fig.width=9-------------------------------------------
barplot(ggoCC, drop=TRUE, showCategory=12)
barplot(ggoMF, drop=TRUE, showCategory=12)
barplot(ggoBP, drop=TRUE, showCategory=12)

barplot(ego, showCategory=8)
dotplot(ggoCC, colour=qvalue)
head(ggoCC)
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
cnetplot(ego, categorySize="pvalue", foldChange=SigGeneList)

## ----fig.height=12, fig.width=8, eval=FALSE------------------------------
goplot(ego)

## ----eval=FALSE----------------------------------------------------------
  ego2 <- enrichGO(gene         = SigGeneList,
                  OrgDb         = org.Hs.eg.db,
#                  keyType       = 'ENTREZID',
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)

## ----eval=FALSE----------------------------------------------------------
  ego2 <- setReadable(ego2, OrgDb = org.Hs.eg.db)

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

