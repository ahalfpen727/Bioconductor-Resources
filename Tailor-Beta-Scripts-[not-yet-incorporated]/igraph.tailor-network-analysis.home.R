library(igraph); library(Hmisc); library(cummeRbund)
library(ellipse); library(limma); library(edgeR)
library(pvclust); library(gplots); library(Hmisc)
library(readxl); library(knitr);library(DOSE);
library(GO.db);library(org.Hs.eg.db);library(GSEABase)
library(clusterProfiler);library(ReactomePA)
data(geneList)
#######################################################################
######Tutorial: Network plot from expression data in R using igraph
###### igraph netowrk analysis for cuffdiff data
############################################################################
#https://www.biostars.org/p/285296/
   # These arguments are passed in on the command line via launchR.eachDiffDir.sh
   # args:
   # diffDir=\"${diffDir}\" inDir=\"${INPUT}/${diffDir}\" outDir=\"${OUTPUT}/${diffDir}\"
   # under=\"${under}\" over=\"${over}\" Rplots=\"$\{OUTPUT}/${diffDir}/${R_PLOTS}\"
   # FPKMmatrix=\"${OUTPUT}/${diffDir}/${FPKM_MATRIX}\" DiffTable=\"${OUTPUT}/${diffDir}/${DIFF_TABLE}\"
#library(devtools)#install_github("genomicsclass/maPooling")
############################################################################
### code chunk number 1: Load Libraries and initialize cuff object
############################################################################
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
diffDir<-file.path("../cuffdiff_results_hg38_default/LUTS-over-CTRL")
setwd(diffDir);dir()
over="LUTS";under="CTRL"
cuff=readCufflinks(dir=".",rebuild=F,verbose = T)
cuff

runInfo(cuff)
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
####################################################################################################
### code chunk number : Step 3: #Hierarchical clustering
####################################################################################################
pdf("Significant-Gene-Subnetworks-by-Condition.pdf")
dev.off()
library(pvclust); library(gplots) # Loads the required packages.
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/dendroCol.R") # Import tree coloring function.
# Import an alternative color scheme for the heatmap function.
sig.g.rep.ma<-as.matrix(sig.g.rep.matrix)
mydatascale <- t(scale(t(sig.g.rep.ma))) # Centers and scales data.
hr <- hclust(as.dist(1-cor(t(mydatascale), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
hc <- hclust(as.dist(1-cor(mydatascale, method="spearman")), method="complete")
# Clusters columns by Spearman correlation.
heatmap(sig.g.rep.ma, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),
        scale="row", main="Clustering of all Sig Diff Expr Genes")
# Plot the data table as heatmap and the cluster results as dendrograms.
mycl <- cutree(hr, h=max(hr$height)/1.5)
mycolhc <- sample(rainbow(256))
mycolhc <- mycolhc[as.vector(mycl)]

heatmap(sig.g.rep.ma, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),
        scale="row", RowSideColors=mycolhc, main="Clustering of all Sig Diff Expr Genes")

pv <- pvclust(scale(sig.g.rep.matrix), method.dist="correlation", method.hclust="complete", nboot=10)
plot(pv, hang=-1, main="Clustering by Sample with all Sig Diff Expr Genes"); pvrect(pv, alpha=0.95)
clsig <- unlist(pvpick(pv, alpha=0.95, pv="au", type="geq", max.only=TRUE)$clusters) # Retrieve members of significant clusters.
dend_colored <- dendrapply(as.dendrogram(pv$hclust), dendroCol, keys=clsig, xPar="edgePar", bgr="black", fgr="red", pch=20)
# Create dendrogram object where the significant clusters are labeled in red.
sig.g.rep.ma<-as.matrix(sig.g.rep.matrix)
# Plot the heatmap from above, but with the significant clusters in red and the cluster bins from the tree cutting step in
# the color bar.
heatmap.2(sig.g.rep.ma, Rowv=dend_colored, Colv=as.dendrogram(hc), col=c("green", "red"), # my.colorFct(),
          scale="row", trace="none") # RowSideColors=mycolhc) # Plot heatmap with heatmap.2() function which scales better for


pv <- pvclust(scale(dist(sig.g.rep.matrix)), use.cor="all.obs",method.dist="correlation", method.hclust="complete", nboot=10)
plot(pv, hang=-1,  main="Bootstrapping the Scaled Distances"); pvrect(pv, alpha=0.95)
clsig <- unlist(pvpick(pv, alpha=0.95, pv="au", type="geq", max.only=TRUE)$clusters) # Retrieve members of significant clusters.
dend_colored <- dendrapply(as.dendrogram(pv$hclust), dendroCol, keys=clsig, xPar="edgePar", bgr="black", fgr="red", pch=20)
heatmap(sig.g.rep.ma, Rowv=dend_colored, Colv=as.dendrogram(hc),scale="row",
        main="Clustering the FPKM values all Sig Diff Expr Genes")#, RowSideColors=mycolhc)
heatmap(sig.g.rep.ma, Rowv=dend_colored, Colv=as.dendrogram(hc),  scale="row",
        main="Clustering the FPKM values all Sig Diff Expr Genes")#, RowSideColors=mycolhc)) #, RowSideColors=mycolhc)
heatmap.2(sig.g.rep.ma, Rowv=dend_colored, Colv=as.dendrogram(hc), col=c("green", "red"), # my.colorFct(),
          scale="row", trace="none",  main="Clustering the FPKM values all the Sig Diff Expr Genes")#, RowSideColors=mycolhc))

pv <- pvclust(scale(t(sig.g.rep.matrix)), method.dist="correlation", method.hclust="complete", nboot=10)
plot(pv, hang=-1, main="Bootstrapping of all Sig Diff Expr Genes"); pvrect(pv, alpha=0.95)
clsig <- unlist(pvpick(pv, alpha=0.95, pv="au", type="geq", max.only=TRUE)$clusters) # Retrieve members of significant clusters.
dend_colored <- dendrapply(as.dendrogram(pv$hclust), dendroCol, keys=clsig, xPar="edgePar", bgr="black", fgr="red", pch=20)
# Create dendrogram object where the significant clusters are labeled in red.
sig.g.rep.ma<-as.matrix(sig.g.rep.matrix)
heatmap(sig.g.rep.ma, Rowv=dend_colored, Colv=as.dendrogram(hc), scale="row",
        main="Clustering of all Sig Diff Expr Genes")#, RowSideColors=mycolhc)
heatmap(sig.g.rep.ma, Rowv=dend_colored, Colv=as.dendrogram(hc), scale="row",
        main="Clustering of all Sig Diff Expr Genes") #, RowSideColors=mycolhc)
heatmap.2(sig.g.rep.ma, Rowv=dend_colored, Colv=as.dendrogram(hc), col=c("green", "red"), # my.colorFct(),
          scale="row", trace="none", main="Clustering the FPKM values of all the Sig Diff Expr Genes") # RowSideColors=mycolhc) #
########################################################################
### code chunk number 4: Write adj mtrix and plot graph for over genes
########################################################################
#Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
all.sig.graph <- graph.adjacency(as.matrix(as.dist(cor(t(sig.g.rep.matrix), method="pearson"))),
                               mode="undirected", weighted=TRUE, diag=FALSE)
#NB - Euclidean distances can also be used instead of correlation, but this changes the tutorial slightly:
#all.sig.graph <- graph.adjacency(as.matrix(dist(sig.g.rep.matrix, method="euclidean")),
#                                 mode="undirected", weighted=TRUE, diag=FALSE)
#g <- graph.adjacency(as.matrix(dist(estrogenMainEffects, method="euclidean")),
#Simplfy the adjacency object
#all.sig.graph <- simplify(sig.g.rep.matrix, remove.multiple=TRUE, remove.loops=TRUE)
#Colour negative correlation edges as blue
E(all.sig.graph)[which(E(all.sig.graph)$weight<0)]$color <- "darkblue"
#Colour positive correlation edges as red
E(all.sig.graph)[which(E(all.sig.graph)$weight>0)]$color <- "darkred"
#Convert edge weights to absolute values
E(all.sig.graph)$weight <- abs(E(all.sig.graph)$weight)
#Change arrow size
#For directed graphs only
#E(all.sig.graph)$arrow.size <- 1.0
#Remove edges below absolute Pearson correlation 0.8
all.sig.graph <- delete_edges(all.sig.graph, E(all.sig.graph)[which(E(all.sig.graph)$weight<0.8)])
#Assign names to the graph vertices (optional)
V(all.sig.graph)$name <- V(all.sig.graph)$name
#Change shape of graph vertices
V(all.sig.graph)$shape <- "sphere"
#Change colour of graph vertices
V(all.sig.graph)$color <- "skyblue"
#Change colour of vertex frames
V(all.sig.graph)$vertex.frame.color <- "white"
#Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
#Multiple scaled vales by a factor of 10
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
vSizes <- (scale01(apply(sig.g.rep.matrix, 1, mean)) + 1.0) * 10
#Amplify or decrease the width of the edges
edgeweights <- E(all.sig.graph)$weight * 2.0
#Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst <- mst(all.sig.graph, algorithm="prim")
#Plot the tree object
plot(mst, layout=layout.fruchterman.reingold, edge.curved=TRUE,
     vertex.size=vSizes, vertex.label.dist=-0.5, vertex.label.color="black",
     asp=FALSE, vertex.label.cex=0.6, edge.width=edgeweights, edge.arrow.mode=0,
     main="Gene Network of all the Significantly Differentially Expressed Genes")

#################################################################################################
### code chunk number 4: Identify communities in the tree object based on 'edge betweenness'
###############################################################################################
mst.communities <- edge.betweenness.community(mst, directed=T)
#mst.communities <- edge.betweenness.community(all.sig.graph, directed=FALSE)
mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
V(mst)$color <- mst.communities$membership + 1
#par(mfrow=c(1,2))
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
   main="Gene network with all Significant Differentially Expressed Genes"
)

plot(
   mst.clustering,mst,
   layout=layout.fruchterman.reingold,
   edge.curved=TRUE,
   vertex.size=vSizes,
   vertex.label.dist=-0.5,
   vertex.label.color="black",
   asp=FALSE,
   vertex.label.cex=0.6,
   edge.width=edgeweights,
   edge.arrow.mode=0,
   main="Gene network with all Significant Differentially Expressed Genes"
)

########################################################################
### code chunk number 4: Write adj mtrix and plot graph for over genes
########################################################################

#Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
over.graph <- graph.adjacency(as.matrix(as.dist(cor(t(over.grp.s.o.g.rep.ma), method="pearson"))),
                              mode="undirected", weighted=TRUE, diag=FALSE)
#NB - Euclidean distances can also be used instead of correlation, but this changes the tutorial slightly:
#g <- graph.adjacency(as.matrix(dist(estrogenMainEffects, method="euclidean")),
#g <- graph.adjacency(as.matrix(dist(s.u.g.rep.matrix, method="euclidean")),mode="undirected", weighted=TRUE, diag=FALSE)
#Simplfy the adjacency object
#over.graph <- simplify(over.graph, remove.multiple=TRUE, remove.loops=TRUE)
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
vSizes <- (scale01(apply(over.grp.s.o.g.rep.ma, 1, mean)) + 1.0) * 10
#Amplify or decrease the width of the edges
edgeweights <- E(over.graph)$weight * 2.0
#Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst <- mst(over.graph, algorithm="prim")
#Plot the tree object
plot(mst, layout=layout.fruchterman.reingold, edge.curved=TRUE,
     vertex.size=vSizes, vertex.label.dist=-0.5, vertex.label.color="black",
     asp=FALSE, vertex.label.cex=0.6, edge.width=edgeweights, edge.arrow.mode=0,
     main="Network of significant highly expressed genes in over group")

#################################################################################################
### code chunk number 4: Identify communities in the tree object based on 'edge betweenness'
###############################################################################################

mst.communities <- edge.betweenness.community(mst)
mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
V(mst)$color <- mst.communities$membership + 1

#par(mfrow=c(1,2))
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
   main="Over Group's Significant Highly Expressed Gene Network"
)

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
   main="Over Group's Significant Highly Expressed Gene Network"
)

########################################################################
### code chunk number 4: Write adj mtrix and plot graph for over genes
########################################################################

#Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
under.graph <- graph.adjacency(as.matrix(as.dist(cor(t(under.grp.s.u.g.rep.ma), method="pearson"))),
                              mode="undirected", weighted=TRUE, diag=FALSE)
#NB - Euclidean distances can also be used instead of correlation, but this changes the tutorial slightly:
#g <- graph.adjacency(as.matrix(dist(estrogenMainEffects, method="euclidean")),
#g <- graph.adjacency(as.matrix(dist(s.u.g.rep.matrix, method="euclidean")),mode="undirected", weighted=TRUE, diag=FALSE)
#Simplfy the adjacency object
#over.graph <- simplify(over.graph, remove.multiple=TRUE, remove.loops=TRUE)
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
under.graph <- delete_edges(under.graph, E(under.graph)[which(E(under.graph)$weight<0.9)])
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
vSizes <- (scale01(apply(under.grp.s.u.g.rep.ma, 1, mean)) + 1.0) * 10
#Amplify or decrease the width of the edges
edgeweights <- E(under.graph)$weight * 2.0
#Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst <- mst(under.graph, algorithm="prim")
#Plot the tree object
plot(mst, layout=layout.fruchterman.reingold, edge.curved=TRUE,
     vertex.size=vSizes, vertex.label.dist=-0.5, vertex.label.color="black",
     asp=FALSE, vertex.label.cex=0.6, edge.width=edgeweights, edge.arrow.mode=0,
     main="Network of significant highly expressed genes in under group")

#################################################################################################
### code chunk number 4: Identify communities in the tree object based on 'edge betweenness'
###############################################################################################

mst.communities <- edge.betweenness.community(mst, directed=FALSE)
mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
V(mst)$color <- mst.communities$membership + 1
#par(mfrow=c(1,2))
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
   main="Under Group's Significant Highly Expressed Gene Network"
)
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
   main="Under Group's Significant Highly Expressed Gene Network"
)
dev.off()
#################################################################################################
### code chunk number 5: Step 4, Further analyses
###############################################################################################
head(sig_genes_exp.diff);dim(sig_genes_exp.diff)

#genes_exp.diff<-as.data.frame(cbind(genes=sig_genes_exp.diff$gene_id,logFC=sig_genes_exp.diff$log2_fold_change,p_value=sig_genes_exp.diff$p_value,
#                                    q_value=sig_genes_exp.diff$q_value))
genes_exp.diff<-as.data.frame(cbind(logFC=sig_genes_exp.diff$log2_fold_change,p_value=sig_genes_exp.diff$p_value,
                                    q_value=sig_genes_exp.diff$q_value),row.names=sig_genes_exp.diff$gene_id)
head(genes_exp.diff)
dim(genes_exp.diff)
#gene_exp.diff<-unique(rownames(genes_exp.diff))
#gene.exp.diff<-genes_exp.diff[gene_exp.diff,]

over.group<-grep(pattern = over, x = colnames(g.rep.matrix),ignore.case = T)
under.group<-grep(pattern = under, x = colnames(g.rep.matrix),ignore.case = T)
under.grp.s.u.g.rep.ma

head(under.grp.s.u.g.rep.ma)

under.geneList<-mapIds(x = org.Hs.eg.db,
                 keys =  rownames(under.grp.s.u.g.rep.ma),
                 column = "ENTREZID",
                 keytype = "SYMBOL",
                 multiVals="first")
under.geneList<-na.exclude(under.geneList)
under.geneList<-as.data.frame(under.geneList, row.names=c(rownames(under.geneList)))
colnames(under.geneList)="ENTREZID"
under.geneList

geneList<-mapIds(x = org.Hs.eg.db,
                 keys =  rownames(genes_exp.diff),
                 column = "ENTREZID",
                 keytype = "SYMBOL",
                 multiVals="first")

sig.g.reps<-which(row.names(g.rep.matrix) %in% sigGenes@ids)
sig.g.rep.matrix<-as.matrix(g.rep.matrix[sig.g.reps,])
head(sig.g.rep.matrix)

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
head(sig.gene.List)

sig.gene.order<-order(fold.change, decreasing = T)
fold<-sig.gene.List[sig.gene.order,]
genes<-rownames(sig.gene.List)

sig_genes<-which(gene_diff_data[,"significant"] == 'yes' & gene_diff_data[,"q_value"] < 0.005 & abs(gene_diff_data[,"log2_fold_change"]) > 2)
length(sig_genes)
head(gene_diff_data[sig_genes,])
logfold=gene_diff_data[sig_genes,c("log2_fold_change")]
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

de <- row.names(SigGeneList)
x <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
head(as.data.frame(x))

## ----fig.height=6, fig.width=12------------------------------------------
barplot(x, showCategory=8)
dotplot(x, showCategory=15)
cnetplot(x, categorySize="pvalue", foldChange=fold.change)

## ------------------------------------------------------------------------
#y <- gsePathway(cuff.fc, nPerm=10000,
#                pvalueCutoff=0.2,
#                pAdjustMethod="BH", verbose=FALSE)
#res <- as.data.frame(y); head(res)

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

####################################################################################################
### code chunk number : Step 3, Statistical tests
####################################################################################################

gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
sh <- apply(golub[,gol.fac=="ALL"], 1, function(x) shapiro.test(x)$p.value)
> sum(sh > 0.05)/nrow(golub) * 100
grch38anov<-anova(lm(rowMeans(grch38.lutsgrp) ~ rowMeans(grch38.ctrlgrp)))
grubbs.test(golub[1042, gol.fac=="ALL"])

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

rep.gene.counts<-repCountMatrix(genes(cuff))
head(rep.gene.counts)
dim(rep.gene.counts)

sig.rep.gene.counts<-rep.gene.counts[which(row.names(rep.gene.counts) %in% sig_gene_data$gene_id),]
head(sig.rep.gene.counts)
dim(sig.rep.gene.counts)

rep.count.file = file.path("repCOUNTmatrix")
write.table(rep.gene.counts,file = rep.count.file, sep = "  ", row.names = T, col.names = T,quote = F)

sig.rep.count.file = file.path("sig.repCOUNTmatrix")
write.table(sig.rep.gene.counts,file = sig.rep.count.file, sep = "  ", row.names = T, col.names = T,quote = F)


