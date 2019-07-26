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

############################################################################
### code chunk number 1: Load Libraries and initialize cuff object
############################################################################
library(igraph); library(Hmisc); library(cummeRbund)
library(ellipse); library(limma); library(edgeR)
library(pvclust); library(gplots); library(Hmisc)
library(readxl); library(knitr);library(DOSE);
library(GO.db);library(org.Hs.eg.db);library(GSEABase)
library(clusterProfiler);library(ReactomePA)
data(geneList)

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
genes_exp.diff<-diffTable(genes(cuff))
genes_exp.diff<-diffData(isoforms(cuff))
head(genes_exp.diff);dim(genes_exp.diff)
g.rep.matrix<-repFpkmMatrix(isoforms(cuff))
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

library(pvclust); library(gplots) # Loads the required packages.
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/dendroCol.R") # Import tree coloring function.
# Import an alternative color scheme for the heatmap function.
sig.g.rep.ma<-as.matrix(sig.g.rep.matrix)
mydatascale <- t(scale(t(sig.g.rep.ma))) # Centers and scales data.
hr <- hclust(as.dist(1-cor(t(mydatascale), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
hc <- hclust(as.dist(1-cor(mydatascale, method="spearman")), method="complete")
# Clusters columns by Spearman correlation.
heatmap(sig.g.rep.ma, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=my.colorFct(),
        scale="row", main="Clustering of all Sig Diff Expr Genes")
# Plot the data table as heatmap and the cluster results as dendrograms.
mycl <- cutree(hr, h=max(hr$height)/1.5)
mycolhc <- sample(rainbow(256))
mycolhc <- mycolhc[as.vector(mycl)]
heatmap(sig.g.rep.ma, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),
        col=my.colorFct(), scale="row", RowSideColors=mycolhc, main="Clustering of all Sig Diff Expr Genes")

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
heatmap(sig.g.rep.ma, Rowv=dend_colored, Colv=as.dendrogram(hc), col=my.colorFct(), scale="row",
        main="Clustering the FPKM values all Sig Diff Expr Genes")#, RowSideColors=mycolhc)
heatmap(sig.g.rep.ma, Rowv=dend_colored, Colv=as.dendrogram(hc), col=my.colorFct(), scale="row",
        main="Clustering the FPKM values all Sig Diff Expr Genes")#, RowSideColors=mycolhc)) #, RowSideColors=mycolhc)
heatmap.2(sig.g.rep.ma, Rowv=dend_colored, Colv=as.dendrogram(hc), col=c("green", "red"), # my.colorFct(),
          scale="row", trace="none",  main="Clustering the FPKM values all the Sig Diff Expr Genes")#, RowSideColors=mycolhc))

pv <- pvclust(scale(t(sig.g.rep.matrix)), method.dist="correlation", method.hclust="complete", nboot=10)
plot(pv, hang=-1, main="Bootstrapping of all Sig Diff Expr Genes"); pvrect(pv, alpha=0.95)
clsig <- unlist(pvpick(pv, alpha=0.95, pv="au", type="geq", max.only=TRUE)$clusters) # Retrieve members of significant clusters.
dend_colored <- dendrapply(as.dendrogram(pv$hclust), dendroCol, keys=clsig, xPar="edgePar", bgr="black", fgr="red", pch=20)
# Create dendrogram object where the significant clusters are labeled in red.
sig.g.rep.ma<-as.matrix(sig.g.rep.matrix)
heatmap(sig.g.rep.ma, Rowv=dend_colored, Colv=as.dendrogram(hc), col=my.colorFct(), scale="row",
        main="Clustering of all Sig Diff Expr Genes")#, RowSideColors=mycolhc)
heatmap(sig.g.rep.ma, Rowv=dend_colored, Colv=as.dendrogram(hc), col=my.colorFct(), scale="row",
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
                                 mode="undirected", weighted=TRUE, diag=FALSE)
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
head(genes_exp.diff);dim(genes_exp.diff)

#genes_exp.diff<-as.data.frame(cbind(genes=sig_genes_exp.diff$gene_id,logFC=sig_genes_exp.diff$log2_fold_change,p_value=sig_genes_exp.diff$p_value,
#                                    q_value=sig_genes_exp.diff$q_value))
sig.genes_exp.diff<-as.data.frame(cbind(logFC=sig_genes_exp.diff$log2_fold_change,p_value=sig_genes_exp.diff$p_value,
                                    q_value=sig_genes_exp.diff$q_value),row.names=sig_genes_exp.diff$gene_id)
head(sig.genes_exp.diff)
dim(sig.genes_exp.diff)
#gene_exp.diff<-unique(rownames(genes_exp.diff))
#gene.exp.diff<-genes_exp.diff[gene_exp.diff,]

over.group<-grep(pattern = over, x = colnames(g.rep.matrix),ignore.case = T)
under.group<-grep(pattern = under, x = colnames(g.rep.matrix),ignore.case = T)

head(sig.g.rep.matrix)
under.grp.s.u.g.rep.ma<-s.u.g.cnt.matrix[,under.group]
over.grp.s.o.g.rep.ma<-s.o.g.cnt.matrix[,over.group]

gene.List<-mapIds(x = org.Hs.eg.db,
                 keys =  rownames(sig.g.rep.matrix),
                 column = "ENTREZID",
                 keytype = "SYMBOL",
                 multiVals="first")
#gene.List<-na.exclude(geneList)
genelist<-as.data.frame(gene.List, row.names=c(rownames(gene.List)))
colnames(genelist)="ENTREZID"
as.data.frame(genelist)
sig.g.rep.eg.ma<-as.data.frame(cbind(genelist, sig.g.rep.matrix))
sig.g.rep.eg.ma<-na.delete(sig.g.rep.eg.ma)
row.names(sig.g.rep.eg.ma)=sig.g.rep.eg.ma$ENTREZID
head(sig.g.rep.eg.ma); dim(sig.g.rep.eg.ma)
sig.g.rep.eg.ma<-sig.g.rep.eg.ma[,-1]
head(sig.g.rep.eg.ma); dim(sig.g.rep.eg.ma)

geneList<-mapIds(x = org.Hs.eg.db,
                 keys =  rownames(sig.genes_exp.diff),
                 column = "ENTREZID",
                 keytype = "SYMBOL",
                 multiVals="first")

genelist<-as.data.frame(geneList, row.names=c(rownames(geneList)))
colnames(genelist)="ENTREZID"
as.data.frame(genelist)
sig.genes_exp.eg.df<-as.data.frame(cbind(genelist, sig.genes_exp.diff))
sig.genes_exp.eg.df<-na.delete(sig.genes_exp.eg.df)
row.names(sig.genes_exp.eg.df)=sig.genes_exp.eg.df$ENTREZID
head(sig.genes_exp.eg.df); dim(sig.genes_exp.eg.df)
sig.genes_exp.eg.df<-sig.genes_exp.eg.df[,-1]
head(sig.genes_exp.eg.df); dim(sig.genes_exp.eg.df)

sig.genes_exp.eg.df<-subset(sig.genes_exp.eg.df, (q_value < 0.05 & abs(logFC) > 2))
head(sig.genes_exp.eg.df)
dim(sig.genes_exp.eg.df)

sig.gene.order<-order(sig.genes_exp.eg.df$logFC, decreasing = T)
fold<-sig.genes_exp.eg.df[sig.gene.order,]
summary(fold$logFC)
head(fold)

maxnum<-max(fold$logFC[fold$logFC!=max(fold$logFC)] )
minnum<-min(fold$logFC[fold$logFC!=min(fold$logFC)] )
max.num<-c(maxnum+2)
min.num<-c(minnum-2)

log.FC<-fold$logFC
head(log.FC)

fold[which(fold == Inf),1] <- max.num
fold[which(fold == -Inf),1] <- min.num
head(fold)

sorted.fold.vals<-order(fold$logFC, decreasing = T)
gene.fold.eg.df<-fold[sorted.fold.vals,]
head(gene.fold.eg.df)
gene.foldchange.eg.df<-as.data.frame(cbind(ENTREZIDS=c(rownames(gene.fold.eg.df)), log2FC=gene.fold.eg.df[,1]))
head(gene.foldchange.eg.df)

x <- enrichPathway(gene=gene.foldchange.eg.df$ENTREZIDS,pvalueCutoff=0.05, readable=T)
head(as.data.frame(x))

## ----fig.height=6, fig.width=12------------------------------------------
barplot(x, showCategory=8)
dotplot(x, showCategory=15)
cnetplot(x, categorySize="pvalue")

## ------------------------------------------------------------------------
gene.foldchange.eg.ma<-gene.foldchange.eg.df$log2FC
names(gene.foldchange.eg.ma)<-gene.foldchange.eg.df$ENTREZIDS
head(gene.foldchange.eg.ma)

y <- gsePathway(gene.foldchange.eg.ma, nPerm=10000,
                pvalueCutoff=0.2, pAdjustMethod="BH", verbose=FALSE)
res <- as.data.frame(y); head(res)

#viewPathway("E2F mediated regulation of DNA replication", readable=TRUE, foldChange=geneList)
ggoCC <- groupGO(gene     = row.names(fold),
                 OrgDb    = org.Hs.eg.db,
                 ont      = "CC",
                 level    = 3,
                 readable = TRUE)

ggoBP <- groupGO(gene     = row.names(fold),
                 OrgDb    = org.Hs.eg.db,
                 ont      = "BP",
                 level    = 3,
                 readable = TRUE)

ggoMF <- groupGO(gene     = row.names(fold),
                 OrgDb    = org.Hs.eg.db,
                 ont      = "MF",
                 level    = 3,
                 readable = TRUE)
head(ggoCC);head(ggoMF);head(ggoBP)

barplot(ggoCC, drop=TRUE, showCategory=10,main="Cellular Component Enrichment")
barplot(ggoMF, drop=TRUE, showCategory=10, main="Molecular Function Enrichment")
barplot(ggoBP, drop=TRUE, showCategory=10, main="Biological Process Enrichment")

barplot(ggoCC, drop=TRUE, showCategory=30,main="Cellular Component Enrichment")
barplot(ggoMF, drop=TRUE, showCategory=30, main="Molecular Function Enrichment")
barplot(ggoBP, drop=TRUE, showCategory=30, main="Biological Process Enrichment")

## ------------------------------------------------------------------------
genes_exp.diff<-diffData(genes(cuff))
head(genes_exp.diff);dim(genes_exp.diff)

over.gene.df <- bitr(row.names(over.grp.s.o.g.rep.ma), fromType = "SYMBOL",
                     toType = c("ENTREZID", "SYMBOL"),OrgDb = org.Hs.eg.db)

under.gene.df <- bitr(row.names(under.grp.s.u.g.rep.ma), fromType = "SYMBOL",
                      toType = c("ENTREZID", "SYMBOL"),OrgDb = org.Hs.eg.db)
head(under.gene.df)
dim(under.grp.s.u.g.rep.ma)
dim(over.grp.s.o.g.rep.ma)

en.go.mf <- enrichGO(gene = row.names(fold),
                universe      = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
goplot(en.go.bp)
barplot(en.go.bp, showCategory=20)
dotplot(en.go.bp, showCategory=20)
dotplot(en.go.bp, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
## remove redundent GO terms
ego2 <- simplify(ego)
cnetplot(en.go.bp, foldChange=fold$logFC)
dev.off()
en.go.cc <- enrichGO(gene     = row.names(fold),
                universe      = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

en.go.bp <- enrichGO(gene         = row.names(fold),
                     universe      = gene.df$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

