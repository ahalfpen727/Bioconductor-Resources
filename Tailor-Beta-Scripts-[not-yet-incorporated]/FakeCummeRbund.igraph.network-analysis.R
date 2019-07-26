#### Import and analyze gene expression data from the Tophat/Cufflinks/Cuffdiff pipeline without cummerbund
##################################################
### LibraryLoading
###################################################
# 000006965* sme  Sinorhizobium meliloti 1021 --> NC_003047 11474104,11481430
### R code from vignette source 'STRINGdb.Rnw'
library(GO.db);library(GOstats);library(topGO)
library(gage);library(pathview);library(org.Hs.eg.db)
library(ReactomePA);library(Rsubread); library(edgeR)
library(Path2PPI); library(STRINGdb)
library(KEGG.db);library(KEGGgraph);library(keggorthology)
library(clusterProfiler);library(DOSE);library(limma)
library(igraph);source("http://igraph.sf.net")
###################################################
################################################################################################################
# Load-libraries ####
################################################################################################################
#Starting from the Cuffdiff-output Directory
#install.packages("ggplot2")
library(ggplot2);library(gplots)
library(edgeR);library(STRINGdb)
library(igraph);source("http://igraph.sf.net")
#If X11 not available, open a pdf device for output of all plots
pdf(file="cummeRbund_output.pdf")
#Clean up workspace - i.e. delete variable created by the graphics demo
#rm(list = ls(all = TRUE))
#List the variables that exist in your current work space
## ----KEGG Download------------------------------------------------------------

library(biomaRt)
Gff2GeneTable("1021_genome.gff3")
load("geneTable.rda")
edb<-geneTable$GeneID
head(geneTable)
######################################################################################################
#### Import the gene expression data from the Tophat/Cufflinks/Cuffdiff
######################################################################################################
#Set working directory where results files exist
working_dir = "~/workspace/rnaseq/de/tophat_cufflinks/ref_only"
setwd(working_dir)
dir()
#Import expression and differential expression results from the Bowtie/Samtools/Tophat/Cufflinks/Cuffdiff pipeline
file1="isoforms.read_group_tracking"
file2="isoform_exp.diff"
file3="isoforms.fpkm_tracking"
file4="cds.read_group_tracking"
file5="cds_exp.diff"
file6="cds.fpkm_tracking"
file7="genes.read_group_tracking"
file8="gene_exp.diff"
file9="genes.fpkm_tracking"
file10="tss_groups.read_group_tracking"
file11="tss_group_exp.diff"
file12="tss_groups.fpkm_tracking"
over="LUTS";under="CTRL"

#Read in tab delimited files and assign the resulting 'dataframe' to a variable
#Use 'as.is' for columns that contain text/character values (i.e. non-numerical values)
all_fpkm = read.table(file1, header=TRUE, sep="\t", as.is=c(1:2,9))
tn_de = read.table(file2, header=TRUE, sep="\t", as.is=c(1:7,14))
tn_fpkm = read.table(file3, header=TRUE, sep="\t", as.is=c(1:9,13,17))
allc_fpkm = read.table(file4, header=TRUE, sep="\t", as.is=c(1:2,9))
tnc_de = read.table(file5, header=TRUE, sep="\t", as.is=c(1:7,14))
tnc_fpkm = read.table(file6, header=TRUE, sep="\t", as.is=c(1:9,13,17))
allg_fpkm = read.table(file7, header=TRUE, sep="\t", as.is=c(1:2,9))
tng_de = read.table(file8, header=TRUE, sep="\t", as.is=c(1:7,14))
tng_fpkm = read.table(file9, header=TRUE, sep="\t", as.is=c(1:9,13,17))
allt_fpkm = read.table(file10, header=TRUE, sep="\t", as.is=c(1:2,9))
tnt_de = read.table(file11, header=TRUE, sep="\t", as.is=c(1:7,14))
tnt_fpkm = read.table(file12, header=TRUE, sep="\t", as.is=c(1:9,13,17))

#View the column names
names(tn_de)

#Determine the dimensions of the dataframe.  'dim()' will return the number of rows and columns
dim(tn_de)
head(tn_de)
#Get the first 3 rows of data and a selection of columns
tn_de[1:3,c(2:4,7,10,12)]

#Do the same thing, but using the column names instead of numbers
tn_de[1:3, c("gene_id","locus","value_1","value_2")]

head(tn_fpkm)
head(tng_fpkm)
head(tnc_fpkm)
head(tnt_fpkm)

#Get ID to gene name mapping
isoform_mapping=tn_fpkm[,"gene_short_name"]
names(isoform_mapping)=tn_fpkm[,"tracking_id"]

gene_mapping=tng_fpkm[,"gene_short_name"]
names(gene_mapping)=tng_fpkm[,"tracking_id"]

cds_mapping=tnc_fpkm[,"gene_short_name"]
names(cds_mapping)=tnc_fpkm[,"tracking_id"]

tss_mapping=tnt_fpkm[,"gene_short_name"]
names(tss_mapping)=tnt_fpkm[,"tracking_id"]


#Reformat per-replicate gene FPKM data into a standard matrix
LUTS_1=allg_fpkm[allg_fpkm[,"condition"]==over & allg_fpkm[,"replicate"]==0,"FPKM"]
LUTS_2=allg_fpkm[allg_fpkm[,"condition"]==over & allg_fpkm[,"replicate"]==1,"FPKM"]
LUTS_3=allg_fpkm[allg_fpkm[,"condition"]==over & allg_fpkm[,"replicate"]==2,"FPKM"]
LUTS_4=allg_fpkm[allg_fpkm[,"condition"]==over & allg_fpkm[,"replicate"]==3,"FPKM"]
LUTS_5=allg_fpkm[allg_fpkm[,"condition"]==over & allg_fpkm[,"replicate"]==4,"FPKM"]
LUTS_6=allg_fpkm[allg_fpkm[,"condition"]==over & allg_fpkm[,"replicate"]==5,"FPKM"]
LUTS_7=allg_fpkm[allg_fpkm[,"condition"]==over & allg_fpkm[,"replicate"]==6,"FPKM"]
LUTS_8=allg_fpkm[allg_fpkm[,"condition"]==over & allg_fpkm[,"replicate"]==7,"FPKM"]

CTRL_1=allg_fpkm[allg_fpkm[,"condition"]==under & allg_fpkm[,"replicate"]==0,"FPKM"]
CTRL_2=allg_fpkm[allg_fpkm[,"condition"]==under & allg_fpkm[,"replicate"]==1,"FPKM"]
CTRL_3=allg_fpkm[allg_fpkm[,"condition"]==under & allg_fpkm[,"replicate"]==2,"FPKM"]
CTRL_4=allg_fpkm[allg_fpkm[,"condition"]==under & allg_fpkm[,"replicate"]==3,"FPKM"]
CTRL_5=allg_fpkm[allg_fpkm[,"condition"]==under & allg_fpkm[,"replicate"]==4,"FPKM"]
CTRL_6=allg_fpkm[allg_fpkm[,"condition"]==under & allg_fpkm[,"replicate"]==5,"FPKM"]
CTRL_7=allg_fpkm[allg_fpkm[,"condition"]==under & allg_fpkm[,"replicate"]==6,"FPKM"]
CTRL_8=allg_fpkm[allg_fpkm[,"condition"]==under & allg_fpkm[,"replicate"]==7,"FPKM"]

ids=unique(allg_fpkm[,"tracking_id"])
gene_names=gene_mapping[ids]
gene.ids=unique(allg_fpkm[,"tracking_id"])
gene_names=gene_mapping[gene.ids]
gene.fpkm_matrix=data.frame(gene.ids, LUTS_1,LUTS_2,LUTS_3,LUTS_4,LUTS_5,LUTS_6,LUTS_7,LUTS_8,CTRL_1,CTRL_2,CTRL_3,CTRL_4,CTRL_5,CTRL_6,CTRL_7,CTRL_8)
row.names(gene.fpkm_matrix)=gene.ids
data_columns=c(2:17)
short_names=c("LUTS_1","LUTS_2","LUTS_3","LUTS_4","LUTS_5","LUTS_6","LUTS_7","LUTS_8","CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5","CTRL_6","CTRL_7","CTRL_8")

head(gene.fpkm_matrix);dim(gene.fpkm_matrix)
head(tng_de);dim(tng_de)

sig_genes_exp.diff<-subset(tng_de, tng_de$significant == "yes")
head(sig_genes_exp.diff);dim(sig_genes_exp.diff)

sig.genes_exp.diff<-as.data.frame(cbind(logFC=sig_genes_exp.diff$log2.fold_change,q_value=sig_genes_exp.diff$q_value),
                                  row.names=sig_genes_exp.diff$gene_id)

head(sig.genes_exp.diff);dim(sig.genes_exp.diff)
gene.fpkm_matrix<-gene.fpkm_matrix[,-1]
head(gene.fpkm_matrix)

head(gene.fpkm_matrix[unique(row.names(gene.fpkm_matrix)),])
g.rep.ma<-gene.fpkm_matrix[unique(row.names(gene.fpkm_matrix)),]
head(g.rep.ma);dim(g.rep.ma)
g.rep.ma<-as.matrix(g.rep.ma)

over.group<-grep(pattern = over, x = colnames(g.rep.ma),ignore.case = T)
under.group<-grep(pattern = under, x = colnames(g.rep.ma),ignore.case = T)

g.under.matrix<-g.rep.ma[,under.group]
head(g.under.matrix)
dim(g.under.matrix)

g.over.matrix<-g.rep.ma[,over.group]
head(g.over.matrix)
dim(g.over.matrix)

sig_over_gene_data<-subset(sig_genes_exp.diff, (significant =='yes') & (log2.fold_change. > 0))
head(sig_over_gene_data);dim(sig_over_gene_data)

sig_under_gene_data<-subset(sig_genes_exp.diff, (significant =='yes') & (log2.fold_change. < 0))
head(sig_under_gene_data);dim(sig_under_gene_data)

sig.genes<-row.names(g.rep.ma) %in% row.names(sig.genes_exp.diff)
head(g.rep.ma[sig.genes,]);dim(g.rep.ma[sig.genes,])
s.g.rep.matrix<-g.rep.ma[sig.genes,]
head(s.g.rep.matrix);dim(s.g.rep.matrix)

min.samps<-(length(colnames(s.g.rep.matrix))/2)*.5
min.expr<-round(min.samps)

isexpr<-rowSums(cpm(s.g.rep.matrix) > 1) >= min.expr
s.g.rep.ma<-s.g.rep.matrix[isexpr,]
dim(s.g.rep.ma)

g.under.genes<-row.names(g.rep.ma) %in% sig_under_gene_data$gene_id
g.u.rep.ma<-g.rep.ma[g.under.genes,]
head(g.u.rep.ma);dim(g.u.rep.ma)

g.over.genes<-row.names(g.rep.ma) %in% sig_over_gene_data$gene_id
g.o.rep.ma<-g.rep.ma[g.over.genes,]
head(g.o.rep.ma);dim(g.o.rep.ma)

s.g.over.matrix<-g.o.rep.ma[,over.group]
head(s.g.over.matrix);dim(s.g.over.matrix)
min.o.samps<-(length(colnames(s.g.over.matrix))*.5)
min.o.expr<-round(min.o.samps)

isexpr<-rowSums(cpm(s.g.over.matrix) > 1) >= min.o.expr
s.g.over.matrix<-s.g.over.matrix[isexpr,]
head(s.g.over.matrix);dim(s.g.over.matrix)

s.g.under.matrix<-g.u.rep.ma[,under.group]
head(s.g.under.matrix);dim(s.g.under.matrix)
min.u.samps<-(length(colnames(s.g.under.matrix))*.5)
min.u.expr<-round(min.u.samps)

isexpr<-rowSums(cpm(s.g.under.matrix) > 1) >= min.u.expr
s.g.under.matrix<-s.g.under.matrix[isexpr,]
head(s.g.under.matrix);dim(s.g.under.matrix)

#under.inf = which(foldchange == "-Inf" & qval < 0.05 | foldchange < -1 & foldchange != "-Inf" & qval < 0.05)
#over.inf = which(foldchange == "Inf" & qval < 0.05 | foldchange > 1 &  foldchange != "Inf" & qval < 0.05)

############################################################################
# adj.matrix.graph-all-over-genes ####
############################################################################

s.g.rep.matrix;s.g.under.matrix;s.g.over.matrix
#Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
over.graph <- graph.adjacency(as.matrix(as.dist(cor(t(s.g.over.matrix), method="pearson"))),
                              mode="undirected", weighted=TRUE, diag=FALSE)
#Simplfy the adjacency object
#over.graph <- simplify(over.graph, remove.multiple=TRUE, remove.loops=TRUE)
#Colour negative correlation edges as blue
E(over.graph)[which(E(over.graph)$weight<0)]$color <- "darkblue"
#Colour positive correlation edges as red
E(over.graph)[which(E(over.graph)$weight>0)]$color <- "darkred"
#Convert edge weights to absolute values
E(over.graph)$weight <- abs(E(over.graph)$weight)
#Change arrow size #For directed graphs only
E(over.graph)$arrow.size <- 1.0
#Remove edges below absolute Pearson correlation 0.7
over.graph <- delete_edges(over.graph, E(over.graph)[which(E(over.graph)$weight<0.7)])
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
vSizes <- (scale01(apply(s.g.over.matrix, 1, mean)) + 1.0) * 10
#Amplify or decrease the width of the edges
edgeweights <- E(over.graph)$weight * 2.0
#Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst <- mst(over.graph, algorithm="prim")
#Plot the tree object
plot(mst, layout=layout.fruchterman.reingold, edge.curved=TRUE,
     vertex.size=vSizes, vertex.label.dist=-0.5, vertex.label.color="black",
     asp=FALSE, vertex.label.cex=0.6, edge.width=edgeweights, edge.arrow.mode=0,
     main="Network of significant genes in Over group")
# code chunk number 4: Identify communities in the tree object based on 'edge betweenness'
#mst.communities <- edge.betweenness.community(mst, directed=T)
mst.communities <- edge.betweenness.community(mst, directed=T)
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
  main="Over Group Significant Gene Network"
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
  main="Over Group Significant Gene Network"
)

############################################################################
# adj.matrix.graph-all-under-genes ####
############################################################################

#under.graph <- graph.adjacency(as.matrix(as.dist(cor(t(s.g.under.matrix), method="pearson"))),
#                               mode="directed", weighted=TRUE, diag=FALSE)
#Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
under.graph <- graph.adjacency(as.matrix(as.dist(cor(t(s.g.under.matrix), method="pearson"))),
                               mode="undirected", weighted=TRUE, diag=FALSE)
#NB - Euclidean distances can also be used instead of correlation, but this changes the tutorial slightly:
#g <- graph.adjacency(as.matrix(dist(estrogenMainEffects, method="euclidean")),
#mode="undirected", weighted=TRUE, diag=FALSE)
#Simplfy the adjacency object
#under.graph <- simplify(under.graph, remove.multiple=TRUE, remove.loops=TRUE)
#Colour negative correlation edges as blue
E(under.graph)[which(E(under.graph)$weight<0)]$color <- "darkblue"
#Colour positive correlation edges as red
E(under.graph)[which(E(under.graph)$weight>0)]$color <- "darkred"
#Convert edge weights to absolute values
E(under.graph)$weight <- abs(E(under.graph)$weight)
#Change arrow size #For directed graphs only
E(over.graph)$arrow.size <- 1.0
#Remove edges below absolute Pearson correlation 0.7
under.graph <- delete_edges(under.graph, E(under.graph)[which(E(under.graph)$weight<0.7)])
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
vSizes <- (scale01(apply(s.g.under.matrix, 1, mean)) + 1.0) * 10
#Amplify or decrease the width of the edges
edgeweights <- E(under.graph)$weight * 2.0
#Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst <- mst(under.graph, algorithm="prim")
#Plot the tree object
plot(mst, layout=layout.fruchterman.reingold, edge.curved=TRUE,
     vertex.size=vSizes, vertex.label.dist=-0.5, vertex.label.color="black",
     asp=FALSE, vertex.label.cex=0.6, edge.width=edgeweights, edge.arrow.mode=0,
     main="Network of Significant Genes in Under Group")
## code chunk number 4: Identify communities in the tree object based on 'edge betweenness'
#mst.communities <- edge.betweenness.community(mst, directed=T)
mst.communities <- edge.betweenness.community(mst, directed=F)
mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
#mst.m <- make_clusters(mst, membership=mst.communities$merges)
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
  main="Under Group Significant Gene Network"
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
  main="Under Group Significant Gene Network"
)

#Check the vertex degree, i.e., number of connections to each vertex
degree(mst)
#Output information for each community, including vertex-to-community assignments and modularity
commSummary <- data.frame(mst.communities$names,mst.communities$membership, mst.communities$modularity)
commSummary <- data.frame(mst.communities$names,mst.communities$membership, mst.clustering$modularity)

colnames(commSummary) <- c("Gene", "Community", "Modularity")
options(scipen=999);commSummary
#Compare community structures using the variance of information (vi) metric (not relevant here)
#Community structures that are identical will have a vi=0
compare(mst.communities, mst.communities, method="vi")

## ----selectOrg2, eval=FALSE----------------------------------------------

g.rep.cor<-cor(t(s.g.under.matrix))
diag(g.rep.cor)<-0
g1<-graph_from_adjacency_matrix(g.rep.cor, weighted = T)
g2<-delete.edges(g1, E(g1)[abs(weight) < 0.8])
write_graph(g1, "significant_gene_FPKM_sub_network.gml", format="gml")

Visualize network in Cytoscape
Import gml file (Import Network From File)
Import annotation file (Import Table from )
Analyze network

###################################################
# STRING.db-Pathway_Analysis ####
###################################################
#source("https://bioconductor.org/biocLite.R")
#biocLite(c("Rsubread","GSAR","STRINGdb"))
library(STRINGdb)
#browseVignettes("STRINGdb")
STRINGdb$help("get_bioc_graph"); STRINGdb$help("get_graph")
#hsa.string.db <- STRINGdb$new( version="10", species=9606,
#                               score_threshold=0, input_directory="" )
#hsa.string.db
# species is a variable that is defined in the settings file
species="Homo sapiens"
hsa.kegg <- search_kegg_organism(species, by='scientific_name')
dim(hsa.kegg); head(hsa.kegg)
hsa.kegg$kegg_code

species.all<-get_STRING_species(version="10", species_name=NULL)
colnames(species.all);head(species.all);dim(species.all)
hsa9606<-grep(pattern=species, species.all$official_name, ignore.case = T)
hsa.taxa.info<-species.all[hsa9606,]
hsa.taxID<-hsa.taxa.info$species_id
string.db.hsa9606 <- STRINGdb$new(version="10", species=hsa.taxID,
                                  score_threshold=0, input_directory="")
string.db.hsa9606
## get_interactions(string_ids)   # returns the interactions in between the input proteins
## get_neighbors(string_ids)      # Get the neighborhoods of a protein (or of a vector of proteins).
## get_subnetwork(string_ids)     # returns a subgraph from the given input proteins

hsa.pwys <- download.kegg.pathways(hsa.kegg$kegg_code)
hsa.kegg.gs <- getGenesets(hsa.kegg$kegg_code)
hsa.kegg.sets<-kegg.gsets(species = hsa.kegg$kegg_code,id.type = "kegg")
hsa.9606.kegg<-download_KEGG(species=hsa.kegg$kegg_code, keggType = "KEGG", keyType = "kegg")
head(hsa.kegg.sets)
head(hsa.pwys)
head(hsa.kegg.gs)
head(hsa.9606.kegg)

hsa.9606.kegg$KEGGPATHID2EXTID[1:10,1]
hsa.9606.kegg$KEGGPATHID2EXTID[1:10,2]
length(hsa.9606.kegg)
names(hsa.9606.kegg)
na.seq.sum<-runInfo(cuff)
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
refgtf<-file.path("/home/drew/umb_triley/RefGenomes/UCSC_hg38/genes.gtf")

#####################################################################
# enrichment
######################################################################

sig.genes.exp.df<-as.data.frame(cbind(pvalue=sig_genes_exp.diff$q_value,
                                      logFC=sig_genes_exp.diff$log2.fold_change,
                                      gene=sig_genes_exp.diff$gene_id), stringsAsFactors=F)

head(sig.genes.exp.df);dim(sig.genes.exp.df)
sig.genes.DE_mapped <- string.db.hsa9606$map( sig.genes.exp.df, "gene", removeUnmappedRows = TRUE )
#write.table(sig.genes.DE_mapped, file="AB.vs.wt1021B.KEGG.difftable")
head(sig.genes.DE_mapped);dim(sig.genes.DE_mapped)

sig.genes.DE.df<-as.data.frame(cbind(gene=sig.genes.DE_mapped$gene,
                                         pvalue=sig.genes.DE_mapped$pvalue,
                                         logFC=sig.genes.DE_mapped$logFC), stringsAsFactors=F)
dim(sig.genes.DE.df)
head(sig.genes.DE.df)

sig.genes.intersected<-string.db.hsa9606$map(sig.genes.DE.df, "gene", removeUnmappedRows=T)
head(sig.genes.intersected)
dim(sig.genes.intersected)
string.db.hsa9606$plot_network(sig.genes.intersected$STRING_id[1:100],)
sig.gene.subnets.str<-string.db.hsa9606$get_subnetwork(sig.genes.intersected$STRING_id[1:100])
sig.gene.subnets.str
sig.gene.mapped_sig<-as.data.frame(cbind(genes=c(sig.genes.intersected$gene[sig.genes.intersected$pvalue < 0.05]),
                                              pvalue=c(sig.genes.intersected$pvalue[sig.genes.intersected$pvalue < 0.05]),
                                              logFC=c(sig.genes.intersected$logFC[sig.genes.intersected$pvalue < 0.05]),
                                              STRING_id=c(sig.genes.intersected$STRING_id[sig.genes.intersected$pvalue < 0.05])),
                                        stringsAsFactors=F, row.names=F)
head(sig.gene.mapped_sig)
dim(sig.gene.mapped_sig)

sig.gene.DE.pv.fc.STRING<-as.data.frame(cbind(gene=sig.genes.DE_mapped$gene, pvalue=sig.genes.DE_mapped$pvalue,
                                                   logFC=sig.genes.DE_mapped$logFC, STRING_id=sig.genes.DE_mapped$STRING_id),
                                                    stringsAsFactors=F, row.names=F, col.names=T)
head(sig.gene.DE.pv.fc.STRING)
dim(sig.gene.DE.pv.fc.STRING)
# post payload information to the STRING server
sig.gene_pval01 <- string.db.hsa9606$post_payload(sig.gene.mapped_sig$STRING_id,colors=sig.gene.mapped_sig["pvalue"]$color )
# display a STRING network png with the "halo"
string.db.hsa9606$plot_network(sig.gene.DE.pv.fc.STRING$STRING_id[1:50],
                                payload_id=sig.gene_pval01,
                                required_score=sig.gene.DE.pv.fc.STRING$logFC[1:50])


sig.under.genes<-sig_genes_exp.diff$gene %in% row.names(s.g.under.matrix)
sig.under.exp.df<-sig_genes_exp.diff[sig.under.genes,]
head(sig.under.exp.df);dim(sig.under.exp.df)
sig.u.exp.df<-as.data.frame(cbind(gene=sig.under.exp.df$gene_id,
                                      logFC=sig.under.exp.df$log2.fold_change,
                                      pvalue=sig.under.exp.df$q_value), stringsAsFactors=F)
head(sig.u.exp.df);dim(sig.u.exp.df)

sig.over.genes<-sig_genes_exp.diff$gene %in% row.names(s.g.over.matrix)
sig.over.exp.df<-sig_genes_exp.diff[sig.over.genes,]
head(sig.over.exp.df);dim(sig.over.exp.df)
sig.o.exp.df<-as.data.frame(cbind(gene=sig.over.exp.df$gene_id,
                                      logFC=sig.over.exp.df$log2.fold_change,
                                      pvalue=sig.over.exp.df$q_value), stringsAsFactors=F)
head(sig.o.exp.df);dim(sig.o.exp.df)

sig.over.DE_mapped <- string.db.hsa9606$map( sig.o.exp.df, "gene", removeUnmappedRows = TRUE )
write.table(sig.over.DE_mapped, file="sig.over.string.difftable")
head(sig.over.DE_mapped);dim(sig.over.DE_mapped)
string.db.hsa9606$plot_network(sig.over.DE_mapped$STRING_id[1:400],)

sig.under.DE_mapped <- string.db.hsa9606$map( sig.u.exp.df, "gene", removeUnmappedRows = TRUE )
write.table(sig.under.DE_mapped, file="sig.over.string.difftable")
head(sig.under.DE_mapped);dim(sig.under.DE_mapped)
string.db.hsa9606$plot_network(sig.under.DE_mapped$STRING_id[1:400],)

#####################################################################
# enrichment
#######################################################################

sig.over.DE.df<-as.data.frame(cbind(gene=sig.over.DE_mapped$gene,
                                         pvalue=sig.over.DE_mapped$pvalue,
                                         logFC=sig.over.DE_mapped$logFC), stringsAsFactors=F)
dim(sig.over.DE.df)
head(sig.over.DE.df)

sig.over.intersected<-string.db.hsa9606$map(sig.over.DE.df, "gene", removeUnmappedRows=T)
head(sig.over.intersected);dim(sig.over.intersected)
class(sig.over.intersected)
string.db.hsa9606$plot_network(sig.over.intersected$STRING_id[1:100],)

sig.over.subnets<-string.db.hsa9606$get_subnetwork(sig.over.intersected$STRING_id)
sig.over.subnets

sig.over.mapped_sig<-as.data.frame(cbind(genes=c(sig.over.intersected$gene[sig.over.intersected$pvalue < 0.05]),
                                              pvalue=c(sig.over.intersected$pvalue[sig.over.intersected$pvalue < 0.05]),
                                              logFC=c(sig.over.intersected$logFC[sig.over.intersected$pvalue < 0.05]),
                                              STRING_id=c(sig.over.intersected$STRING_id[sig.over.intersected$pvalue < 0.05])),
                                        stringsAsFactors=F, row.names=F)
head(sig.over.mapped_sig)

sig.over.DE.pv.fc.STRING<-as.data.frame(cbind(gene=sig.over.DE_mapped$gene,
                                                   pvalue=sig.over.DE_mapped$pvalue,
                                                   logFC=sig.over.DE_mapped$logFC,
                                                   STRING_id=sig.over.DE_mapped$STRING_id), stringsAsFactors=F, row.names=F, col.names=T)
head(sig.over.DE.pv.fc.STRING)
# post payload information to the STRING server
sig.over.payload<-string.db.hsa9606$post_payload(sig.over.mapped_sig$STRING_id,
                                                       colors=sig.over.mapped_sig["logFC"]$color, )

# display a STRING network png with the "halo"
string.db.hsa9606$plot_network( sig.over.DE.pv.fc.STRING$STRING_id[1:50],
                                payload_id=sig.over.payload,
                                required_score=sig.over.DE.pv.fc.STRING$logFC[1:50])

string.db.hsa9606$plot_network( sig.over.DE.pv.fc.STRING$STRING_id[1:100],
                                payload_id=sig.over.payload,
                                required_score=sig.over.DE.pv.fc.STRING$logFC[1:100])

# plot the enrichment for the best 100 genes
sig.over.top100<-string.db.hsa9606$plot_ppi_enrichment( sig.over.DE_mapped$STRING_id[1:100], quiet=F,)
hits<-sig.over.DE_mapped$STRING_id[1:200]
s.o.backgroundV <- sig.over.DE_mapped$STRING_id  # as an example, we use the first 2000 genes
string.db.hsa9606$set_background(s.o.backgroundV)
string.db.bg <- STRINGdb$new( version="10", score_threshold=0, backgroundV = s.o.backgroundV, species=9606)
eh <- string.db.hsa9606$enrichment_heatmap(s.o.enrichmentGOmf, title="Significantly upregulated gene networks for each condition" )
eh <- string.db.hsa9606$enrichment_heatmap(list( hits[1:100], hits[101:200]),
                                           list("over-group 1st tier","under-group 2nd tier"), title="Significantly upregulated gene networks for each condition" )
clustersList <- string.db.hsa9606$get_clusters(sig.over.DE_mapped$STRING_id[1:200])

for(i in seq(1:5)){
  string.db.hsa9606$plot_network(clustersList[[i]])
}


#####################################################################
# enrichment
#######################################################################

s.o.enrichmentGObp <- string.db.hsa9606$get_enrichment( hits, category = "Process", methodMT = "fdr", iea = TRUE )
s.o.enrichmentGOmf <- string.db.hsa9606$get_enrichment( hits, category = "Function", methodMT = "fdr", iea = TRUE )
s.o.enrichmentGOcc <- string.db.hsa9606$get_enrichment( hits, category = "Component", methodMT = "fdr", iea = TRUE )
s.o.enrichmentKEGG <- string.db.hsa9606$get_enrichment( hits, category = "KEGG", methodMT = "fdr", iea = TRUE )
head(s.o.enrichmentGObp, n=7); dim(s.o.enrichmentGObp)
head(s.o.enrichmentGOmf, n=7); dim(s.o.enrichmentGOmf)
head(s.o.enrichmentGOcc, n=7); dim(s.o.enrichmentGOcc)
head(s.o.enrichmentKEGG, n=7); dim(s.o.enrichmentKEGG)

## ----fig.height=12, fig.width=8------------------------------------------
sig.over.mkegg <- enrichMKEGG(gene = sig.over.DE_mapped$gene, organism = 'hsa')
gene.df <- bitr_kegg(sig.over.DE_mapped$gene, fromType = "SYMBOL",toType = c("ENTREZID", "KEGG"),organism=species, drop=T)

AB.vs.wt1021.kegg <- enrichKEGG(gene = AB.vs.wt1021.DE_mapped$GeneSymbol,
                                organism = 'sme')

barplot(AB.vs.wt1021.mkegg, drop=TRUE, showCategory=12)
barplot(AB.vs.wt1021.kegg, showCategory=8)
dotplot(AB.vs.wt1021.mkegg)
dotplot(AB.vs.wt1021.kegg)
cnetplot(AB.vs.wt1021.mkegg, categorySize="pvalue") # ,wt1021.vs.wt1021B.kegg
enrichMap(AB.vs.wt1021.mkegg)
cnetplot(AB.vs.wt1021.kegg, categorySize="pvalue",foldChange=AB.vs.wt1021.DE_mapped$LogFoldChange)
enrichMap(AB.vs.wt1021.kegg)
# GO analysis adjusting for gene length bias
# (assuming that y$genes$Length contains gene lengths)
library(EnrichmentBrowser)
go.abund <- goanna(sig.over.DE.df, geneid = "gene", trend = T)
go.abund

go.len <- goanna(A.vs.wt1021.DE, geneid = "GeneID", trend = "Length")

topGO(s.o.enrichmentGObp) #, sort =s.o.enrichmentGObp$pvalue)
topGO(go.len, sort = "Pval")
#Avswt1021.up_reg_genes<-topGO(go.abund, sort = "Qval")
#tAvswt1021.down_reg_genes<-topGO(go.abund, sort = "Qval")


## Default usage with a list of gene sets:

go.de <- goana(list(DE1 = EG.DE1, DE2 = EG.DE2, DE3 = EG.DE3))
topGO(go.de, sort = "DE1")
topGO(go.de, sort = "DE2")
topGO(go.abund, ontology = "BP")
topGO(go.de, ontology = "CC", sort = "DE3")
topGO(go.de, ontology = "MF", sort = "DE3")

## Standard KEGG analysis

sig.over.DE.kegg <- kegga(sig.over.DE.df$gene, species.KEGG="hsa") # equivalent to previous
sig.over.DE.kegg

barplot(sig.over.DE.kegg$DE,drop=T, showCategory=8)
## ------------------------------------------------------------------------
dotplot(AB.vs.wt1021.DE.kegg)
## ----fig.cap="enrichment map of enrichment result", fig.align="center", fig.height=16, fig.width=16, eval=FALSE----
enrichMap(mkk)
## ## categorySize can be scaled by 'pvalue' or 'geneNum'
cnetplot(sig.over.DE_mapped$STRING_id, categorySize="pvalue", foldChange=sig.over.DE_mapped$logFC)

## ----fig.height=12, fig.width=8------------------------------------------
plotGOgraph(s.o.enrichmentGOmf)
## ----fig.height=5, fig.width=8-------------------------------------------
barplot(ego, showCategory=8)

## ------------------------------------------------------------------------
dotplot(mkk)

## ----
#fig.cap="enrichment map of enrichment result", fig.align="center", fig.height=16, fig.width=16, eval=FALSE----
# enrichMap(ego)

## ----fig.height=14, fig.width=14, eval=FALSE-----------------------------
## ## categorySize can be scaled by 'pvalue' or 'geneNum'
## cnetplot(ego, categorySize="pvalue", foldChange=geneList)

## ----fig.height=12, fig.width=8------------------------------------------
plotGOgraph(ego)

## ----fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8----
gseaplot(kk2, geneSetID = "hsa04145")

## ----eval=FALSE----------------------------------------------------------
## browseKEGG(kk, 'hsa04110')
## ----fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8----
gseaplot(AB.vs.wt1021.DE.kegg, geneSetID = "sme")

## ------------------------------------------------------------------------
ego <- enrichGO(gene=entrezgenes,
                universe=names(geneList),
                OrgDb= org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

## ----eval=FALSE----------------------------------------------------------
## ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
##                 OrgDb         = org.Hs.eg.db,
## 		keytype       = 'ENSEMBL',
##                 ont           = "CC",
##                 pAdjustMethod = "BH",
##                 pvalueCutoff  = 0.01,
##                 qvalueCutoff  = 0.05)

## ----eval=FALSE----------------------------------------------------------
## ego2 <- setReadable(ego2, OrgDb = org.Hs.eg.db)

## ----eval=FALSE----------------------------------------------------------
## ego3 <- gseGO(geneList     = geneList,
##               OrgDb        = org.Hs.eg.db,
##               ont          = "CC",
##               nPerm        = 1000,
##               minGSSize    = 100,
##               maxGSSize    = 500,
##               pvalueCutoff = 0.05,
##               verbose      = FALSE)


## ----KEGG Download------------------------------------------------------------
sme.kegg.code<-search_kegg_organism('sme', by='kegg_code')
sme.kegg.code
go.abund
Smeliloti.kegg <- search_kegg_organism('Sinorhizobium meliloti 1021', by='scientific_name')
dim(Smeliloti.kegg)
head(Smeliloti.kegg)
sme.1021.kegg<-download_KEGG(species="sme", keggType = "KEGG", keyType = "kegg")
sme.1021.kegg$KEGGPATHID2EXTID[1:10,1]
sme.1021.kegg$KEGGPATHID2EXTID[1:10,2]
length(sme.1021.kegg)
names(sme.1021.kegg)
bitr_kegg
sme.1021.kegg
## ------------------------------------------------------------------------
barcodeplot(AB.vs.wt1021B.DE[,8], index = AB.vs.wt1021B.DE[,7],index2 = AB.vs.wt1021B.DE[,8], col.bars = "dodgerblue",alpha=.01,
            labels = "LogFoldChange",xlab="FoldChange")

barcodeplot(A.vs.wt1021.DE[,8], index = A.vs.wt1021.DE[,7],index2 = A.vs.wt1021.DE[,8], col.bars = "dodgerblue",alpha=.01,
            labels = "LogFoldChange",xlab="FoldChange")

barcodeplot(wt1021B.vs.wt1021.DE[,8], index = wt1021B.vs.wt1021.DE[,7],index2 = wt1021B.vs.wt1021.DE[,8], col.bars = "dodgerblue",alpha=.01,
            labels = "LogFoldChange",xlab="FoldChange")

barcodeplot(A.vs.AB.DE[,8], index = A.vs.AB.DE[,7],index2 = A.vs.AB.DE[,8], col.bars = "dodgerblue",alpha=.01,
            labels = "LogFoldChange",xlab="FoldChange")

barcodeplot(AB.vs.wt1021.DE[,8], index = AB.vs.wt1021.DE[,7],index2 = AB.vs.wt1021.DE[,8], col.bars = "dodgerblue",alpha=.01,
            labels = "LogFoldChange",xlab="FoldChange")

## ------------------------------------------------------------------------
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)

## ----eval = FALSE--------------------------------------------------------
sme.genes <- enrichMKEGG(gene = geneList,
                         organism = 'sme')

## ----eval=FALSE----------------------------------------------------------
## mkk2 <- gseMKEGG(geneList = geneList,
##                  species = 'hsa')

## ----eval=FALSE----------------------------------------------------------
## david <- enrichDAVID(gene = gene,
##                      idType = "ENTREZ_GENE_ID",
##                      listType = "Gene",
##                      annotation = "KEGG_PATHWAY",
##                      david.user = "clusterProfiler@hku.hk")

## ------------------------------------------------------------------------
gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)

egmt <- enricher(gene, TERM2GENE=c5)
head(egmt)

egmt2 <- GSEA(geneList, TERM2GENE=c5, verbose=FALSE)
head(egmt2)

## ----fig.height=5, fig.width=9-------------------------------------------
barplot(mkk, drop=TRUE, showCategory=12)

## ----fig.height=5, fig.width=8-------------------------------------------
barplot(ego, showCategory=8)

## ------------------------------------------------------------------------
dotplot(mkk)

## ----
#fig.cap="enrichment map of enrichment result", fig.align="center", fig.height=16, fig.width=16, eval=FALSE----
# enrichMap(ego)

## ----fig.height=14, fig.width=14, eval=FALSE-----------------------------
## ## categorySize can be scaled by 'pvalue' or 'geneNum'
## cnetplot(ego, categorySize="pvalue", foldChange=geneList)

## ----fig.height=12, fig.width=8------------------------------------------
plotGOgraph(ego)

## ----fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8----
gseaplot(kk2, geneSetID = "hsa04145")

#####################################################################
# enrichment
###################################################

sig.under.DE.df<-as.data.frame(cbind(gene=sig.under.DE_mapped$gene,
                                    pvalue=sig.under.DE_mapped$pvalue,
                                    logFC=sig.under.DE_mapped$logFC), stringsAsFactors=F)
dim(sig.under.DE.df)
head(sig.under.DE.df)

sig.under.intersected<-string.db.hsa9606$map(sig.under.DE.df, "gene", removeUnmappedRows=T)
head(sig.under.intersected);dim(sig.under.intersected)

string.db.hsa9606$plot_network(sig.under.intersected$STRING_id[1:100],)

sig.under.subnets<-string.db.hsa9606$get_subnetwork(sig.under.intersected)
sig.under.subnets

sig.under.mapped_sig<-as.data.frame(cbind(genes=c(sig.under.intersected$gene[sig.under.intersected$pvalue < 0.05]),
                                         pvalue=c(sig.under.intersected$pvalue[sig.under.intersected$pvalue < 0.05]),
                                         logFC=c(sig.under.intersected$logFC[sig.under.intersected$pvalue < 0.05]),
                                         STRING_id=c(sig.under.intersected$STRING_id[sig.under.intersected$pvalue < 0.05])),
                                   stringsAsFactors=F, row.names=F)
dim(sig.under.mapped_sig)

sig.under.DE.pv.fc.STRING<-as.data.frame(cbind(gene=sig.under.DE_mapped$gene,
                                              pvalue=sig.under.DE_mapped$pvalue,
                                              logFC=sig.under.DE_mapped$logFC,
                                              STRING_id=sig.under.DE_mapped$STRING_id), stringsAsFactors=F, row.names=F, col.names=T)
head(sig.under.DE.pv.fc.STRING)

# post payload information to the STRING server
sig.under.payload<-string.db.hsa9606$post_payload(sig.under.mapped_sig$STRING_id,
                                                 colors=sig.under.mapped_sig["pvalue"]$color, )

# display a STRING network png with the "halo"

string.db.hsa9606$plot_network( sig.under.DE.pv.fc.STRING$STRING_id[1:50],
                                payload_id=sig.under.payload,
                                required_score=sig.under.DE.pv.fc.STRING$logFC[1:50])

string.db.hsa9606$plot_network( sig.under.DE.pv.fc.STRING$STRING_id,
                                payload_id=sig.under.payload,
                                required_score=sig.under.DE.pv.fc.STRING$logFC)

# plot the enrichment for the best 100 genes
sig.under.top100<-string.db.hsa9606$plot_ppi_enrichment( sig.under.intersected$STRING_id[1:100], quiet=F,)

#############################################################################


sig.under.DE_mapped <- string.db.hsa9606$map( sig.u.exp.df, "gene", removeUnmappedRows = TRUE )
write.table(sig.under.DE_mapped, file="sig.under.string.difftable")
head(sig.under.DE_mapped);dim(sig.under.DE_mapped)
string.db.hsa9606$plot_network(sig.under.DE_mapped$STRING_id[1:400],)


AB.vs.wt1021B.subnets<-string.db.sme1021$get_subnetwork(AB.vs.wt1021B.intersected)
AB.vs.wt1021B.subnets

AB.vs.wt1021B.mapped_sig<-as.data.frame(cbind(genes=c(AB.vs.wt1021B.intersected$gene[AB.vs.wt1021B.intersected$pvalue < 0.05]),
                                              pvalue=c(AB.vs.wt1021B.intersected$pvalue[AB.vs.wt1021B.intersected$pvalue < 0.05]),
                                              logFC=c(AB.vs.wt1021B.intersected$logFC[AB.vs.wt1021B.intersected$pvalue < 0.05]),
                                              STRING_id=c(AB.vs.wt1021B.intersected$STRING_id[AB.vs.wt1021B.intersected$pvalue < 0.05])),
                                        stringsAsFactors=F, row.names=F)
head(AB.vs.wt1021B.mapped_sig)

sig.under.DE.pv.fc.STRING<-as.data.frame(cbind(gene=sig.under.DE_mapped$gene,
                                                   pvalue=sig.under.DE_mapped$pvalue,
                                                   logFC=sig.under.DE_mapped$logFC,
                                                   STRING_id=sig.under.DE_mapped$STRING_id),
                                             stringsAsFactors=F, row.names=F, col.names=T)
head(sig.under.DE.pv.fc.STRING)

# post payload information to the STRING server
sig.under.payload <- string.db.hsa9606$post_payload(sig.under.DE_mapped$STRING_id,
                                                    colors=sig.under.DE_mapped["pvalue"]$color )

# display a STRING network png with the "halo"

string.db.hsa9606$plot_network( sig.under.DE_mapped$STRING_id[1:50],
                                payload_id=sig.under.payload,
                                required_score=sig.under.DE_mapped$logFC[1:50])

# plot the enrichment for the best 100 genes
sig.under.top100<-string.db.hsa9606$plot_ppi_enrichment( sig.under.DE_mapped$STRING_id[1:100], quiet=TRUE )

## ----eval = FALSE--------------------------------------------------------
all.sig.kegg.rich <- enrichKEGG(gene = sig.genes.DE_mapped$gene,organism='hsa',pvalueCutoff = 0.05)
head(all.sig.kegg.rich)
dim(all.sig.kegg.rich)
wt1021.vs.wt1021B.mkegg.rich <- enrichMKEGG(gene = wt1021.vs.wt1021B.DE_mapped$GeneSymbol,organism='sme',pvalueCutoff = 0.05)
head(wt1021.vs.wt1021B.mkegg.rich)
dim(wt1021.vs.wt1021B.mkegg.rich)

barplot(wt1021.vs.wt1021B.mkegg.rich, drop=TRUE, showCategory=12)
barplot(wt1021.vs.wt1021B.kegg.rich,drop=T, showCategory=12)
dotplot(wt1021.vs.wt1021B.mkegg.rich)
dotplot(wt1021.vs.wt1021B.kegg)
cnetplot(wt1021.vs.wt1021B.mkegg.rich, categorySize="pvalue") # ,wt1021.vs.wt1021B.kegg
enrichMap(wt1021.vs.wt1021B.mkegg.rich)
cnetplot(wt1021.vs.wt1021B.kegg.rich, categorySize="pvalue") # ,wt1021.vs.wt1021B.kegg
enrichMap(wt1021.vs.wt1021B.kegg.rich)
cnetplot(wt1021.vs.wt1021B.kegg.rich,categorySize="pvalue", foldChange=,wt1021.vs.wt1021B.DE_mapped$LogFoldChange, )
cnetplot(wt1021.vs.wt1021B.kegg.rich, categorySize="pvalue") # ,wt1021.vs.wt1021B.kegg
enrichMap(wt1021.vs.wt1021B.kegg.rich)


## ----fig.height=12, fig.width=8------------------------------------------
AB.vs.wt1021.mkegg <- enrichMKEGG(gene = AB.vs.wt1021.DE_mapped$GeneSymbol,
                                  organism = 'sme')
AB.vs.wt1021.kegg <- enrichKEGG(gene = AB.vs.wt1021.DE_mapped$GeneSymbol,
                                organism = 'sme')

barplot(AB.vs.wt1021.mkegg, drop=TRUE, showCategory=12)
barplot(AB.vs.wt1021.kegg, showCategory=8)
dotplot(AB.vs.wt1021.mkegg)
dotplot(AB.vs.wt1021.kegg)
cnetplot(AB.vs.wt1021.mkegg, categorySize="pvalue") # ,wt1021.vs.wt1021B.kegg
enrichMap(AB.vs.wt1021.mkegg)
cnetplot(AB.vs.wt1021.kegg, categorySize="pvalue",foldChange=AB.vs.wt1021.DE_mapped$LogFoldChange)
enrichMap(AB.vs.wt1021.kegg)
# GO analysis adjusting for gene length bias
# (assuming that y$genes$Length contains gene lengths)
library(EnrichmentBrowser)
go.abund <- goa(A.vs.wt1021.DE, geneid = "GeneID", trend = T)
go.abund

go.len <- goanna(A.vs.wt1021.DE, geneid = "GeneID", trend = "Length")

topGO(go.len, sort = "Qval")
topGO(go.len, sort = "Pval")
#Avswt1021.up_reg_genes<-topGO(go.abund, sort = "Qval")
#tAvswt1021.down_reg_genes<-topGO(go.abund, sort = "Qval")


## Default usage with a list of gene sets:

go.de <- goana(list(DE1 = EG.DE1, DE2 = EG.DE2, DE3 = EG.DE3))
topGO(go.de, sort = "DE1")
topGO(go.de, sort = "DE2")
topGO(go.abund, ontology = "BP")
topGO(go.de, ontology = "CC", sort = "DE3")
topGO(go.de, ontology = "MF", sort = "DE3")

## Standard KEGG analysis

AB.vs.wt1021.DE.kegg <- kegga(AB.vs.wt1021.DE$GeneSymbol, species.KEGG="sme") # equivalent to previous
AB.vs.wt1021.DE.kegg

barplot(AB.vs.wt1021.DE.kegg$DE,drop=T, showCategory=8)
## ------------------------------------------------------------------------
dotplot(AB.vs.wt1021.DE.kegg)
## ----fig.cap="enrichment map of enrichment result", fig.align="center", fig.height=16, fig.width=16, eval=FALSE----
enrichMap(mkk)
## ## categorySize can be scaled by 'pvalue' or 'geneNum'
cnetplot(AB.vs.wt1021.DE_mapped$STRING_id, categorySize="pvalue", foldChange=AB.vs.wt1021.DE_mapped$LogFoldChange)

## ----fig.height=12, fig.width=8------------------------------------------
plotGOgraph(ego)

## ----fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8----
gseaplot(AB.vs.wt1021.DE.kegg, geneSetID = "sme")

head(ggo)
AB.vs.wt1021.DE_gse<-as.data.frame(AB.vs.wt1021.DE_mapped$Pval, row.names=c(AB.vs.wt1021.DE_mapped$GeneSymbol), stringsAsFactors=F)
names(AB.vs.wt1021.DE_sorted.gse)<-AB.vs.wt1021.DE_mapped$GeneSymbol
AB.vs.wt1021.DE_sorted.gse<-sort(AB.vs.wt1021.DE.kegg$P.DE, decreasing=T)
kk2 <- gseKEGG(geneList     = AB.vs.wt1021.DE_sorted.gse,
               organism     = 'sme',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)
##############################################################################
### DOSE/ClusterProfile - KEGG and Gene Ontology analysis
###############################################################################
wt1021.vs.wt1021B.mkegg.rich <- enrichMKEGG(gene = wt1021.vs.wt1021B.DE_mapped$GeneSymbol,organism = 'sme')
wt1021.vs.wt1021B.kegg.rich <- enrichKEGG(gene = wt1021.vs.wt1021B.DE_mapped$GeneSymbol,organism='sme',pvalueCutoff = 0.05)
head(wt1021.vs.wt1021B.kegg.rich)
dim(wt1021.vs.wt1021B.kegg.rich)
wt1021.vs.wt1021B.DE.kegg.rich <- enrichMKEGG(gene = wt1021.vs.wt1021B.DE_mapped$GeneSymbol,organism='sme')
head(wt1021.vs.wt1021B.DE.kegg.rich)
dim(wt1021.vs.wt1021B.DE.kegg.rich)


AB.vs.wt1021.kegg.rich <- enrichMKEGG(gene = AB.vs.wt1021.DE_mapped$GeneSymbol,organism='sme',pvalueCutoff = 0.05)
head(AB.vs.wt1021.kegg.rich)
dim(AB.vs.wt1021.kegg.rich)
AB.vs.wt1021B.kegg.rich <- enrichKEGG(gene = AB.vs.wt1021B.DE_mapped$GeneSymbol,organism='sme',pvalueCutoff = 0.05)
head(AB.vs.wt1021B.kegg.rich)
dim(AB.vs.wt1021B.kegg.rich)
AB.vs.wt1021B.kegg.rich
A.vs.wt1021.kegg.rich <- enrichKEGG(gene = A.vs.wt1021.DE_mapped$GeneSymbol,organism='sme',pvalueCutoff = 0.05)
head(A.vs.wt1021.kegg.rich)
dim(A.vs.wt1021.kegg.rich)
A.vs.wt1021.kegg.rich

gene.df <- bitr(names(AB.vs.wt1021.DE_mapped$GeneSymbol), fromType = "SYMBOL",toType = c("ENTREZID", "KEGG"),OrgDb=sme.1021.kegg)


library("EnrichmentBrowser")
vignette("EnrichmentBrowser")
sme.pwys <- download.kegg.pathways("sme")
kegg.gs <- get.kegg.genesets("sme")
head(sme.pwys)
head(kegg.gs)

library(gage)
data(gse16873)
sme.kegg.sets<-kegg.gsets(species = "sme",id.type = "kegg")
sme.kegg.sets

## ------------------------------------------------------------------------
A.vs.wt1021.DE_mapped <- string.db.sme1021$map( A.vs.wt1021.DE, "GeneSymbol", removeUnmappedRows = TRUE )
write.table(A.vs.wt1021.DE_mapped, file="A.vs.wt1021.KEGG.difftable")
head(A.vs.wt1021.DE_mapped)

AB.vs.wt1021B.DE_mapped <- string.db.sme1021$map( AB.vs.wt1021B.DE, "GeneSymbol", removeUnmappedRows = TRUE )
write.table(AB.vs.wt1021B.DE_mapped, file="AB.vs.wt1021B.KEGG.difftable")
head(AB.vs.wt1021B.DE_mapped)

AB.vs.wt1021.DE_mapped <- string.db.sme1021$map( AB.vs.wt1021.DE, "GeneSymbol", removeUnmappedRows = TRUE )
write.table(AB.vs.wt1021.DE_mapped, file="AB.vs.wt1021.KEGG.difftable")
head(AB.vs.wt1021.DE_mapped)

A.vs.AB.DE_mapped <- string.db.sme1021$map( A.vs.AB.DE, "GeneSymbol", removeUnmappedRows = TRUE )
write.table(A.vs.AB.DE_mapped, file="A.vs.AB.KEGG.difftable")
head(A.vs.AB.DE_mapped)

wt1021B.vs.wt1021.DE_mapped <- string.db.sme1021$map( wt1021B.vs.wt1021.DE, "GeneSymbol", removeUnmappedRows = TRUE )
write.table(wt1021B.vs.wt1021.DE_mapped, file="wt1021B.vs.wt1021.KEGG.difftable")
head(tw1021B.vs.wt1021.DE_mapped)

# GO analysis adjusting for gene length bias
# (assuming that y$genes$Length contains gene lengths)
library(EnrichmentBrowser)
go.abund <- goa(A.vs.wt1021.DE, geneid = "GeneID", trend = T)
go.abund

go.len <- goanna(A.vs.wt1021.DE, geneid = "GeneID", trend = "Length")

topGO(go.len, sort = "Qval")
topGO(go.len, sort = "Pval")
#Avswt1021.up_reg_genes<-topGO(go.abund, sort = "Qval")
#tAvswt1021.down_reg_genes<-topGO(go.abund, sort = "Qval")


## Default usage with a list of gene sets:

go.de <- goana(list(DE1 = EG.DE1, DE2 = EG.DE2, DE3 = EG.DE3))
topGO(go.de, sort = "DE1")
topGO(go.de, sort = "DE2")
topGO(go.abund, ontology = "BP")
topGO(go.de, ontology = "CC", sort = "DE3")
topGO(go.de, ontology = "MF", sort = "DE3")

## Standard KEGG analysis

AB.vs.wt1021.DE.kegg <- kegga(AB.vs.wt1021.DE$GeneSymbol, species.KEGG="sme") # equivalent to previous
AB.vs.wt1021.DE.kegg

barplot(AB.vs.wt1021.DE.kegg$DE,drop=T, showCategory=8)
## ------------------------------------------------------------------------
dotplot(AB.vs.wt1021.DE.kegg)
## ----fig.cap="enrichment map of enrichment result", fig.align="center", fig.height=16, fig.width=16, eval=FALSE----
enrichMap(mkk)
## ## categorySize can be scaled by 'pvalue' or 'geneNum'
cnetplot(AB.vs.wt1021.DE_mapped$STRING_id, categorySize="pvalue", foldChange=AB.vs.wt1021.DE_mapped$LogFoldChange)

## ----fig.height=12, fig.width=8------------------------------------------
plotGOgraph(ego)

## ----fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8----
gseaplot(AB.vs.wt1021.DE.kegg, geneSetID = "sme")

head(ggo)
#####################################################################
# enrichment
###################################################

AB.vs.wt1021B.DE.df<-as.data.frame(cbind(gene=AB.vs.wt1021B.DE_mapped$GeneSymbol,
                                         pvalue=AB.vs.wt1021B.DE_mapped$Pval,
                                         logFC=AB.vs.wt1021B.DE_mapped$LogFoldChange), stringsAsFactors=F)
dim(AB.vs.wt1021B.DE.df)
head(AB.vs.wt1021B.DE.df)

AB.vs.wt1021B.intersected<-string.db.sme1021$map(AB.vs.wt1021B.DE.df, "gene", removeUnmappedRows=T)
head(AB.vs.wt1021B.intersected)
class(AB.vs.wt1021B.intersected)
string.db.sme1021$plot_network(AB.vs.wt1021B.intersected$STRING_id[1:400],)

AB.vs.wt1021B.subnets<-string.db.sme1021$get_subnetwork(AB.vs.wt1021B.intersected)
AB.vs.wt1021B.subnets

AB.vs.wt1021B.mapped_sig<-as.data.frame(cbind(genes=c(AB.vs.wt1021B.intersected$gene[AB.vs.wt1021B.intersected$pvalue < 0.05]),
                                              pvalue=c(AB.vs.wt1021B.intersected$pvalue[AB.vs.wt1021B.intersected$pvalue < 0.05]),
                                              logFC=c(AB.vs.wt1021B.intersected$logFC[AB.vs.wt1021B.intersected$pvalue < 0.05]),
                                              STRING_id=c(AB.vs.wt1021B.intersected$STRING_id[AB.vs.wt1021B.intersected$pvalue < 0.05])),
                                        stringsAsFactors=F, row.names=F)
head(AB.vs.wt1021B.mapped_sig)

AB.vs.wt1021B.DE.pv.fc.STRING<-as.data.frame(cbind(gene=AB.vs.wt1021B.DE_mapped$GeneSymbol,
                                                   pvalue=AB.vs.wt1021B.DE_mapped$Pval,
                                                   logFC=AB.vs.wt1021B.DE_mapped$LogFoldChange,
                                                   STRING_id=AB.vs.wt1021B.DE_mapped$STRING_id), stringsAsFactors=F, row.names=F, col.names=T)
head(AB.vs.wt1021B.DE.pv.fc.STRING)

# post payload information to the STRING server
AB.vs.wt1021B_pval01 <- string.db.sme1021$post_payload(AB.vs.wt1021B.mapped_sig$STRING_id,
                                                       colors=AB.vs.wt1021B.mapped_sig["pvalue"]$color )

# display a STRING network png with the "halo"

string.db.sme1021$plot_network( AB.vs.wt1021B.DE.pv.fc.STRING$STRING_id[1:50],
                                payload_id=AB.vs.wt1021B_pval01,
                                required_score=AB.vs.wt1021B.DE.pv.fc.STRING$logFC[1:50])

# plot the enrichment for the best 100 genes
ab.wt1021.top100<-string.db.sme1021$plot_ppi_enrichment( AB.vs.wt1021B.intersected$STRING_id[1:500], quiet=TRUE )

#############################################################################
#####################################################################
# enrichment
###################################################

A.vs.wt1021.DE.df<-as.data.frame(cbind(gene=A.vs.wt1021.DE_mapped$GeneSymbol,
                                       pvalue=A.vs.wt1021.DE_mapped$Pval,
                                       logFC=A.vs.wt1021.DE_mapped$LogFoldChange), stringsAsFactors=F)
dim(A.vs.wt1021.DE.df)
head(A.vs.wt1021.DE.df)

A.vs.wt1021.intersected<-string.db.sme1021$map(A.vs.wt1021.DE.df, "gene", removeUnmappedRows=T)
head(A.vs.wt1021.intersected)
class(A.vs.wt1021.intersected)
string.db.sme1021$plot_network(A.vs.wt1021.intersected$STRING_id[1:400],)

A.vs.wt1021.subnets<-string.db.sme1021$get_subnetwork(A.vs.wt1021.intersected)
A.vs.wt1021.subnets

A.vs.wt1021.mapped_sig<-as.data.frame(cbind(genes=c(A.vs.wt1021.intersected$gene[A.vs.wt1021.intersected$pvalue < 0.05]),
                                            pvalue=c(A.vs.wt1021.intersected$pvalue[A.vs.wt1021.intersected$pvalue < 0.05]),
                                            logFC=c(A.vs.wt1021.intersected$logFC[A.vs.wt1021.intersected$pvalue < 0.05]),
                                            STRING_id=c(A.vs.wt1021.intersected$STRING_id[A.vs.wt1021.intersected$pvalue < 0.05])),
                                      stringsAsFactors=F, row.names=F)
head(A.vs.wt1021.mapped_sig)

A.vs.wt1021.DE.pv.fc.STRING<-as.data.frame(cbind(gene=A.vs.wt1021.DE_mapped$GeneSymbol,
                                                 pvalue=A.vs.wt1021.DE_mapped$Pval,
                                                 logFC=A.vs.wt1021.DE_mapped$LogFoldChange,
                                                 STRING_id=A.vs.wt1021.DE_mapped$STRING_id), stringsAsFactors=F, row.names=F, col.names=T)
head(A.vs.wt1021.DE.pv.fc.STRING)
# post payload information to the STRING server
A.vs.wt1021_pval01 <- string.db.sme1021$post_payload(A.vs.wt1021.mapped_sig$STRING_id,
                                                     colors=A.vs.wt1021.mapped_sig["pvalue"]$color )

# display a STRING network png with the "halo"

string.db.sme1021$plot_network( A.vs.wt1021.DE.pv.fc.STRING$STRING_id[1:50],
                                payload_id=A.vs.wt1021_pval01,
                                required_score=A.vs.wt1021.DE.pv.fc.STRING$logFC[1:50])

# plot the enrichment for the best 100 genes
ab.wt1021.top100<-string.db.sme1021$plot_ppi_enrichment( A.vs.wt1021.intersected$STRING_id[1:500], quiet=TRUE )

#####################################################################
# enrichment
###################################################

AB.vs.wt1021.DE.df<-as.data.frame(cbind(gene=AB.vs.wt1021.DE_mapped$GeneSymbol,
                                        pvalue=AB.vs.wt1021.DE_mapped$Pval,
                                        logFC=AB.vs.wt1021.DE_mapped$LogFoldChange), stringsAsFactors=F)
dim(AB.vs.wt1021.DE.df)
head(AB.vs.wt1021.DE.df)

AB.vs.wt1021.intersected<-string.db.sme1021$map(AB.vs.wt1021.DE.df, "gene", removeUnmappedRows=T)
head(AB.vs.wt1021.intersected)
class(AB.vs.wt1021.intersected)
string.db.sme1021$plot_network(AB.vs.wt1021.intersected$STRING_id[1:400],)

AB.vs.wt1021.subnets<-string.db.sme1021$get_subnetwork(AB.vs.wt1021.intersected)
AB.vs.wt1021.subnets

AB.vs.wt1021.mapped_sig<-as.data.frame(cbind(genes=c(AB.vs.wt1021.intersected$gene[AB.vs.wt1021.intersected$pvalue < 0.05]),
                                             pvalue=c(AB.vs.wt1021.intersected$pvalue[AB.vs.wt1021.intersected$pvalue < 0.05]),
                                             logFC=c(AB.vs.wt1021.intersected$logFC[AB.vs.wt1021.intersected$pvalue < 0.05]),
                                             STRING_id=c(AB.vs.wt1021.intersected$STRING_id[AB.vs.wt1021.intersected$pvalue < 0.05])),
                                       stringsAsFactors=F, row.names=F)
head(AB.vs.wt1021.mapped_sig)

AB.vs.wt1021.DE.pv.fc.STRING<-as.data.frame(cbind(gene=AB.vs.wt1021.DE_mapped$GeneSymbol,
                                                  pvalue=AB.vs.wt1021.DE_mapped$Pval,
                                                  logFC=AB.vs.wt1021.DE_mapped$LogFoldChange,
                                                  STRING_id=AB.vs.wt1021.DE_mapped$STRING_id), stringsAsFactors=F, row.names=F, col.names=T)
head(AB.vs.wt1021.DE.pv.fc.STRING)
# post payload information to the STRING server
AB.vs.wt1021_pval01 <- string.db.sme1021$post_payload(AB.vs.wt1021.mapped_sig$STRING_id,
                                                      colors=AB.vs.wt1021.mapped_sig["pvalue"]$color )

# display a STRING network png with the "halo"

string.db.sme1021$plot_network( AB.vs.wt1021.DE.pv.fc.STRING$STRING_id[1:50],
                                payload_id=AB.vs.wt1021_pval01,
                                required_score=AB.vs.wt1021.DE.pv.fc.STRING$logFC[1:50])

# plot the enrichment for the best 100 genes
ab.wt1021.top100<-string.db.sme1021$plot_ppi_enrichment( AB.vs.wt1021.intersected$STRING_id[1:500], quiet=TRUE )

#############################################################################
#####################################################################
# enrichment
###################################################

A.vs.AB.DE.df<-as.data.frame(cbind(gene=A.vs.AB.DE_mapped$GeneSymbol,
                                   pvalue=A.vs.AB.DE_mapped$Pval,
                                   logFC=A.vs.AB.DE_mapped$LogFoldChange), stringsAsFactors=F)
dim(A.vs.AB.DE.df)
head(A.vs.AB.DE.df)

A.vs.AB.intersected<-string.db.sme1021$map(A.vs.AB.DE.df, "gene", removeUnmappedRows=T)
head(A.vs.AB.intersected)
class(A.vs.AB.intersected)
string.db.sme1021$plot_network(A.vs.AB.intersected$STRING_id[1:400],)

A.vs.AB.subnets<-string.db.sme1021$get_subnetwork(A.vs.AB.intersected)
A.vs.AB.subnets

A.vs.AB.mapped_sig<-as.data.frame(cbind(genes=c(A.vs.AB.intersected$gene[A.vs.AB.intersected$pvalue < 0.05]),
                                        pvalue=c(A.vs.AB.intersected$pvalue[A.vs.AB.intersected$pvalue < 0.05]),
                                        logFC=c(A.vs.AB.intersected$logFC[A.vs.AB.intersected$pvalue < 0.05]),
                                        STRING_id=c(A.vs.AB.intersected$STRING_id[A.vs.AB.intersected$pvalue < 0.05])),
                                  stringsAsFactors=F, row.names=F)
head(A.vs.AB.mapped_sig)

A.vs.AB.DE.pv.fc.STRING<-as.data.frame(cbind(gene=A.vs.AB.DE_mapped$GeneSymbol,
                                             pvalue=A.vs.AB.DE_mapped$Pval,
                                             logFC=A.vs.AB.DE_mapped$LogFoldChange,
                                             STRING_id=A.vs.AB.DE_mapped$STRING_id), stringsAsFactors=F, row.names=F, col.names=T)
head(A.vs.AB.DE.pv.fc.STRING)

# post payload information to the STRING server
A.vs.AB_pval01 <- string.db.sme1021$post_payload(A.vs.AB.mapped_sig$STRING_id,
                                                 colors=A.vs.AB.mapped_sig["pvalue"]$color )

# display a STRING network png with the "halo"

string.db.sme1021$plot_network( A.vs.AB.DE.pv.fc.STRING$STRING_id[1:50],
                                payload_id=A.vs.AB_pval01,
                                required_score=A.vs.AB.DE.pv.fc.STRING$logFC[1:50])

# plot the enrichment for the best 100 genes
ab.wt1021.top100<-string.db.sme1021$plot_ppi_enrichment( A.vs.AB.intersected$STRING_id[1:500], quiet=TRUE )

#############################################################################
#####################################################################
# enrichment
###################################################

wt1021.vs.wt1021B.DE.df<-as.data.frame(cbind(gene=wt1021.vs.wt1021B.DE_mapped$GeneSymbol,
                                             pvalue=wt1021.vs.wt1021B.DE_mapped$Pval,
                                             logFC=wt1021.vs.wt1021B.DE_mapped$LogFoldChange), stringsAsFactors=F)
dim(wt1021.vs.wt1021B.DE.df)
head(wt1021.vs.wt1021B.DE.df)

wt1021.vs.wt1021B.intersected<-string.db.sme1021$map(wt1021.vs.wt1021B.DE.df, "gene", removeUnmappedRows=T)
head(wt1021.vs.wt1021B.intersected)
class(wt1021.vs.wt1021B.intersected)
string.db.sme1021$plot_network(wt1021.vs.wt1021B.intersected$STRING_id[1:400],)

wt1021.vs.wt1021B.subnets<-string.db.sme1021$get_subnetwork(wt1021.vs.wt1021B.intersected)
wt1021.vs.wt1021B.subnets

wt1021.vs.wt1021B.mapped_sig<-as.data.frame(cbind(genes=c(wt1021.vs.wt1021B.intersected$gene[wt1021.vs.wt1021B.intersected$pvalue < 0.05]),
                                                  pvalue=c(wt1021.vs.wt1021B.intersected$pvalue[wt1021.vs.wt1021B.intersected$pvalue < 0.05]),
                                                  logFC=c(wt1021.vs.wt1021B.intersected$logFC[wt1021.vs.wt1021B.intersected$pvalue < 0.05]),
                                                  STRING_id=c(wt1021.vs.wt1021B.intersected$STRING_id[wt1021.vs.wt1021B.intersected$pvalue < 0.05])),
                                            stringsAsFactors=F, row.names=F)
head(wt1021.vs.wt1021B.mapped_sig)

wt1021.vs.wt1021B.DE.pv.fc.STRING<-as.data.frame(cbind(gene=wt1021.vs.wt1021B.DE_mapped$GeneSymbol,
                                                       pvalue=wt1021.vs.wt1021B.DE_mapped$Pval,
                                                       logFC=wt1021.vs.wt1021B.DE_mapped$LogFoldChange,
                                                       STRING_id=wt1021.vs.wt1021B.DE_mapped$STRING_id), stringsAsFactors=F, row.names=F, col.names=T)
head(wt1021.vs.wt1021B.DE.pv.fc.STRING)
# post payload information to the STRING server
wt1021.vs.wt1021B_pval01 <- string.db.sme1021$post_payload(wt1021.vs.wt1021B.mapped_sig$STRING_id,
                                                           colors=wt1021.vs.wt1021B.mapped_sig["pvalue"]$color )

# display a STRING network png with the "halo"

string.db.sme1021$plot_network( wt1021.vs.wt1021B.DE.pv.fc.STRING$STRING_id[1:50],
                                payload_id=wt1021.vs.wt1021B_pval01,
                                required_score=wt1021.vs.wt1021B.DE.pv.fc.STRING$logFC[1:50])

# plot the enrichment for the best 100 genes
wt1021.wt1021b.top100<-string.db.sme1021$plot_ppi_enrichment( wt1021.vs.wt1021B.intersected$STRING_id[1:500], quiet=TRUE )


## ----eval = FALSE--------------------------------------------------------
wt1021.vs.wt1021B.kegg.rich <- enrichKEGG(gene = wt1021.vs.wt1021B.DE_mapped$GeneSymbol,organism='sme',pvalueCutoff = 0.05)
head(wt1021.vs.wt1021B.kegg.rich)
dim(wt1021.vs.wt1021B.kegg.rich)
wt1021.vs.wt1021B.mkegg.rich <- enrichMKEGG(gene = wt1021.vs.wt1021B.DE_mapped$GeneSymbol,organism='sme',pvalueCutoff = 0.05)
head(wt1021.vs.wt1021B.mkegg.rich)
dim(wt1021.vs.wt1021B.mkegg.rich)

barplot(wt1021.vs.wt1021B.mkegg.rich, drop=TRUE, showCategory=12)
barplot(wt1021.vs.wt1021B.kegg.rich,drop=T, showCategory=12)
dotplot(wt1021.vs.wt1021B.mkegg.rich)
dotplot(wt1021.vs.wt1021B.kegg)
cnetplot(wt1021.vs.wt1021B.mkegg.rich, categorySize="pvalue") # ,wt1021.vs.wt1021B.kegg
enrichMap(wt1021.vs.wt1021B.mkegg.rich)
cnetplot(wt1021.vs.wt1021B.kegg.rich, categorySize="pvalue") # ,wt1021.vs.wt1021B.kegg
enrichMap(wt1021.vs.wt1021B.kegg.rich)
cnetplot(wt1021.vs.wt1021B.kegg.rich,categorySize="pvalue", foldChange=,wt1021.vs.wt1021B.DE_mapped$LogFoldChange, )
cnetplot(wt1021.vs.wt1021B.kegg.rich, categorySize="pvalue") # ,wt1021.vs.wt1021B.kegg
enrichMap(wt1021.vs.wt1021B.kegg.rich)


## ----fig.height=12, fig.width=8------------------------------------------
AB.vs.wt1021.mkegg <- enrichMKEGG(gene = AB.vs.wt1021.DE_mapped$GeneSymbol,
                                  organism = 'sme')
AB.vs.wt1021.kegg <- enrichKEGG(gene = AB.vs.wt1021.DE_mapped$GeneSymbol,
                                organism = 'sme')

barplot(AB.vs.wt1021.mkegg, drop=TRUE, showCategory=12)
barplot(AB.vs.wt1021.kegg, showCategory=8)
dotplot(AB.vs.wt1021.mkegg)
dotplot(AB.vs.wt1021.kegg)
cnetplot(AB.vs.wt1021.mkegg, categorySize="pvalue") # ,wt1021.vs.wt1021B.kegg
enrichMap(AB.vs.wt1021.mkegg)
cnetplot(AB.vs.wt1021.kegg, categorySize="pvalue",foldChange=AB.vs.wt1021.DE_mapped$LogFoldChange)
enrichMap(AB.vs.wt1021.kegg)
# GO analysis adjusting for gene length bias
# (assuming that y$genes$Length contains gene lengths)
library(EnrichmentBrowser)
go.abund <- goa(A.vs.wt1021.DE, geneid = "GeneID", trend = T)
go.abund

go.len <- goanna(A.vs.wt1021.DE, geneid = "GeneID", trend = "Length")

topGO(go.len, sort = "Qval")
topGO(go.len, sort = "Pval")
#Avswt1021.up_reg_genes<-topGO(go.abund, sort = "Qval")
#tAvswt1021.down_reg_genes<-topGO(go.abund, sort = "Qval")


## Default usage with a list of gene sets:

go.de <- goana(list(DE1 = EG.DE1, DE2 = EG.DE2, DE3 = EG.DE3))
topGO(go.de, sort = "DE1")
topGO(go.de, sort = "DE2")
topGO(go.abund, ontology = "BP")
topGO(go.de, ontology = "CC", sort = "DE3")
topGO(go.de, ontology = "MF", sort = "DE3")

## Standard KEGG analysis

AB.vs.wt1021.DE.kegg <- kegga(AB.vs.wt1021.DE$GeneSymbol, species.KEGG="sme") # equivalent to previous
AB.vs.wt1021.DE.kegg

barplot(AB.vs.wt1021.DE.kegg$DE,drop=T, showCategory=8)
## ------------------------------------------------------------------------
dotplot(AB.vs.wt1021.DE.kegg)
## ----fig.cap="enrichment map of enrichment result", fig.align="center", fig.height=16, fig.width=16, eval=FALSE----
enrichMap(mkk)
## ## categorySize can be scaled by 'pvalue' or 'geneNum'
cnetplot(AB.vs.wt1021.DE_mapped$STRING_id, categorySize="pvalue", foldChange=AB.vs.wt1021.DE_mapped$LogFoldChange)

## ----fig.height=12, fig.width=8------------------------------------------
plotGOgraph(ego)

## ----fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8----
gseaplot(AB.vs.wt1021.DE.kegg, geneSetID = "sme")

head(ggo)
AB.vs.wt1021.DE_gse<-as.data.frame(AB.vs.wt1021.DE_mapped$Pval, row.names=c(AB.vs.wt1021.DE_mapped$GeneSymbol), stringsAsFactors=F)
names(AB.vs.wt1021.DE_sorted.gse)<-AB.vs.wt1021.DE_mapped$GeneSymbol
AB.vs.wt1021.DE_sorted.gse<-sort(AB.vs.wt1021.DE.kegg$P.DE, decreasing=T)
kk2 <- gseKEGG(geneList     = AB.vs.wt1021.DE_sorted.gse,
               organism     = 'sme',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)
##############################################################################
### DOSE/ClusterProfile - KEGG and Gene Ontology analysis
###############################################################################
wt1021.vs.wt1021B.mkegg.rich <- enrichMKEGG(gene = wt1021.vs.wt1021B.DE_mapped$GeneSymbol,organism = 'sme')
wt1021.vs.wt1021B.kegg.rich <- enrichKEGG(gene = wt1021.vs.wt1021B.DE_mapped$GeneSymbol,organism='sme',pvalueCutoff = 0.05)
head(wt1021.vs.wt1021B.kegg.rich)
dim(wt1021.vs.wt1021B.kegg.rich)
wt1021.vs.wt1021B.DE.kegg.rich <- enrichMKEGG(gene = wt1021.vs.wt1021B.DE_mapped$GeneSymbol,organism='sme')
head(wt1021.vs.wt1021B.DE.kegg.rich)
dim(wt1021.vs.wt1021B.DE.kegg.rich)


AB.vs.wt1021.kegg.rich <- enrichMKEGG(gene = AB.vs.wt1021.DE_mapped$GeneSymbol,organism='sme',pvalueCutoff = 0.05)
head(AB.vs.wt1021.kegg.rich)
dim(AB.vs.wt1021.kegg.rich)
AB.vs.wt1021B.kegg.rich <- enrichKEGG(gene = AB.vs.wt1021B.DE_mapped$GeneSymbol,organism='sme',pvalueCutoff = 0.05)
head(AB.vs.wt1021B.kegg.rich)
dim(AB.vs.wt1021B.kegg.rich)
AB.vs.wt1021B.kegg.rich
A.vs.wt1021.kegg.rich <- enrichKEGG(gene = A.vs.wt1021.DE_mapped$GeneSymbol,organism='sme',pvalueCutoff = 0.05)
head(A.vs.wt1021.kegg.rich)
dim(A.vs.wt1021.kegg.rich)
A.vs.wt1021.kegg.rich

gene.df <- bitr(names(AB.vs.wt1021.DE_mapped$GeneSymbol), fromType = "SYMBOL",toType = c("ENTREZID", "KEGG"),OrgDb=sme.1021.kegg)



# GO analysis adjusting for gene length bias
# (assuming that y$genes$Length contains gene lengths)
library(EnrichmentBrowser)
go.abund <- goanna(sig.o.exp.df, geneid = "gene", trend = T)
go.abund

go.len <- goanna(A.vs.wt1021.DE, geneid = "GeneID", trend = "Length")

topGO(go.len, sort = "Qval")
topGO(go.len, sort = "Pval")
#Avswt1021.up_reg_genes<-topGO(go.abund, sort = "Qval")
#tAvswt1021.down_reg_genes<-topGO(go.abund, sort = "Qval")


## Default usage with a list of gene sets:

go.de <- goana(list(DE1 = EG.DE1, DE2 = EG.DE2, DE3 = EG.DE3))
topGO(go.de, sort = "DE1")
topGO(go.de, sort = "DE2")
topGO(go.abund, ontology = "BP")
topGO(go.de, ontology = "CC", sort = "DE3")
topGO(go.de, ontology = "MF", sort = "DE3")

## Standard KEGG analysis

sig.o..DE.kegg <- kegga(sig.over.DE_mapped$gene, species.KEGG="hsa") # equivalent to previous
sig.o..DE.kegg

barplot(sig.o..DE.kegg$DE,drop=T, showCategory=8)
## ------------------------------------------------------------------------
dotplot(sig.o..DE.kegg)
## ----fig.cap="enrichment map of enrichment result", fig.align="center", fig.height=16, fig.width=16, eval=FALSE----
enrichMap(sig.o..DE.kegg)
## ## categorySize can be scaled by 'pvalue' or 'geneNum'
cnetplot(sig.over.DE_mapped$STRING_id, categorySize="pvalue", foldChange=sig.over.DE_mapped$logFC)

## ----fig.height=12, fig.width=8------------------------------------------
plotGOgraph(ego)

## ----fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8----
gseaplot(sig.o..DE.kegg, geneSetID = "hsa")

head(ggo)
AB.vs.wt1021.DE_gse<-as.data.frame(AB.vs.wt1021.DE_mapped$Pval, row.names=c(AB.vs.wt1021.DE_mapped$GeneSymbol), stringsAsFactors=F)
names(AB.vs.wt1021.DE_sorted.gse)<-AB.vs.wt1021.DE_mapped$GeneSymbol
AB.vs.wt1021.DE_sorted.gse<-sort(AB.vs.wt1021.DE.kegg$P.DE, decreasing=T)
kk2 <- gseKEGG(geneList     = AB.vs.wt1021.DE_sorted.gse,
               organism     = 'sme',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)
## ------------------------------------------------------------------------
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)


## ------------------------------------------------------------------------
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)

## ------------------------------------------------------------------------
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)

## ----eval = FALSE--------------------------------------------------------
## mkk <- enrichMKEGG(gene = gene,
##                    organism = 'hsa')

## ----eval=FALSE----------------------------------------------------------
## mkk2 <- gseMKEGG(geneList = geneList,
##                  species = 'hsa')

## ----eval=FALSE----------------------------------------------------------
## david <- enrichDAVID(gene = gene,
##                      idType = "ENTREZ_GENE_ID",
##                      listType = "Gene",
##                      annotation = "KEGG_PATHWAY",
##                      david.user = "clusterProfiler@hku.hk")

## ----fig.height=5, fig.width=9-------------------------------------------
barplot(ggo, drop=TRUE, showCategory=12)

## ----fig.height=5, fig.width=8-------------------------------------------
barplot(mkk, showCategory=8)
## ------------------------------------------------------------------------
dotplot(mkk)
## ----fig.cap="enrichment map of enrichment result", fig.align="center", fig.height=16, fig.width=16, eval=FALSE----
enrichMap(mkk)
## ## categorySize can be scaled by 'pvalue' or 'geneNum'
cnetplot(mkk, categorySize="pvalue", foldChange=geneList)

## ----fig.height=12, fig.width=8------------------------------------------
plotGOgraph(ego)

## ----fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8----
gseaplot(kk, geneSetID = "sme")

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

## ----fig.height=6, fig.width=10------------------------------------------
dotplot(formula_res)
dotplot(formula_res, x=~group) + ggplot2::facet_grid(~othergroup)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"
# [8] "EVIDENCE"     "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"        "IPI"          "MAP"
#[15] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"         "PROSITE"
#[22] "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIGENE"      "UNIPROT"
eg = bitr(sig.over.DE_mapped$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg)
uniprot_ids <- bitr(glist, fromType="SYMBOL", toType=c("UNIPROT"), OrgDb="org.Hs.eg.db")
head(uniprot_ids)
refseq_ids <- bitr(glist, fromType="SYMBOL", toType=c("REFSEQ"), OrgDb="org.Hs.eg.db")
head(refseq_ids)
go_ids <- bitr(glist, fromType="SYMBOL", toType=c("UCSCKG"), OrgDb="org.Hs.eg.db")
head(go_ids)
go_ids <- bitr(glist, fromType="SYMBOL", toType=c("GOALL"), OrgDb="org.Hs.eg.db")
head(go_ids)
#eg2np <- bitr_kegg(glist, fromType='ncbi-geneid', toType='kegg', organism='hsa')
#bitr_kegg("Z5100", fromType="kegg", toType='ncbi-proteinid', organism='ece')
#bitr_kegg("Z5100", fromType="kegg", toType='uniprot', organism='ece')
library(DOSE)
na.omit(genelist)
gene <- names(glist)
gene.df <- bitr(glist, fromType = "SYMBOL",toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
str(gene.df)
entrezgenes<-gene.df[,"ENTREZID"]
ggo <- groupGO(gene=entrezgenes, OrgDb=org.Hs.eg.db, ont="CC",
               level    = 3,readable = TRUE)
head(ggo)
kk <- enrichKEGG(gene = entrezgenes,organism='hsa',pvalueCutoff = 0.05)
head(kk)

gene.df <- bitr(AB.vs.wt1021.DE$GeneSymbol, fromType = "SYMBOL",toType = c("ENTREZID", "KEGG"),OrgDb=Org.Hs.egOMIM2EG@datacache)

S.me1021 <- enrichKEGG(gene = geneList,organism='sme',pvalueCutoff = 0.05)
head(S.me1021)

#############################################################################


sig.genes.npid <- bitr_kegg(row.names(sig.genes_exp.diff), fromType='kegg', toType='ncbi-proteinid', organism='hsa',drop=T)
head(sig.genes.npid)
dim(AB.vs.wt1021.npid)

AB.vs.wt1021.geneid <- bitr_kegg(AB.vs.wt1021.DE$GeneSymbol, fromType='kegg', toType='ncbi-geneid', organism='sme')
dim(AB.vs.wt1021.geneid)


library(package = affyLib, character.only = TRUE)
## the distribution of the adjusted p-values
hist(geneList, 100)

## how many differentially expressed genes are:
sum(topDiffGenes(geneList))

## build the topGOdata class
GOdata <- new("topGOdata",ontology = "BP",
              allGenes = geneList,geneSel = topDiffGenes,
              annot = annFUN.db,affylib = affyLib)

## display the GOdata object
GOdata


############################################################################
# adj.matrix.graph-all-sig-genes ####
############################################################################

o.group<-grep(pattern = over, x = colnames(g.count.matrix),ignore.case = T)
u.group<-grep(pattern = under, x = colnames(g.count.matrix),ignore.case = T)

over.group<-grep(pattern = over, x = colnames(g.rep.matrix),ignore.case = T)
under.group<-grep(pattern = under, x = colnames(g.rep.matrix),ignore.case = T)


library(igraph)

#Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
g <- graph.adjacency(as.matrix(as.dist(cor(t(MyData), method="pearson"))), mode="undirected", weighted=TRUE, diag=FALSE)

You can also use Euclidean distances instead of correlation, with:

  g <- graph.adjacency(as.matrix(dist(MyData, method="euclidean")), mode="undirected", weighted=TRUE, diag=FALSE)

#########################################################################################################
# GOplots-genelists ####
#########################################################################################################

geneList<-mapIds(x = org.Hs.eg.db,
                 keys =  rownames(g.rep.matrix),
                 column = "ENTREZID",
                 keytype = "SYMBOL",
                 multiVals="first")

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

head(sig.genes_exp.diff)
sig_gene_fold.df<-sig_genes_exp.diff[order(sig_genes_exp.diff$log2_fold_change,decreasing = T),]
sig.gene.logfold.df<-as.data.frame(sig_gene_fold.df)
rownames(sig.gene.logfold.df) = sigENTREZsiggenes
symbols<-names(ENTREZsiggenes)

#over.hilog<-c(sig_genes_exp.diff[,"log2_fold_change"])
#names(over.hilog)<-(ENTREZQvalOVERhi)
#over.hiqval<-sig_genes_exp.diff[,"q_value"]
#names(over.hiqval)<-ENTREZQvalOVERhi
#over.hipval<-sig_genes_exp.diff[,"p_value"]
#names(over.hipval)<-ENTREZQvalOVERhi
# rownames(sig_gene_fold.df)<-
#names(under.hilog)<-ENTREZsiggenes
#under.hiqval<-sig_genes_exp.diff[,"q_value"]
#names(under.hiqval)<-ENTREZQvalUNDERhi
#under.hipval<-(sig_genes_exp.diff[,"p_value"])
#names(under.hipval) <- ENTREZQvalUNDERhi

xx<-annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
topDiffGenes <- function(allScore) {
  return(allScore < 0.01)}

#########################################################################################################
# GOplots-go-objects ####
#########################################################################################################

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

#########################################################################################################
# GOplots-under-plots ####
#########################################################################################################

pdf("GOplots_grch38.CTRL_highly_expressed.pdf")

showSigOfNodes(HIunder.BP.GOdata, score(HIunder.BPFisher), firstSigNodes = 5, useInfo = "all")
showSigOfNodes(HIunder.BP.GOdata, score(HIunder.BPtKS), firstSigNodes = 5, useInfo = "all")
showSigOfNodes(HIunder.BP.GOdata, score(HIunder.BPtKS.elim), firstSigNodes = 5, useInfo = "all")
#printGraph(HIunder.BP.GOdata, HIunder.BPFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
#printGraph(HIunder.BP.GOdata, HIunder.BPtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
#printGraph(HIunder.BP.GOdata, HIunder.BPtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(HIunder.MF.GOdata, score(HIunder.MFFisher), firstSigNodes = 5, useInfo = "all")
showSigOfNodes(HIunder.MF.GOdata, score(HIunder.MFtKS), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(HIunder.MF.GOdata, score(HIunder.MFFtKS.elim), firstSigNodes = 5, useInfo = "all")
#printGraph(HIunder.MF.GOdata, HIunder.MFFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
#printGraph(HIunder.MF.GOdata, HIunder.MFFtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
#printGraph(HIunder.MF.GOdata, HIunder.MFtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(HIunder.CC.GOdata, score(HIunder.CCFisher), firstSigNodes = 5, useInfo = "all")
showSigOfNodes(HIunder.CC.GOdata, score(HIunder.CCtKS), firstSigNodes = 5, useInfo = "all")
showSigOfNodes(HIunder.CC.GOdata, score(HIunder.CCtKS.elim), firstSigNodes = 5, useInfo = "all")
#printGraph(HIunder.CC.GOdata, HIunder.CCFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
#printGraph(HIunder.CC.GOdata, HIunder.CCtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
#printGraph(HIunder.CC.GOdata, HIunder.CCtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

sigGOcc <- groupGO(gene=ENTREZQval.hi_exprUNDER, OrgDb=org.Hs.eg.db, ont="CC",
                   level    = 3,readable = TRUE)
sigGOmf <- groupGO(gene=ENTREZQval.hi_exprUNDER, OrgDb=org.Hs.eg.db, ont="MF",
                   level    = 3,readable = TRUE)
sigGObp <- groupGO(gene=ENTREZQval.hi_exprUNDER, OrgDb=org.Hs.eg.db, ont="BP",
                   level    = 3,readable = TRUE)
barplot(sigGOcc, drop=TRUE, showCategory=12,args.legend="Sig Diff Expr Genes Grouping by GO-CC")
barplot(sigGOmf, drop=TRUE, showCategory=12,args.legend="Sig Diff Expr Genes Grouping by GO-MF")
barplot(sigGObp, drop=TRUE, showCategory=12,args.legend="Sig Diff Expr Genes Grouping by GO-BP")

dev.off()
#########################################################################################################
# GOplots-over-plots ####
#########################################################################################################

HIover.BP.GOdata <- new("topGOdata",ontology = "BP",allGenes = Qval.hi_exprOVER, nodeSize = 5,
                        annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                        ID = "symbol")
HIover.MF.GOdata <- new("topGOdata",ontology = "MF",allGenes = Qval.hi_exprOVER, nodeSize = 5,
                        annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                        ID = "symbol")

HIover.CC.GOdata <- new("topGOdata",ontology = "CC",allGenes = Qval.hi_exprOVER, nodeSize = 5,
                        annot = annFUN.org,mapping = "org.Hs.eg.db",geneSel = topDiffGenes,
                        ID = "symbol")


HIover.BPtKS <- runTest(HIover.BP.GOdata, algorithm = "classic", statistic = "ks")
HIover.BPtKS
HIover.BPFisher <- runTest(HIover.BP.GOdata, algorithm = "classic", statistic = "fisher")
HIover.BPFisher
HIover.BPtKS.elim <- runTest(HIover.BP.GOdata, algorithm = "elim", statistic = "ks")
HIover.BPtKS.elim

HIover.MFtKS <- runTest(HIover.MF.GOdata, algorithm = "classic", statistic = "ks")
HIover.MFtKS
HIover.MFFisher <- runTest(HIover.MF.GOdata, algorithm = "classic", statistic = "fisher")
HIover.MFFisher
HIoverMFtKS.elim <- runTest(HIover.MF.GOdata, algorithm = "elim", statistic = "ks")
HIoverMFtKS.elim

HIover.CCtKS <- runTest(HIover.CC.GOdata, algorithm = "classic", statistic = "ks")
HIover.CCtKS
HIover.CCFisher <- runTest(HIover.CC.GOdata, algorithm = "classic", statistic = "fisher")
HIover.CCFisher
HIover.CCtKS.elim <- runTest(HIover.CC.GOdata, algorithm = "elim", statistic = "ks")
HIover.CCtKS.elim

pdf("GOplots_grch38.LUTS_highly_expressed.pdf")

showSigOfNodes(HIover.BP.GOdata, score(HIover.BPFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(HIover.BP.GOdata, score(HIover.BPtKS), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(HIover.BP.GOdata, score(HIover.BPtKS.elim), firstSigNodes = 5, useInfo = "all" )
#printGraph(HIover.BP.GOdata, HIover.BPFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
#printGraph(HIover.BP.GOdata, HIover.BPtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
#printGraph(HIover.BP.GOdata, HIover.BPtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(HIover.MF.GOdata, score(HIover.MFFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(HIover.MF.GOdata, score(HIover.MFtKS), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(HIover.MF.GOdata, score(HIoverMFtKS.elim), firstSigNodes = 5, useInfo = "all" )
#printGraph(HIover.MF.GOdata, HIover.MFFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
#printGraph(HIover.MF.GOdata, HIoverMFtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
#printGraph(HIover.MF.GOdata, HIover.MFtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

showSigOfNodes(HIover.CC.GOdata, score(HIover.CCFisher), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(HIover.CC.GOdata, score(HIover.CCtKS), firstSigNodes = 5, useInfo = "all" )
showSigOfNodes(HIover.CC.GOdata, score(HIover.CCtKS.elim), firstSigNodes = 5, useInfo = "all" )
#printGraph(HIover.CC.GOdata, HIover.CCFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
#printGraph(HIover.CC.GOdata, HIover.CCtKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
#printGraph(HIover.CC.GOdata, HIover.CCtKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

sigGOcc <- groupGO(gene=ENTREZQval.hi_exprOVER, OrgDb=org.Hs.eg.db, ont="CC",
                   level    = 3,readable = TRUE)
sigGOmf <- groupGO(gene=ENTREZQval.hi_exprOVER, OrgDb=org.Hs.eg.db, ont="MF",
                   level    = 3,readable = TRUE)
sigGObp <- groupGO(gene=ENTREZQval.hi_exprOVER, OrgDb=org.Hs.eg.db, ont="BP",
                   level    = 3,readable = TRUE)

barplot(sigGOcc, drop=TRUE, showCategory=12,args.legend="Sig Diff Expr Genes Grouping by GO-CC")
barplot(sigGOmf, drop=TRUE, showCategory=12,args.legend="Sig Diff Expr Genes Grouping by GO-MF")
barplot(sigGObp, drop=TRUE, showCategory=12,args.legend="Sig Diff Expr Genes Grouping by GO-BP")

dev.off()
