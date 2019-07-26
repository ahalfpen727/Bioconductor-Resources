### R code from vignette source 'NPMVignette.Rnw'
### Encoding: UTF-8
library(genomeIntervals)
library(qpgraph)
library(GGBase)
library(GSRI)
library(ArrayExpress)
library(keggorthology)
library(MLInterfaces)
library(biclust)
library(eisa)
library(ALL)
library(biclust)
library(ggnetwork)
library(IRanges)
library(GO.db)
library(GenomeGraphs)
library(ggbio)
library(EGSEA)
###################################################
### code chunk number 1: no.nonsense
###################################################
rm(list=ls())


library(GGtools)
library(GeneR)
library(Path2enet)
library(gene2pathway)
library(ggstatsplot)
library(enrichplot)
library("iterativeBMA")
###################################################
### code chunk number 2: Load_package
###################################################
library(NetPathMiner)
source("https://bioconductor.org/biocLite.R")
biocLite(c("Rsubread","rBiopaxParser", "psichomics"))
biocLite(c("CellNOptR", "NetPathMiner", "airway", "paxtoolsr", "Path2enet", "qpgraph", "gaggle", "genomeIntervals", "GGBase", "ArrayExpress",
           "MLInterfaces", "keggorthology", "geneLenDataBase", "biclust", "eisa", "ALL","GenomeGraphs", "enrichplot", "ggstatsplot", "GSRI" ))

getwd()
install.packages( "/home/drew/R/x86_64-pc-linux-gnu-library/3.4/gogadget_2.1.tar.gz", repos = NULL, type="source")
library(gogadget)
library(ggnetwork)
library(topGO)
library(GenomicFeatures)
library(GenomicRanges)
library(KEGGgraph)
###################################################
### code chunk number 3: NPMVignette.Rnw:235-237 (eval = FALSE)
###################################################
## graph <- KGML2igraph(filename = file)
## graph <- SBML2igraph(filename = file)
library("systemPipeR")
library(Rsubread)
###################################################
### code chunk number 4: NPMVignette.Rnw:243-246 (eval = FALSE)
###################################################
## require(rBiopaxParser)
## biopax = readBiopax(file)
## graph <- BioPAX2igraph(biopax = biopax)


###################################################
### code chunk number 5: NPMVignette.Rnw:252-253 (eval = FALSE)
###################################################
## graph <- KGML2igraph(filename = c(file1, file2))


###################################################
### code chunk number 6: NPMVignette.Rnw:257-258 (eval = FALSE)
###################################################
## graph <- KGML2igraph(filename = ".")


###################################################
### code chunk number 7: NPMVignette.Rnw:263-268 (eval = FALSE)
###################################################
## # Extract all MIRIAM identifiers from an SBML file.
## graph <- SBML2igraph(filename = file, miriam = "all")
##
## # Extract all MIRIAM identifiers from an SBML file.
## graph <- BioPAX2igraph(biopax = biopax, miriam = "go")


###################################################
### code chunk number 8: NPMVignette.Rnw:275-276
###################################################
file <- file.path(find.package("NetPathMiner"), "extdata", "hsa00860.xml")


###################################################
### code chunk number 9: NPMVignette.Rnw:278-282 (eval = FALSE)
###################################################
## graph <- KGML2igraph(filename = file, parse.as = "signaling")
##
## graph <- KGML2igraph(filename = file, parse.as = "signaling",
## 	expand.complexes = TRUE)


###################################################
### code chunk number 10: NPMVignette.Rnw:288-291
###################################################
data("ex_sbml")
graph <- ex_sbml
graph


###################################################
### code chunk number 11: NPMVignette.Rnw:300-301
###################################################
head( V(graph) )


###################################################
### code chunk number 12: NPMVignette.Rnw:304-305
###################################################
head( E(graph) )


###################################################
### code chunk number 13: NPMVignette.Rnw:308-309
###################################################
head( V(graph)[ reactions ] )


###################################################
### code chunk number 14: NPMVignette.Rnw:314-315
###################################################
V(graph)[ "reaction_71850" ]$attr


###################################################
### code chunk number 15: NPMVignette.Rnw:322-323
###################################################
getAttrNames(graph)


###################################################
### code chunk number 16: NPMVignette.Rnw:330-331
###################################################
getAttrStatus(graph, pattern = "^miriam.")


###################################################
### code chunk number 17: NPMVignette.Rnw:337-345 (eval = FALSE)
###################################################
## require("RCurl")
## # Fetch uniprot annotation
## graph <- fetchAttribute(graph, organism = "Homo sapiens",
## target.attr = "miriam.ncbigene" , source.attr = "miriam.uniprot")
##
## # Fetch ChEBI annotation.
## graph <- fetchAttribute(graph, target.attr = "miriam.chebi",
## source.attr = "miriam.kegg.compound")


###################################################
### code chunk number 18: NPMVignette.Rnw:355-357
###################################################
rgraph <- makeReactionNetwork(graph, simplify=FALSE)
rgraph


###################################################
### code chunk number 19: NPMVignette.Rnw:362-364 (eval = FALSE)
###################################################
## rgraph <- simplifyReactionNetwork(rgraph)
## rgraph <- makeReactionNetwork(graph, simplify=TRUE)


###################################################
### code chunk number 20: NPMVignette.Rnw:369-375
###################################################
# Expand complexes of gene network.
ggraph <- expandComplexes(rgraph, v.attr = "miriam.uniprot",
		keep.parent.attr= c("^pathway", "^compartment"))

# Convert reaction network to gene network.
ggraph <- makeGeneNetwork(rgraph)


###################################################
### code chunk number 21: NPMVignette.Rnw:387-389
###################################################
data(ex_microarray)



###################################################
### code chunk number 22: NPMVignette.Rnw:390-395 (eval = FALSE)
###################################################
## # Assign weights to edges.
## if(require("RCurl") && url.exists( NPMdefaults("bridge.web") ))
## 	rgraph <- fetchAttribute(rgraph, organism = "Homo sapiens",
## 						target.attr = "miriam.affy.probeset",
## 						source.attr = "miriam.uniprot")


###################################################
### code chunk number 23: NPMVignette.Rnw:405-409 (eval = FALSE)
###################################################
## library(ALL)
## data(ALL)
## rgraph <- assignEdgeWeights(microarray = exprs(ALL), graph = rgraph,
## weight.method = "cor", use.attr="miriam.affy.probeset", y=ALL$mol.bio, bootstrap = FALSE)


###################################################
### code chunk number 24: NPMVignette.Rnw:413-416
###################################################
data(ex_microarray)
rgraph <- assignEdgeWeights(microarray = ex_microarray, graph = rgraph,
weight.method = "cor", use.attr="miriam.uniprot", y=colnames(ex_microarray), bootstrap = FALSE)


###################################################
### code chunk number 25: NPMVignette.Rnw:420-422
###################################################
rgraph$y.labels
head( E(rgraph)$edge.weights )


###################################################
### code chunk number 26: NPMVignette.Rnw:431-433
###################################################
ranked.p <- pathRanker(rgraph, method = "prob.shortest.path",
	K = 25, minPathSize = 6)


###################################################
### code chunk number 27: NPMVignette.Rnw:438-443 (eval = FALSE)
###################################################
## pathsample <- samplePaths(rgraph, max.path.length = vcount(rgraph),
## num.samples = 1000, num.warmup = 10)
##
## ranked.p <- pathRanker(rgraph, method = "pvalue",
## sampledpaths = pathsample ,alpha=0.1)


###################################################
### code chunk number 28: NPMVignette.Rnw:448-450
###################################################
# Get paths as edge IDs.
eids <- getPathsAsEIDs(paths = ranked.p, graph = rgraph)


###################################################
### code chunk number 29: NPMVignette.Rnw:455-457
###################################################
# Convert paths to other networks.
eids <- getPathsAsEIDs(paths = ranked.p, graph = ggraph)


###################################################
### code chunk number 30: NPMVignette.Rnw:464-467
###################################################
# Clustering.
ybinpaths <- pathsToBinary(ranked.p)
p.cluster <- pathCluster(ybinpaths, M = 2)


###################################################
### code chunk number 31: NPMVignette.Rnw:469-470
###################################################
plotClusters(ybinpaths, p.cluster)


###################################################
### code chunk number 32: NPMVignette.Rnw:475-476
###################################################
p.class <- pathClassifier(ybinpaths, target.class = "BCR/ABL", M = 2)


###################################################
### code chunk number 33: NPMVignette.Rnw:478-479 (eval = FALSE)
###################################################
## plotClassifierROC(p.class)


###################################################
### code chunk number 34: NPMVignette.Rnw:487-488
###################################################
plotClusters(ybinpaths, p.class)


###################################################
### code chunk number 35: NPMVignette.Rnw:496-497
###################################################
plotNetwork(rgraph, vertex.color="compartment.name")


###################################################
### code chunk number 36: NPMVignette.Rnw:502-506 (eval = FALSE)
###################################################
## plotPaths(ranked.p, rgraph)
##
## # With clusters
## plotPaths(ranked.p, graph, path.clusters=p.class)


###################################################
### code chunk number 37: NPMVignette.Rnw:511-513
###################################################
plotAllNetworks(ranked.p, metabolic.net = graph, reaction.net = rgraph,
		path.clusters=p.class, vertex.label = "", vertex.size = 4)


###################################################
### code chunk number 38: NPMVignette.Rnw:518-522 (eval = FALSE)
###################################################
## layout.c <- clusterVertexByAttr(rgraph, "pathway", cluster.strength = 3)
## v.color <- colorVertexByAttr(rgraph, "pathway")
## plotPaths(ranked.p , rgraph, clusters=p.class,
## 	layout = layout.c, vertex.color = v.color)


###################################################
### code chunk number 39: NPMVignette.Rnw:527-529 (eval = FALSE)
###################################################
## plotCytoscapeGML(graph, file="example.gml", layout = layout.c,
## 				vertex.size = 5, vertex.color = v.color)


###################################################
### code chunk number 40: NPMVignette.Rnw:536-537
###################################################
getGeneSets(graph, use.attr="compartment", gene.attr="miriam.uniprot")


###################################################
### code chunk number 41: NPMVignette.Rnw:542-543
###################################################
getGeneSetNetworks(graph, use.attr="compartment")


###################################################
### code chunk number 42: NPMVignette.Rnw:549-550 (eval = FALSE)
###################################################
## graphNEL <- toGraphNEL(graph, export.attr="^miriam.")

###################################
## Basic usage of GO information ##
###################################
library(GOstats); library(GO.db); library(ath1121501.db); library(annotate) # Loads the required libraries.
goann <- as.list(GOTERM) # Retrieves full set of GO annotations.
zz <- eapply(GOTERM, function(x) x@Ontology); table(unlist(zz)) # Calculates the number of annotations for each ontology category.
?GOTERM # To find out, how to access the different GO components.
GOTERM$"GO:0003700"; GOMFPARENTS$"GO:0003700"; GOMFCHILDREN$"GO:0003700"
# Shows how to print out the GO annotations for one entry and how to retrieve its direct parents and children.
GOMFANCESTOR$"GO:0003700"; GOMFOFFSPRING$"GO:0003700" # Prints out complete lineages of parents and children for a GO ID.
goterms <- unlist(eapply(GOTERM, function(x) x@Term)); goterms[grep("molecular_function", goterms)]
# Retrieves all GO terms and prints out only those matching a search string given in the grep function. The same can
# be done for the definition field with 'x@Definition'. A set of GO IDs can be provided as well: goterms[GOMFANCESTOR$"GO:0005507"]
go_df <- data.frame(GOID=unlist(eapply(GOTERM, function(x) x@GOID)), Term=unlist(eapply(GOTERM, function(x) x@Term)), Ont=unlist(eapply(GOTERM, function(x) x@Ontology)))
# Generates data frame of the commonly used GO components: GOID, GO Term and Ontology Type.
affyGO <- eapply(ath1121501GO, getOntology, "MF"); table(sapply(affyGO, length))
# Retrieves MF GO terms for all probe IDs of a chosen Affy chip and calculates how many probes have multiple GO terms
# associated. Use "BP" and "CC" arguments to retrieve BP/CC GO terms.
affyGOdf <- data.frame(unlist(affyGO)); affyGOdf <- data.frame(AffyID=row.names(affyGOdf), GOID=affyGOdf[,1]); affyGOdf <- merge(affyGOdf, go_df, by.x="GOID", by.y="GOID", all.x=T)
# Converts above MF list object into a data frame. The AffyID occurence counts are appended to AffyIDs. The last step
# merges the two data frames: 'affyGOdf' and 'go_df'.
unique(lookUp("GO:0004713", "ath1121501", "GO2ALLPROBES")) # Retrieves all Affy IDs that are associated with a GO node.
z <- affyGO[c("254759_at", "260744_at")]; as.list(GOTERM)[z[[1]]]
# Retrieves GO IDs for set of Affy IDs and then the corresponding GO term for first Affy ID.
a <- data.frame(unlist(z)); a <- data.frame(ID=row.names(a), a); b <- data.frame(goterms[as.vector(unlist(z))]); b <- data.frame(ID=row.names(b), b); merge(b, a, by.x = "ID", by.y="unlist.z.")
# Merges Affy ID, GO ID and GO annotation information.
affyEv <- eapply(ath1121501GO, getEvidence); table(unlist(affyEv, use.names = FALSE))
# Provides evidence code information for each gene and summarizes the result.
test1 <- eapply(ath1121501GO, dropECode, c("IEA", "NR")); table(unlist(sapply(test1, getEvidence), use.names = FALSE))
# This example shows how one can remove certain evidence codes (e.g. IEA, IEP) from the analysis.

##############################################
## GO term enrichment analysis with GOstats ##
##############################################
## Example of how to test a sample set of probe set keys for over-representation of GO terms using a hypergeometric distribution
## test with the function hyperGTest(). For more information, read the GOstatsHyperG manual.
library(ath1121501.db); library(ath1121501cdf)
affySample <- c("266592_at", "266703_at", "266199_at", "246949_at", "267370_at", "267115_s_at", "266489_at", "259845_at", "266295_at", "262632_at")
geneSample <- as.vector(unlist(mget(affySample, ath1121501ACCNUM, ifnotfound=NA)))
affyUniverse <- ls(ath1121501cdf)
geneUniverse <- as.vector(unlist(mget(affyUniverse, ath1121501ACCNUM, ifnotfound=NA)))
params <- new("GOHyperGParams", geneIds = geneSample, universeGeneIds = geneUniverse, annotation="ath1121501", ontology = "MF", pvalueCutoff = 0.5, conditional = FALSE, testDirection = "over")
hgOver <- hyperGTest(params)
summary(hgOver)
htmlReport(hgOver, file = "MyhyperGresult.html")
######################################################################################################
### GOHyperGAll: Global Hypergeometric Test Using Custom Gene-to-GO Mappings Plus GO Slim Analysis ###
######################################################################################################
## Author: Thomas Girke
## Last update: Feb 7, 2008
## Utility: To test a sample population of genes for over-representation of GO terms, the
## function 'GOHyperGAll' computes for all GO nodes a hypergeometric distribution test and
## returns the corresponding raw and Bonferroni corrected p-values. A subsequent filter function
## performs a GO Slim analysis using default or custom GO Slim categories.
## The associated publication is available in Plant Physiol (2008) 147, 41-57.
## Note: GOHyperGAll provides similar utilities as the GOHyperG function in the GOstats package
## from BioConductor. The main difference is that GOHyperGAll simplifies the usage of custom
## chip-to-gene and gene-to-GO mappings.
##
## How it works:
## (A) Generate the required data objects (slow, but needs to be done only once)
## (B) Define GOhyperG_All function
## (C) Subsetting and plotting of results by assigned nodes or goSlim categories
## To demo the script and import all required functions, run the following source() command:
##         source("http://bioinfo.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/GOHyperGAll.txt")
## HTML Instructions:
##       http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/R_BioCondManual.html#GOHyperGAll

#########################################
## (A) Generate the required data objects
#########################################
## (A.1) Generate sample data frames with assigned gene-to-GO mappings,
## one for MF, one for BP and one for CC mappings
## custom mappings can be used here, but need to have the same format as GO_XX_DFs in the following examples
## (A.1.1) Obtain mappings from geneontology.org
readGOorg <- function(myfile = "gene_association.tair", colno = c(5,11,9), org) {
   go_org <- read.delim(myfile, na.strings = "", header=F, comment.char = "!", sep="\t")
   go_org <- go_org[ , colno]
   names(go_org) <- c("GOID", "GeneID", "GOCAT")
   if(org == "Arabidopsis") {
      go_org[,"GeneID"] <- gsub(".*(AT.G\\d\\d\\d\\d\\d).*", "\\1", as.character(go_org[,2]), perl=T)
      go_org <- go_org[grep("^AT.G\\d\\d\\d\\d\\d", as.character(go_org$GeneID), perl=T),]
      go_org <- go_org[!duplicated(paste(go_org[,"GOID"], gsub("\\.\\d{1,}", "", as.character(go_org[,"GeneID"]), perl=T), sep="_")),]
   }
   go_org <- na.omit(go_org)
   GO_MF_DF <<- go_org[go_org[,3]=="F",]
   write.table(GO_MF_DF, file="GO_MF_DF", quote=T, sep="\t")
   cat("\n", "Object 'GO_MF_DF' created containing assigned gene-to-MFGO mappings. To use custom mappings, generate data frame with the same structure in col 1-2.", "\n")
   GO_BP_DF <<- go_org[go_org[,3]=="P",]
   write.table(GO_BP_DF, file="GO_BP_DF", quote=T, sep="\t")
   cat("\n", "Object 'GO_BP_DF' created containing assigned gene-to-MFGO mappings. To use custom mappings, generate data frame with the same structure in col 1-2.", "\n")
   GO_CC_DF <<- go_org[go_org[,3]=="C",]
   write.table(GO_CC_DF, file="GO_CC_DF", quote=T, sep="\t")
   cat("\n", "Object 'GO_CC_DF' created containing assigned gene-to-MFGO mappings. To use custom mappings, generate data frame with the same structure in col 1-2.", "\n")

   ## Generates "go_df" data frame containing the commonly used components for all GO nodes: GOID, GO Term and Ontology Type.
   require(GOstats); require(GO.db)
   go_df <- data.frame(GOID=unlist(eapply(GOTERM, function(x) x@GOID)), Term=unlist(eapply(GOTERM, function(x) x@Term)), Ont=unlist(eapply(GOTERM, function(x) x@Ontology)))
   go_df <- na.omit(go_df)
   go_df <<- go_df
   write.table(go_df, file="go_df", quote=T, sep="\t")
   cat("\n", "Object 'go_df' created containing for all GO nodes the commonly used components: GOID, GO Term and Ontology Type", "\n")
}

## (A.1.1b) Convert GO objects from GOHyperGAll into GSEA format
## This step is not required for for GOHyperGAll approach. It simply converts XX_node_affy_list or GO_XX_DF files into *.gmt
## formatted files that can be imported into GSEA from the Broad Institute.
GOhyper2GSEA <- function(myfile=c("MF_node_affy_list", "BP_node_affy_list", "CC_node_affy_list"), type="all") {
   if(type=="all"){
      for(i in 1:length(myfile)) {
         mynames <- gsub("(..)_.*", "GO_\\1_ALL", myfile)
         load(file=myfile[i])
         GO_List <- eval(parse(text=myfile[i]))
         myindex <- sapply(GO_List, length)
         GO_List <- GO_List[myindex!=0]
         GO_List <- lapply(GO_List, paste, collapse="\t")
         exportDF <- data.frame(GS_ID=names(GO_List), Desc=rep("NA", length(GO_List)), GeneID=as.vector(unlist(GO_List)))
         write.table(exportDF, file=paste(mynames[i], ".gmt", sep=""), row.names=F, col.names=F, quote=F, sep="\t")
         cat("Saved file:", mynames[i], "\n")
      }
   }
   if(type=="terminal") {
      for(i in 1:length(myfile)) {
         mynames <- gsub("_DF", "_TERM", myfile)
         GO_DF <- read.delim(file=myfile[i])
         GO_List <- tapply(as.vector(GO_DF[,2]), as.factor(as.vector(GO_DF[ ,1])), as.vector)
         GO_List <- lapply(GO_List, paste, collapse="\t")
         exportDF <- data.frame(GS_ID=names(GO_List), Desc=rep("NA", length(GO_List)), GeneID=as.vector(unlist(GO_List)))
         write.table(exportDF, file=paste(mynames[i], ".gmt", sep=""), row.names=F, col.names=F, quote=F, sep="\t")
         cat("Saved file:", mynames[i], "\n")
      }
   }
}

## (A.1.2) Obtain mappings from BioC
sampleDFgene2GO <- function(lib="ath1121501.db") {
   require(GOstats); require(GO.db); require(annotate); require(lib, character.only=T)
   mylibbase <- gsub(".db", "", lib)
   affyGOMF <- eapply(get(paste(mylibbase, "GO", sep="")), getOntology, "MF") # generates list with GeneID components containing MFGOs
   GO_MF_DF <<- data.frame(GOID=unlist(affyGOMF), GeneID=rep(names(affyGOMF), as.vector(sapply(affyGOMF, length))), Count=rep(as.vector(sapply(affyGOMF, length)), as.vector(sapply(affyGOMF, length))))
   write.table(GO_MF_DF, file="GO_MF_DF", quote=F, sep="\t")
   cat("\n", "Object 'GO_MF_DF' created containing assigned gene-to-MFGO mappings. To use custom mappings, generate data frame with the same structure in col 1-2.", "\n")
   affyGOBP <- eapply(get(paste(mylibbase, "GO", sep="")), getOntology, "BP") # generates list with GeneID components containing BPGOs
   GO_BP_DF <<- data.frame(GOID=unlist(affyGOBP), GeneID=rep(names(affyGOBP), as.vector(sapply(affyGOBP, length))), Count=rep(as.vector(sapply(affyGOBP, length)), as.vector(sapply(affyGOBP, length))))
   write.table(GO_BP_DF, file="GO_BP_DF", quote=F, sep="\t")
   cat("\n", "Object 'GO_BP_DF' created containing assigned gene-to-BPGO mappings. To use custom mappings, generate data frame with the same structure in col 1-2.", "\n")
   affyGOCC <- eapply(get(paste(mylibbase, "GO", sep="")), getOntology, "CC") # generates list with GeneID components containing CCGOs
   GO_CC_DF <<- data.frame(GOID=unlist(affyGOCC), GeneID=rep(names(affyGOCC), as.vector(sapply(affyGOCC, length))), Count=rep(as.vector(sapply(affyGOCC, length)), as.vector(sapply(affyGOCC, length))))
   write.table(GO_CC_DF, file="GO_CC_DF", quote=F, sep="\t")
   cat("\n", "Object 'GO_CC_DF' created containing assigned gene-to-CCGO mappings. To use custom mappings, generate data frame with the same structure in col 1-2.", "\n")

   ## Generates "go_df" data frame containing the commonly used components for all GO nodes: GOID, GO Term and Ontology Type.
   require(GOstats); require(GO.db)
   go_df <- data.frame(GOID=unlist(eapply(GOTERM, function(x) x@GOID)), Term=unlist(eapply(GOTERM, function(x) x@Term)), Ont=unlist(eapply(GOTERM, function(x) x@Ontology)))
   go_df <- na.omit(go_df)
   go_df <<- go_df
   write.table(go_df, file="go_df", quote=T, sep="\t")
   cat("\n", "Object 'go_df' created containing for all GO nodes the commonly used components: GOID, GO Term and Ontology Type", "\n")
}
cat("\n", "(A.1) To use the GOHyperGAll() function, one needs 4 data frames containing the gene-to-GO mappings MF, BP, CC and the GO terms. \n       Demo data sets from GO.org or BioC can be created with one of these commands: \n \t (A.1.1) For annotations from geneontology.org: \n \t readGOorg(myfile = \"gene_association.tair\", colno = c(5,11,9), org = \"Arabidopsis\") \n \t \t myfile: download annotation table from geneontology.org and unzip it. Then point function to file name. \n \t \t colno: required column numbers; default 'c(5,11,9)' should work in most cases \n \t \t org: \"Arabidopsis\" or any string for other organisms \n \t (A.1.2) For annotations from BioC: \n \t sampleDFgene2GO(lib=\"ath1121501\") \n \t \t lib: defines annotation library, default is \"ath1121501\"", "\n")

## (A.2) Generate list containing gene-to-GO-OFFSPRING associations including assiged nodes
## This is very slow (3x3 minutes), but needs to be done only once!
gene2GOlist <- function(rootUK=T) { # If the argument 'rootUK' is set to TRUE then the root nodes are treated as terminal nodes to account for the new unknown terms
   require(GOstats); require(GO.db)
   for(i in c("MF","BP","CC")) {
      if(i=="MF") {
         go_offspr_list <- as.list(GOMFOFFSPRING) }
      if(i=="BP") {
         go_offspr_list <- as.list(GOBPOFFSPRING) }
      if(i=="CC") {
         go_offspr_list <- as.list(GOCCOFFSPRING) }
      go_offspr_list <- lapply(go_offspr_list, unlist); go_offspr_list <- lapply(go_offspr_list, as.vector) # clean-up step for the list
      go_offspr_list_temp <- lapply(names(go_offspr_list), function(x) c(x, go_offspr_list[[x]]) ) # include list component (GOID) names in corresponding (GOID) vectors
      names(go_offspr_list_temp) <- names(go_offspr_list) # names list components after go_offspr_list
      go_offspr_list <- go_offspr_list_temp
      go_offspr_list <- lapply(go_offspr_list, function(x) x[!is.na(x)]) # remove NAs in vectors

      ## Treat root nodes as terminal nodes to account for the new unknown terms. This step removes the offspring information from the root nodes.
      if(rootUK==T) {
         if(i=="MF") { go_offspr_list[["GO:0003674"]] <- c("GO:0003674") }
         if(i=="BP") { go_offspr_list[["GO:0008150"]] <- c("GO:0008150") }
         if(i=="CC") { go_offspr_list[["GO:0005575"]] <- c("GO:0005575") }
      }

      ## Retrieve gene/affy IDs for GOID vectors
      if(i=="MF") {
         MF_node_affy_list <<- lapply(go_offspr_list, function(x) unique(as.vector(GO_MF_DF[GO_MF_DF$GOID %in% x, 2])))
         save(MF_node_affy_list, file="MF_node_affy_list") }
      if(i=="BP") {
         BP_node_affy_list <<- lapply(go_offspr_list, function(x) unique(as.vector(GO_BP_DF[GO_BP_DF$GOID %in% x, 2])))
         save(BP_node_affy_list, file="BP_node_affy_list") }
      if(i=="CC") {
         CC_node_affy_list <<- lapply(go_offspr_list, function(x) unique(as.vector(GO_CC_DF[GO_CC_DF$GOID %in% x, 2])))
         save(CC_node_affy_list, file="CC_node_affy_list") }
      cat("\n", paste("Object '", i, "_node_affy_list'", sep=""), "with gene-to-GO-OFFSPRING associations created and saved in your working directory.", "\n")
   }
}
cat("\n", "(A.2) The corresponding gene-to-GO-OFFSPRING associations are created from the three data frames with the following \n       command. This creates 3 list objects with the required MF, BP and CC associations. \n \t gene2GOlist(rootUK=T)", "\n")

## (A.3) Generate AffyID-to-GeneID mappings when working with chip feature IDs
## This function creates a AffyID-to-GeneID mapping data frame using by default the TAIR mappings for the Arabidopsis ATH1 chip.
## Once the decoding data frame 'affy2locusDF' is created, the function returns for a query set of AffyIDs the corresponding GeneIDs.
## To use the function for the mappings of other chips, one needs to create the corresponding decoding data frame 'affy2locusDF'.
AffyID2GeneID <- function(map = "ftp://ftp.arabidopsis.org/home/tair/Microarrays/Affymetrix/affy_ATH1_array_elements-2008-5-29.txt", affyIDs, probe2gene=1) {
   if(!exists("affy2locusDF")) {
      cat("\n", "Downloading AffyID-to-GeneID mappings, creating object 'affy2locusDF' and saving it in your working directory", "\n")
      affy2locus <- read.delim(map, na.strings = "", fill=TRUE, header=T, sep="\t")[,-c(2:4,7:9)]
      names(affy2locus) <- c("AffyID", "AGI", "Desc")
      row.names(affy2locus) <- as.vector(affy2locus[,1])
      my_list <- apply(affy2locus[,-c(3)], 1, list); my_list <- lapply(my_list, unlist)
      my_list <- lapply(my_list, function(x) as.vector(x[-1]))
      my_list <- lapply(my_list, strsplit, ";"); my_list <- lapply(my_list, unlist)
      affy2locusDF <- data.frame(unlist(my_list))
      affy2locusDF <- data.frame(rep(names(unlist(lapply(my_list, length))), as.vector(unlist(lapply(my_list, length)))), affy2locusDF)
      names(affy2locusDF) <- c("AffyID", "GeneID")
      affy2locusDF <<- affy2locusDF
      write.table(affy2locusDF, file="affy2locusDF", quote=F, sep="\t")
   }
   if(!missing(affyIDs)) {
      if(probe2gene==1) { # For probe sets that match several loci, only the first locus ID will be used
         affy2locusDF <- affy2locusDF[!duplicated(affy2locusDF$AffyID),]
      }
      GeneIDs <- unique(as.vector(affy2locusDF[affy2locusDF[,1] %in% affyIDs, 2]))
      return(GeneIDs)
   }
}
cat("\n", "(A.3) To work with AffyIDs, the function AffyID2GeneID() can be used to import custom AffyID-to-GeneID mappings. \n \t AffyID2GeneID(map = \"ftp://ftp.arabidopsis.org/home/tair/Microarrays/Affymetrix/affy_ATH1_array_elements-2006-07-14.txt\") \n \t \t map: location of custom AffyID-to-GeneID mappings", "\n")

## (A.4) Next time things are much faster by reading the 6 data objects from file
loadData <- function() {
   need_affy2gene <- dir()
   if(length(need_affy2gene[need_affy2gene=="affy2locusDF"]) > 0) {
      affy2locusDF <<- read.table(file="affy2locusDF", header=T, colClasses = "character")
   }
   GO_MF_DF <<- read.table(file="GO_MF_DF", header=T, colClasses = "character")
   GO_BP_DF <<- read.table(file="GO_BP_DF", header=T, colClasses = "character")
   GO_CC_DF <<- read.table(file="GO_CC_DF", header=T, colClasses = "character")
   if(any(dir() %in% "go_df")) { go_df <<- read.table(file="go_df", header=T, colClasses = "character") }
}
cat("\n", "(A.4) In future R sessions one can can omit the previous 3 steps (A.1-A.3) by importing all 6 (7) data objects like this: \n \t loadData(); load(file=\"MF_node_affy_list\"); load(file=\"BP_node_affy_list\"); load(file=\"CC_node_affy_list\")", "\n")

############################
## (B) GOhyperG_All function
############################
## (B.1) Define GOhyperG_All function
GOHyperGAll <- function(gocat="MF", sample, Nannot=2) {
   ## Generates data frame (go_df) containing the commonly used components for all GO nodes: GOID, GO Term and Ontology Type. This step is only required if "go_df" hasn't been imported with the above load() function.
   require(GOstats); require(GO.db)
   if(!exists("go_df")) {
      go_df <- data.frame(GOID=unlist(eapply(GOTERM, function(x) x@GOID)), Term=unlist(eapply(GOTERM, function(x) x@Term)), Ont=unlist(eapply(GOTERM, function(x) x@Ontology)))
      go_df <<- na.omit(go_df)
      cat("\n", "Object 'go_df' created containing for all GO nodes the commonly used components: GOID, GO Term and Ontology Type", "\n")
   }
   ## (m): Obtain for every node in GO tree their number of associated genes or chip features
   if(gocat=="MF") {node_affy_list <- MF_node_affy_list}
   if(gocat=="BP") {node_affy_list <- BP_node_affy_list}
   if(gocat=="CC") {node_affy_list <- CC_node_affy_list}
   node_stats_df <- data.frame(NodeSize=sapply(node_affy_list, length))
   node_stats_df <- data.frame(GOID=row.names(node_stats_df), node_stats_df)
   row.names(node_stats_df) <- 1:length(node_stats_df[,1])
   m <- as.vector(node_stats_df$NodeSize)

   ## (x): Obtain for every node in GO tree the number of matching genes in sample set
   node_sample_stats <- sapply(node_affy_list, function(x) { sum(unlist(x) %in% sample) } )
   node_sample_stats <- as.vector(node_sample_stats)
   x <- node_sample_stats

   ## (n): Obtain the number of unique genes at GO nodes with direct annotations
   if(gocat=="MF") { GO_DF <- GO_MF_DF }
   if(gocat=="BP") { GO_DF <- GO_BP_DF }
   if(gocat=="CC") { GO_DF <- GO_CC_DF }
   n <- length(unique(GO_DF[, 2]))

   ## (k): Obtain number of unique genes in test sample that have GO mappings
   k <- length(unique(GO_DF[GO_DF[,2] %in% sample, 2]))

   ## Obtain gene/chip keys matching at GO nodes
   match_key <- sapply(node_affy_list, function(x) { x[unlist(x) %in% sample] } )
   match_key <- sapply(match_key, function(x) { paste(x, collapse=" ") } )
   match_key <- as.vector(match_key)
   key <- match_key; key[key==""] <- "NA"

   ## Apply phyper function
   phyp_v <- phyper(x-1, m, n-m , k, lower.tail = FALSE)

   ## P-value correction according to Bioinformatics, 20, 3710-3715
   Ncorrect <- table(GO_DF[GO_DF$GeneID %in% sample, 1]) # Obtain the GO nodes with direct annotations from sample set
   Ncorrect <- sum(Ncorrect >= Nannot) # Count only those that have 2 or more annotations from sample set
   if(Ncorrect<=1) {
      adj_phyp_v <- phyp_v # no adjustment necessary if Ncorrect <= 1
   } else {
      adj_phyp_v <- phyp_v * Ncorrect # Calculates simple Bonferroni correction.
      adj_phyp_v[adj_phyp_v >= 1] <- 1
      # adj_phyp_v <- sapply(phyp_v, p.adjust, method=Padj, n = Ncorrect) # Runs p.adjust(). This is disabled because most adjustment methods require that the length of the p-value vector is >= n.
   }

   ## Generate output data format
   result_df <- data.frame(node_stats_df, SampleMatch=x, Phyper=phyp_v, Padj=adj_phyp_v, SampleKeys=key)
   result_df <- merge(result_df, go_df, x.by="GOID", y.by="GOID", all.x=T)
   result_df <- result_df[order(result_df$Phyper), ]
   result_df <- result_df[,c(1:5,7:8,6)]
   result_df
}
cat("\n", "(B.1) The function GOHyperGAll() runs the phyper test against all nodes in the GO network. \n Usage: \n \t GOHyperGAll(gocat=\"MF\", sample=test_sample, Nannot=2)[1:20,] \n \t \t gocat: \"MF\", \"BP\" or \"CC\" \n \t \t Nannot: minimum number of direct annotations for p-value adjustment \n \t \t test_sample <- unique(as.vector(GO_MF_DF[1:40,2])) # for GeneIDs\n \t \t test_sample <- AffyID2GeneID(affyIDs=affy_sample, probe2gene=1) # for AffyIDs \n \t \t affy_sample <- c(\"266592_at\", \"266703_at\", \"266199_at\", \"246949_at\", \"267370_at\", \"267115_s_at\", \"266489_at\", \"259845_at\", \"266295_at\", \"262632_at\")", "\n")

####################################################################################
## (C) Subsetting of results from GOHyperGAll by assigned nodes or goSlim categories
####################################################################################
## (C.1) Define subsetting function
GOHyperGAll_Subset <- function(GOHyperGAll_result, sample=test_sample, type="goSlim", myslimv) { # type: "goSlim" or "assigned"; optional argument "myslimv" to privde custom goSlim vector
   if(type=="goSlim") {
      if(missing(myslimv)) {
         slimv <- c("GO:0003674", "GO:0008150", "GO:0005575", "GO:0030246","GO:0008289","GO:0003676","GO:0000166","GO:0019825","GO:0005515","GO:0003824","GO:0030234","GO:0003774","GO:0004871","GO:0005198","GO:0030528","GO:0045182","GO:0005215","GO:0006519","GO:0007154","GO:0016043","GO:0006412","GO:0006464","GO:0006810","GO:0007275","GO:0007049","GO:0005975","GO:0006629","GO:0006139","GO:0019748","GO:0015979","GO:0005618","GO:0005829","GO:0005783","GO:0005768","GO:0005794","GO:0005739","GO:0005777","GO:0009536","GO:0005840","GO:0005773","GO:0005764","GO:0005856","GO:0005634","GO:0005886","GO:0005576") # contains new unknown terms: "GO:0003674", "GO:0008150", "GO:0005575"
         # slimv <- c("GO:0005554", "GO:0000004", "GO:0008372", "GO:0030246","GO:0008289","GO:0003676","GO:0000166","GO:0019825","GO:0005515","GO:0003824","GO:0030234","GO:0003774","GO:0004871","GO:0005198","GO:0030528","GO:0045182","GO:0005215","GO:0006519","GO:0007154","GO:0016043","GO:0006412","GO:0006464","GO:0006810","GO:0007275","GO:0007049","GO:0005975","GO:0006629","GO:0006139","GO:0019748","GO:0015979","GO:0005618","GO:0005829","GO:0005783","GO:0005768","GO:0005794","GO:0005739","GO:0005777","GO:0009536","GO:0005840","GO:0005773","GO:0005764","GO:0005856","GO:0005634","GO:0005886","GO:0005576") # contains old unknown terms: "GO:0005554", "GO:0000004", "GO:0008372"
      } else {
         slimv <- myslimv }
      GOHyperGAll_subset <- GOHyperGAll_result[GOHyperGAll_result[,1] %in% slimv, ]
   }
   if(type=="assigned") {
      termGO <- c(as.vector(GO_MF_DF[GO_MF_DF$GeneID %in% sample, 1]),
                  as.vector(GO_BP_DF[GO_BP_DF$GeneID %in% sample, 1]),
                  as.vector(GO_CC_DF[GO_CC_DF$GeneID %in% sample, 1]))
      subset_v <- unique(termGO)
      GOHyperGAll_subset <- GOHyperGAll_result[GOHyperGAll_result[,1] %in% subset_v, ]
   }
   GOHyperGAll_subset
}
cat("\n", "(C.1) The function GOHyperGAll_Subset() allows subsetting of the GOHyperGAll() results by assigned GO nodes or custom goSlim categories.", "\n", "Usage:", "\n", "\t GOHyperGAll_result <- GOHyperGAll(gocat=\"MF\", sample=test_sample, Nannot=2)", "\n", "\t GOHyperGAll_Subset(GOHyperGAll_result, sample=test_sample, type=\"goSlim\")", "\n", "\t \t type: \"goSlim\" or \"assigned\"", "\n", "\t \t myslimv: optional argument allows usage of a custom goSlim vector", "\n")

## Apply subsetting function
# GOHyperGAll_Subset(GOHyperGAll_result, sample=test_sample, type="goSlim")

## (C.2) Plotting of subsetted results
# subset <- GOHyperGAll_Subset(GOHyperGAll_result, sample=test_sample, type="goSlim")
# pie(subset[subset$SampleMatch>0 ,3], labels=as.vector(subset[subset$SampleMatch>0 ,1]), main=unique(as.vector(subset[subset$SampleMatch>0 ,6])))
cat("\n", "(C.2) Plot pie chart of subsetted results: \n \t subset <- GOHyperGAll_Subset(GOHyperGAll_result, sample=test_sample, type=\"goSlim\") \n \t pie(subset[subset$SampleMatch>0 ,3], labels=as.vector(subset[subset$SampleMatch>0 ,1]), main=unique(as.vector(subset[subset$SampleMatch>0, 7])))", "\n")

#########################################################
## (D) Reduce GO Term Redundancy in 'GOHyperGAll_results'
#########################################################
## (D.1) The function 'GOHyperGAll_Simplify' subsets the data frame 'GOHyperGAll_result' by a user
## specified adjusted p-value cutoff and removes from it all GO nodes with overlapping children sets
## (OFFSPRING). Only the best scoring nodes remain in the data frame.
## The argument 'correct' is experimental. It aims to favor the selection of distal (information rich)
## GO terms that have at the same time a large number of sample matches. The following calculation is used
## for this adjustment: phyper x Number_of_children / SampleMatch
## Define GOHyperGAll_Simplify()
GOHyperGAll_Simplify <- function(GOHyperGAll_result, gocat="MF", cutoff=0.001, correct=T) { # gocat: "MF", "BP" or "CC"; cutoff: p-value cutoff; correct: TRUE or FALSE
   if(gocat!=as.vector(GOHyperGAll_result$Ont[!is.na(GOHyperGAll_result$Ont)])[1]) { stop("The GO categories in GOHyperGAll_Simplify() and GOHyperGAll_result need to match") }
   testDF <- GOHyperGAll_result[GOHyperGAll_result$Padj<=cutoff,]
   testDF <- data.frame(testDF, test=rep(0, times=length(testDF[,1])))
   testDF <- testDF[!is.na(testDF$Ont),]
   GOIDv <- NULL
   GO_OL_Matchv <- NULL
   while(sum(testDF$test==0)>0) {
      clusterv <- NULL
      test <- as.vector(testDF[,1])
      for(j in 1:length(test)) {
         if(gocat=="MF") { mymatch <- sum(unique(na.omit(c(test[j], as.list(GOMFOFFSPRING)[[test[j]]])) %in% na.omit(c(test[1], as.list(GOMFOFFSPRING)[[test[1]]]))))
         if(mymatch==1) { mymatch <- length(as.list(GOMFOFFSPRING)[[test[j]]])  }
         }
         if(gocat=="BP") { mymatch <- sum(unique(na.omit(c(test[j], as.list(GOBPOFFSPRING)[[test[j]]])) %in% na.omit(c(test[1], as.list(GOBPOFFSPRING)[[test[1]]]))))
         if(mymatch==1) { mymatch <- length(as.list(GOBPOFFSPRING)[[test[j]]])  }
         }
         if(gocat=="CC") { mymatch <- sum(unique(na.omit(c(test[j], as.list(GOCCOFFSPRING)[[test[j]]])) %in% na.omit(c(test[1], as.list(GOCCOFFSPRING)[[test[1]]]))))
         if(mymatch==1) { mymatch <- length(as.list(GOCCOFFSPRING)[[test[j]]])  }
         }
         clusterv <- c(clusterv, mymatch)
      }
      clusterv[clusterv==0] <- NA
      testDF <- data.frame(testDF[,-9], test=clusterv)
      if(correct==T) {
         testDF <- data.frame(testDF, decide=testDF$Padj * (testDF$test/testDF$SampleMatch))
      } else {
         testDF <- data.frame(testDF, decide=testDF$Padj) }
      GOIDv <- c(GOIDv, as.vector(testDF[order(testDF[,10]),][1,1]))
      GO_OL_Matchv <- c(GO_OL_Matchv, length(unique(unlist(strsplit(as.vector(testDF[!is.na(testDF$test),8]), " ")))))
      testDF <- testDF[is.na(testDF$test),]
      testDF <- testDF[order(testDF[,5]),-c(9,10)]
      testDF <- data.frame(testDF, test=rep(0, times=length(testDF[,1])))
      cat(GOIDv, "\n")
   }
   simplifyDF <- data.frame(GOID=GOIDv, GO_OL_Match=GO_OL_Matchv)
   simplifyDF
}

## Apply GOHyperGAll_Simplify
## simplifyDF <- GOHyperGAll_Simplify(GOHyperGAll_result, gocat="MF", cutoff=0.001, correct=T)
## data.frame(GOHyperGAll_result[GOHyperGAll_result[,1] %in% simplifyDF[,1], -8], GO_OL_Match=simplifyDF[,2])

########################################
## (D.2) Batch Analysis of Gene Clusters
########################################
## The function 'GOCluster_Report' performs the three GO analyses in batch mode: 'GOHyperGAll',
## 'GOHyperGAll_Subset' or 'GOHyperGAll_Simplify'. It processes many groups of genes (e.g.
## gene expression clusters) and organizes the results in a single data frame.
## The gene sets need to be provided in a data frame of this format:
## 	probeID/geneID	ClusterID	ClusterSize
##	id1		CL1		2
##	id2		CL1		2
##	id3		CL2		1
##	...		...		...
##
## Define 'GOCluster_Report()'

GOCluster_Report <- function(CL_DF=CL_DF, id_type="affy", method="all", CLSZ=10, cutoff=0.001, gocats=c("MF", "BP", "CC"), myslimv="default", correct=TRUE, recordSpecGO=NULL, ...) { # CLSZ: minimum cluster size; method: "all", "slim" or "simplify"; gocat: "MF", "BP" or "CC"; cutoff: adjusted p-value cutoff; recordSpecGO: argument to include one specific GOID in each of the 3 ontologies, e.g: recordSpecGO=c("GO:0003674", "GO:0008150", "GO:0005575")
   cluster_loop <- unique(as.vector(CL_DF[CL_DF[,3]>=CLSZ,2]))
   if(length(cluster_loop[grep("CL", cluster_loop)])>0) {
      cluster_loop <- paste("CL", sort(as.numeric(gsub("CL","", as.character(cluster_loop)))), sep="")
   }
   if(method=="all") {
      containerDF <- data.frame(CLID=NULL, CLSZ=NULL, GOID=NULL, NodeSize=NULL, SampleMatch=NULL, Phyper=NULL, Padj=NULL, Term=NULL, Ont=NULL, SampleKeys=NULL)
      for(i in cluster_loop) {
         cat("\n", "Processing cluster no", i, "with method: \"all\" (GOHyperGAll) \n")
         if(id_type=="affy") {
            affy_sample <- CL_DF[CL_DF[,2]==i, 1]
            test_sample <- AffyID2GeneID(affyIDs=affy_sample, probe2gene=1)
         }
         if(id_type=="gene") {
            test_sample <- as.vector(CL_DF[CL_DF[,2]==i, 1])
         }
         containerDF2 <- data.frame(GOID=NULL, NodeSize=NULL, SampleMatch=NULL, Phyper=NULL, Padj=NULL, Term=NULL, Ont=NULL, SampleKeys=NULL)
         count <- 0
         for(j in gocats) {
            count <- count+1
            GOHyperGAll_result <- GOHyperGAll(gocat=j, sample=test_sample, ...)
            tempDF <- GOHyperGAll_result[GOHyperGAll_result$Padj <= cutoff, ]
            if(length(tempDF[,1])==0) { # If filter returns empty data frame, then include at least the first two best scoring GO entries
               tempDF <- GOHyperGAll_result[1:2,]
            }
            containerDF2 <- rbind(containerDF2, tempDF)
            if(length(recordSpecGO)>0) {
               containerDF2 <- rbind(containerDF2, GOHyperGAll_result[GOHyperGAll_result[,1]==recordSpecGO[count],])
            }
            no_annot <- test_sample[!test_sample %in% eval(parse(text=paste("GO_", j, "_DF", sep="")))[,2]]
            no_annot <- no_annot[no_annot!="no_match"]
            containerDF2 <- rbind(containerDF2, data.frame(GOID=paste("no_annot_", j, sep=""), NodeSize=NA, SampleMatch=length(no_annot), Phyper=NA, Padj=NA, Term=NA, Ont=j, SampleKeys=paste(no_annot, collapse=", ")))
         }
         tempDF2 <- data.frame(CLID=rep(i, times=length(containerDF2[,1])), CLSZ=rep(unique(as.vector(CL_DF[CL_DF[,2]==i,3])), times=length(containerDF2[,1])), containerDF2)
         containerDF <- rbind(containerDF, tempDF2)
      }
      return(containerDF)
   }
   if(method=="slim") {
      containerDF <- data.frame(CLID=NULL, CLSZ=NULL, GOID=NULL, NodeSize=NULL, SampleMatch=NULL, Phyper=NULL, Padj=NULL, Term=NULL, Ont=NULL, SampleKeys=NULL)
      for(i in cluster_loop) {
         cat("\n", "Processing cluster no", i, "with method: \"slim\" (GOHyperGAll_Subset) \n")
         if(id_type=="affy") {
            affy_sample <- CL_DF[CL_DF[,2]==i, 1]
            test_sample <- AffyID2GeneID(affyIDs=affy_sample, probe2gene=1)
         }
         if(id_type=="gene") {
            test_sample <- as.vector(CL_DF[CL_DF[,2]==i, 1])
         }
         containerDF2 <- data.frame(GOID=NULL, NodeSize=NULL, SampleMatch=NULL, Phyper=NULL, Padj=NULL, Term=NULL, Ont=NULL, SampleKeys=NULL)
         count <- 0
         for(j in gocats) {
            count <- count+1
            GOHyperGAll_result <- GOHyperGAll(gocat=j, sample=test_sample, ...)
            if(any(myslimv == "default")) {
               slimv <- c("GO:0003674", "GO:0008150", "GO:0005575", "GO:0030246","GO:0008289","GO:0003676","GO:0000166","GO:0019825","GO:0005515","GO:0003824","GO:0030234","GO:0003774","GO:0004871","GO:0005198","GO:0030528","GO:0045182","GO:0005215","GO:0006519","GO:0007154","GO:0016043","GO:0006412","GO:0006464","GO:0006810","GO:0007275","GO:0007049","GO:0005975","GO:0006629","GO:0006139","GO:0019748","GO:0015979","GO:0005618","GO:0005829","GO:0005783","GO:0005768","GO:0005794","GO:0005739","GO:0005777","GO:0009536","GO:0005840","GO:0005773","GO:0005764","GO:0005856","GO:0005634","GO:0005886","GO:0005576") # contains new unknown terms: "GO:0003674", "GO:0008150", "GO:0005575"
            } else {
               slimv <- myslimv
            }
            tempDF <- GOHyperGAll_Subset(GOHyperGAll_result, sample=test_sample, type="goSlim", myslimv=slimv)
            containerDF2 <- rbind(containerDF2, tempDF)
            if(length(recordSpecGO)>0) {
               containerDF2 <- rbind(containerDF2, GOHyperGAll_result[GOHyperGAll_result[,1]==recordSpecGO[count],])
            }
            no_annot <- test_sample[!test_sample %in% eval(parse(text=paste("GO_", j, "_DF", sep="")))[,2]]
            no_annot <- no_annot[no_annot!="no_match"]
            containerDF2 <- rbind(containerDF2, data.frame(GOID=paste("no_annot_", j, sep=""), NodeSize=NA, SampleMatch=length(no_annot), Phyper=NA, Padj=NA, Term=NA, Ont=j, SampleKeys=paste(no_annot, collapse=", ")))
         }
         tempDF2 <- data.frame(CLID=rep(i, times=length(containerDF2[,1])), CLSZ=rep(unique(as.vector(CL_DF[CL_DF[,2]==i,3])), times=length(containerDF2[,1])), containerDF2)
         containerDF <- rbind(containerDF, tempDF2)
      }
      return(containerDF)
   }
   if(method=="simplify") {
      containerDF <- data.frame(CLID=NULL, CLSZ=NULL, GOID=NULL, NodeSize=NULL, SampleMatch=NULL, Phyper=NULL, Padj=NULL, Term=NULL, Ont=NULL, SampleKeys=NULL, GO_OL_Match=NULL)
      for(i in cluster_loop) {
         cat("\n", "Processing cluster no", i, "with method: \"simplify\" (GOHyperGAll_Simplify) \n")
         if(id_type=="affy") {
            affy_sample <- CL_DF[CL_DF[,2]==i, 1]
            test_sample <- AffyID2GeneID(affyIDs=affy_sample, probe2gene=1)
         }
         if(id_type=="gene") {
            test_sample <- as.vector(CL_DF[CL_DF[,2]==i, 1])
         }
         containerDF2 <- data.frame(GOID=NULL, NodeSize=NULL, SampleMatch=NULL, Phyper=NULL, Padj=NULL, Term=NULL, Ont=NULL, SampleKeys=NULL, GO_OL_Match=NULL)
         count <- 0
         for(j in gocats) {
            count <- count+1
            GOHyperGAll_result <- GOHyperGAll(gocat=j, sample=test_sample, ...)
            simplifyDF <- GOHyperGAll_Simplify(GOHyperGAll_result, gocat=j, cutoff=cutoff, correct=correct)
            if(length(simplifyDF)==0) { # If simplifyDF() returns empty data frame, then include at least the first two best scoring GO entries
               simplifyDF <- GOHyperGAll_Simplify(GOHyperGAll_result[1:2,], gocat=j, cutoff=1, correct=T)
            }
            tempDF <- data.frame(GOHyperGAll_result[GOHyperGAll_result[,1] %in% simplifyDF[,1], ], GO_OL_Match=simplifyDF[,2])
            containerDF2 <- rbind(containerDF2, tempDF)
            if(length(recordSpecGO)>0) {
               containerDF2 <- rbind(containerDF2, data.frame(GOHyperGAll_result[GOHyperGAll_result[,1]==recordSpecGO[count],], GO_OL_Match=GOHyperGAll_result[GOHyperGAll_result[,1]==recordSpecGO[count],3]))
            }
            no_annot <- test_sample[!test_sample %in% eval(parse(text=paste("GO_", j, "_DF", sep="")))[,2]]
            no_annot <- no_annot[no_annot!="no_match"]
            containerDF2 <- rbind(containerDF2, data.frame(GOID=paste("no_annot_", j, sep=""), NodeSize=NA, SampleMatch=length(no_annot), Phyper=NA, Padj=NA, Term=NA, Ont=j, SampleKeys=paste(no_annot, collapse=", "), GO_OL_Match=length(no_annot)))
         }
         tempDF2 <- data.frame(CLID=rep(i, times=length(containerDF2[,1])), CLSZ=rep(unique(as.vector(CL_DF[CL_DF[,2]==i,3])), times=length(containerDF2[,1])), containerDF2)
         containerDF <- rbind(containerDF, tempDF2)
      }
      containerDF <- containerDF[, c(1:9,11,10)]
      return(containerDF)
   }
}

## Apply GOCluster_Report
## BatchResult <- GOCluster_Report(CL_DF=CL_DF, method="all", id_type="gene", CLSZ=10, cutoff=0.001, gocats=c("MF", "BP", "CC"), recordSpecGO=c("GO:0003674", "GO:0008150", "GO:0005575"))
cat("\n", "(C.3) Batch analysis of many gene clusters: \n \t BatchResult <- GOCluster_Report(CL_DF=CL_DF, method=\"all\", id_type=\"gene\", CLSZ=10, cutoff=0.001, gocats=c(\"MF\", \"BP\", \"CC\"), recordSpecGO=c(\"GO:0003674\", \"GO:0008150\", \"GO:0005575\"))", "\n")


