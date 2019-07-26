### R code from vignette source 'graphite.Rnw'
library(graph)
library(graphite)
###################################################
### code chunk number 2: startup
###################################################
pathwayDatabases()

#hsapiens biocarta
#hsapiens humancyc
#hsapiens     kegg
#hsapiens      nci
#hsapiens  panther
#hsapiens pharmgkb
#hsapiens reactome
#hsapiens    smpdb
hsanci <- pathways("hsapiens", "nci")
names(hsanci)[1:10]
p <- hsanci[[5]]
p

###################################################
### code chunk number 3: base1
###################################################
humanReactome <- pathways("hsapiens", "reactome")
names(humanReactome)[1:10]
p <- humanReactome[["ABC-family proteins mediated transport"]]
p
pathwayTitle(p)
head(nodes(p))
head(edges(p))

head(nodes(p), which = "mixed")
head(edges(p), which = "mixed")

g <- pathwayGraph(p)
g

edgeData(g)[1]

g <- pathwayGraph(p, which = "mixed")
g

edgeData(g)[1]

pSymbol <- convertIdentifiers(p, "SYMBOL")
pSymbol
head(nodes(pSymbol))

###################################################
### code chunk number 14: ident2
###################################################
reactomeSymbol <- convertIdentifiers(humanReactome[1:5], "SYMBOL")
hsanciSymbol <- convertIdentifiers(hsanci[1:5], "SYMBOL")

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

###################################################
### code chunk number 15: spia1
###################################################
library(SPIA)
data(colorectalcancer)

library(hgu133plus2.db)
top$ENTREZ <- mapIds(hgu133plus2.db,
                     top$ID, "ENTREZID", "PROBEID", multiVals = "first")
top <- top[!is.na(top$ENTREZ) & !duplicated(top$ENTREZ), ]
head(top)
top$ENTREZ <- paste("ENTREZID", top$ENTREZ, sep = ":")
tg1 <- top[top$adj.P.Val < 0.05, ]
head(tg1)

DE_Colorectal = tg1$logFC
names(DE_Colorectal) <- tg1$ENTREZ
ALL_Colorectal <- top$ENTREZ
head(DE_Colorectal)
head(ALL_Colorectal)

biocarta <- pathways("hsapiens", "biocarta")[1:10]
biocarta <- convertIdentifiers(biocarta, "ENTREZID")
prepareSPIA(biocarta, "biocartaEx")
res <- runSPIA(de = DE_Colorectal, all = ALL_Colorectal, "biocartaEx")
res[1:5,]


###################################################
### code chunk number 16: tg1
###################################################
library(topologyGSA)
data(examples)
colnames(y1) <- paste("SYMBOL", colnames(y1), sep = ":")
colnames(y2) <- paste("SYMBOL", colnames(y2), sep = ":")

kegg <- pathways("hsapiens", "kegg")
p <- convertIdentifiers(kegg[["Fc epsilon RI signaling pathway"]], "SYMBOL")
runTopologyGSA(p, "var", y1, y2, 0.05)


###################################################
### code chunk number 17: tg2
###################################################
library(ALL)
library(DEGseq)
library(a4Preproc)
library(clipper)

data(ALL)
pheno <- as(phenoData(ALL), "data.frame")
samples <- unlist(lapply(c("NEG", "BCR/ABL"), function(t) {
  which(grepl("^B\\d*", pheno$BT) & (pheno$mol.biol == t))[1:10]
}))
classes <- c(rep(1,10), rep(2,10))

expr <- exprs(ALL)[,samples]
rownames(expr) <- paste("ENTREZID",
                        featureData(addGeneInfo(ALL))$ENTREZID,
                        sep = ":")

k <- as.list(pathways("hsapiens", "kegg"))
selected <- k[c("Chronic myeloid leukemia",
                "Bladder cancer",
                "Cytosolic DNA-sensing pathway")]

clipped <- runClipper(selected, expr, classes, "mean", pathThr = 0.1)
resClip <- do.call(rbind, clipped$results)
resClip[, c("startIdx", "endIdx", "maxIdx", "lenght",
           "maxScore", "aScore", "involvedGenes")]


###################################################
### code chunk number 18: bldp
###################################################
edges <- data.frame(src_type = "ENTREZID", src="672",
                    dest_type = "ENTREZID", dest="7157",
                    direction="undirected", type="binding")
pathway <- buildPathway("#1", "example", "hsapiens", "database",
                        proteinEdges = edges)


###################################################
### code chunk number 19: bldp2
###################################################
edges <- data.frame(src_type = "ENTREZID", src="672",
                    dest_type = "ENTREZID", dest="7157",
                    direction="undirected", type="binding")
edgemix <- data.frame(src_type = "CHEBI", src="77750",
                    dest_type = "ENTREZID", dest="7157",
                    direction="undirected", type="biochemicalReaction")
edgemet <- data.frame(src_type = "CHEBI", src="15351",
                    dest_type = "CHEBI", dest="77750",
                    direction="undirected", type="biochemicalReaction")
pathway <- buildPathway("#1", "example", "hsapiens", "database",
                        proteinEdges = edges,
                        mixedEdges = edgemix,
                        metaboliteEdges = edgemet)


###################################################
### code chunk number 20: para1 (eval = FALSE)
###################################################
## options(Ncpus = 6)


###################################################
### code chunk number 21: para1 (eval = FALSE)
###################################################
## original <- pathways("hsapiens", "reactome")
##
## # Do (will exploit parallelism)
## converted <- convertIdentifiers(original, "SYMBOL")
##
## # Don't (no parallelism here)
## converted <- lapply(original, convertIdentifiers, "SYMBOL")


