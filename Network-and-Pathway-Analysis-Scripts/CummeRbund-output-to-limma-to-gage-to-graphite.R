#source("https://bioconductor.org/biocLite.R")
# biocLite()

library(limma)
exprSet=read.table("ex_matrix_g1g2_h.txt", sep="\t", header=T, row.names = 1, stringsAsFactors = F)
group_list <- factor(x = c(rep("WT",6), rep("ABX",6)),levels=c("WT","ABX"))
# Set WT to be the first level
design <- model.matrix(~group_list)          # Remove the zero

v <- voom(exprSet,design,normalize="quantile", plot = T)
fit <- lmFit(v,design)
fit2 <- eBayes(fit)

tempOutput = topTable(fit2, n=Inf,adjust.method = 'BH',coef='group_listABX')
# Easier to read in 6 months than coef=2.

design <- model.matrix (~x.t$AGE+x.t$Gender+x.t$SMRIN+x.t$SMTSISCH, data=x.t)

# fit the linear model
fit <- lmFit(mydata3, design)

#apply eBayes
fit <- eBayes(fit)

#print summary results
summary(decideTests(fit))

#sort by pvalue
Topmodel=topTable(fit,n=Inf,sort="p", coef=2)

data <- read.csv(file="affy.csv", header=TRUE, sep=",", row.names =1)
res <- as.data.frame(data)

res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="SYMBOL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="SYMBOL",
                     multiVals="first")
na.omit(res)
res[complete.cases(res), ]
library(tidyr)
res %>% drop_na(entrez)

## Load required libraries
library("DESeq2");library("limma")
library("pathview");require(gage)
data(kegg.gs)

## Combine count files into dataframe
# Import data from featureCounts
countdata <- read.table("wt_CEvsRT.txt", header=TRUE, row.names=1)

# Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

# Assign condition
sampleCondition <- c("RT", "RT", "RT", "CE", "CE", "CE")

# Analysis with DESeq2 ----------------------------------------------------
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), sampleCondition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~sampleCondition)

## Run DESeq normalization
dds<-DESeq(dds)

##from GAGE

deseq2.res <- results(dds)
deseq2.fc=deseq2.res$log2FoldChange
names(deseq2.fc)=rownames(deseq2.res)
exp.fc=deseq2.fc
out.suffix="deseq2"

#get the annotation files for mouse

kg.mouse<- kegg.gsets("mouse")
kegg.gs<- kg.mouse$kg.sets[kg.mouse$sigmet.idx]

#convert gene symbol to entrez ID

gene.symbol.eg<- id2eg(ids=names(exp.fc), category='SYMBOL', org='Mm')

names(exp.fc)<- gene.symbol.eg[,2]

fc.kegg.p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
sel <- fc.kegg.p$greater[, "q.val"] < 0.2 & !is.na(fc.kegg.p$greater[, "q.val"])
path.ids <- rownames(fc.kegg.p$greater)[sel]
sel.l <- fc.kegg.p$less[, "q.val"] < 0.2 & !is.na(fc.kegg.p$less[,"q.val"])
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
require(pathview)
#view first 3 pathways as demo
pv.out.list <- sapply(path.ids2[1:3], function(pid) pathview(gene.data = exp.fc, pathway.id = pid,species = "hsa", out.suffix=out.suffix))
library(graphite);library(pathview)
graphite::
HsaReactome <- pathways()
MouseReactome <- pathways("mmusculus", "reactome")
names(MouseReactome)[1:10]
trans_path <- MouseReactome[grep("transmission",names(MouseReactome))]
trans_path
trans_path <- trans_path[[1]]
head(nodes(trans_path))
head(edges(trans_path))
trans_path
trans_path_g <- pathwayGraph(trans_path)
trans_path_g
augmented_pathway <- integrate_mir(trans_path_g, mirTarBase)
augmented_pathways_mmu_04727 <- integrate_mir(kegg_pathways_mmu_04727, mirTarBase)
augmented_pathways_mmu_04727

library(igraph)
library(network);library(sna);library(ndtv)
#Genematrix <- data.matrix(df)
m <- as.matrix(Da.Fr)
adj.m <- t(m) %*% m
diag(adj.m) <- 0
g <- network(Genematrix, directed=FALSE)
 summary(g)
plot(g)

library(pathview);library(gage);library(gageData)
data(kegg.sets.hs);data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)

kg.eco=kegg.gsets()
kg.eco.eg=kegg.gsets("eco", id.type = "entrez")
headkg.eco.eg$kg.sets, 2)
keggres = gage(p3, gsets=kg.eco.eg$kg.sets, ref = ref.idx, samp = samp.idx, same.dir = F)
lapply(keggres, head)
$greater
data(go.sets.hs)
data(go.subs.hs)
keggres = gage(p3, gsets=go.sets.hs[go.subs.hs$BP], same.dir = F)
lapply(keggres, head)


library(edgeR)
rst<-topTags(qlfTest, n=nrow(d));
drst <- data.frame("GENE_ID"=rownames(rst), rst);
write.table(drst, file="mytreatment.DEG.txt", sep="\t", row.names=F, quote=F);
## the above "drst" returns a list of differentially expressed genes with pvalues, FDR, etc., and they look alright.
# get the differentially expressed genes list, throw away genes with FDR >= 0.25
degs <- as.vector(subset(drst, FDR < 0.25)[[1]])
# do GO analysis
go <- goana(degs, species="Hs");
tgo <- topGO(go);


library(GenomicFeatures)
gff_path <- '/path/to/gff3'
makeTxDbFromGFF(file = gff_path, format = 'gff3')

biocM("stephenturner/annotables")
library(annotables)
annotables::grch38
grch38_tx2gene
grch38 %>%
   dplyr::filter(biotype == "protein_coding" & chr == "1") %>%
   dplyr::select(ensgene, symbol, chr, start, end, description) %>%
   head %>%
   knitr::kable(.)

library(DESeq2);library(airway)
data(airway)
airway <- DESeqDataSet(airway, design = ~cell + dex)
airway <- DESeq(airway)
res <- results(airway)
# tidy results with biobroom
library(biobroom)
res_tidy <- tidy.DESeqResults(res)
head(res_tidy)
res_tidy %>%
   dplyr::arrange(p.adjusted) %>%
   head(20) %>%
   dplyr::inner_join(grch38, by = c("gene" = "ensgene")) %>%
   dplyr::select(gene, estimate, p.adjusted, symbol, description) %>%
   knitr::kable(.)

library(Homo.sapiens);library(org.Hs.eg.db);library(annotate)
de <- select(Homo.sapiens, keys = degs, keytype = "SYMBOL", columns="ENTREZID")
## remove duplicate
mat <- match(de$SYMBOL, degs)
de <- de[mat,]
#edgeR code here
#in summary the final two steps are:
et <- exactTest(dge, pair=c("Control", "treated"))
etp <- topTags(et, n= 100000)
#Gage and pathview
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)
res = etp$table
foldchanges = res$logFC
names(foldchanges) = res$symbol
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
lapply(keggres, head)
keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>%
   tbl_df() %>%
   filter(row_number()<=10) %>%
   .$id %>%
   as.character()
keggrespathways
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids
data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
lapply(gobpres, head)

# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges,gene.idtype="SYMBOL", pathway.id=pid, species="hsa")
# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges,gene.idtype="SYMBOL", pathway.id=pid, species="hsa"))

select(org.Hs.eg.db, rownames(res), c("ENTREZID","SYMBOL"), "ALIAS")
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ALIAS",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ALIAS",
                     multiVals="first")

## Run DESeq normalization
countdata <- read.table("B6_vs_PPARA_count.txt", header=TRUE, row.names=1)
countdata <- countdata[ ,6:ncol(countdata)]
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))
countdata <- as.matrix(countdata)
head(countdata)
(condition <- factor(c(rep("ctl", 3), rep("exp", 3))))
library(DESeq2)
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds
dds <- DESeq(dds)

##from GAGE

deseq2.res <- results(dds)
deseq2.fc=deseq2.res$log2FoldChange
names(deseq2.fc)=rownames(deseq2.res)
exp.fc=deseq2.fc
out.suffix="deseq2"


res <- nbinomTest( cds, 'control, 'treat' )

resSig <- res[ res$padj < 0.01 & (res$log2FoldChange >1| res$log2FoldChange < -1), ]

resSig <- na.omit(resSig)

require(gage)
datakegg.gs)
deseq.fc<- resSig$log2FoldChange
names(deseq.fc)<- resSig$id
sum(is.infinite(deseq.fc))  # there are some infinite numbers, if use DESeq2, no such problem.
deseq.fc[deseq.fc>10]=10
deseq.fc[deseq.fc<-10]=-10
exp.fc<- deseq.fc

#kegg.gsets works with 3000 KEGG speicies
data(korg)
head(korg[,1:3], n=20)


#let's get the annotation files for mouse and convert the gene set to gene symbol format
                   kg.mouse<- kegg.gsets("mouse")
                   kegg.gs<- kg.mouse$kg.sets[kg.mouse$sigmet.idx]
                   lapplykegg.gs[1:3],head)

# egSymb is only for human data, so eg2sym and sym2eg functions are only for human data.
#data(egSymb)
#kegg.gs.sym<- lapplykegg.gs, eg2sym)
#lapply(kegg.gs.sym[1:3],head)

# to convert IDs among gene/transcript ID to Entrez GeneID or reverse, use eg2id and id2eg in the pathview package #written by the same person.
library(pathview)
data(bods)
bods

gene.symbol.eg<- id2eg(ids=names(exp.fc), category='SYMBOL', org='Mm') # convert the gene symbol to Entrez Gene ID
headgene.symbol.eg, n=100)
headgene.symbol.eg[,2], n=10)

names(exp.fc)<- gene.symbol.eg[,2]

fc.kegg.p<- gage(exp.fc, gsets= kegg.gs, ref=NULL, samp=NULL)
sel<- fc.kegg.p$greater[,"q.val"] < 0.1 & !is.na(fc.kegg.p$greater[,"q.val"])
table(sel)

sel.l<- fc.kegg.p$less[,"q.val"] < 0.1 & !is.na(fc.kegg.p$greater[,"q.val"])
table(sel.l)



require(gage)
datakegg.gs)

#get the annotation files for mouse

kg.mouse<- kegg.gsets("mouse")
kegg.gs<- kg.mouse$kg.sets[kg.mouse$sigmet.idx]

#convert gene symbol to entrez ID

gene.symbol.eg<- id2eg(ids=names(exp.fc), category='SYMBOL', org='Mm')

names(exp.fc)<- gene.symbol.eg[,2]

fc.kegg.p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
sel <- fc.kegg.p$greater[, "q.val"] < 0.2 & !is.na(fc.kegg.p$greater[, "q.val"])
path.ids <- rownames(fc.kegg.p$greater)[sel]
sel.l <- fc.kegg.p$less[, "q.val"] < 0.2 & !is.na(fc.kegg.p$less[,"q.val"])
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
require(pathview)
#view first 3 pathways as demo
pv.out.list <- sapply(path.ids2[1:3], function(pid) pathview(gene.data = exp.fc, pathway.id = pid,species = "mmu", out.suffix=out.suffix))


pdf("pathwayheatmap.pdf")
head(rownames(keggres) )
for (gs in rownames(keggres$greater) ){
   outname = gsub(" |:|/", "_", substr(gs, 10, 100))
   geneData(genes = kegg.gs[[gs]], exprs =keggres$greater,
            heatmap= T, Colv = F, Rowv = F, dendrogram = "none",  scatterplot = T)
}
dev.off()

But the code all is below:

   res$entrez <- mapIds(org.Hs.eg.db,
                        keys=row.names(res),
                        column="ENTREZID",
                        keytype="ALIAS",
                        multiVals="first")


foldchanges = res$logFC
names(foldchanges) = res$entrez
#========================================
#---------- Use Kegg and gage to get upregulated and downregulated pathways
#========================================

datakegg.gs)
keggres = gage(foldchanges, gsets =kegg.gs, same.dir =TRUE, compare="unpaired")
lapply(keggres, head)
write.csv(keggres,"keggres.csv")
keggrespathwaysup = data.frame(id=rownames(keggres$greater), keggres$greater) %>%
   tbl_df() %>%
   filter(row_number()<=5) %>%
   .$id %>%
   as.character()
keggrespathwaysdn = data.frame(id=rownames(keggres$less), keggres$less) %>%
   tbl_df() %>%
   filter(row_number()<=5) %>%
   .$id %>%
   as.character()

#----------------------
keggresidsup = substr(keggrespathwaysup, start=1, stop=8)
keggresidsup
keggresidsdn = substr(keggrespathwaysdn, start=1, stop=8)

gobpres = gage(foldchanges, gsets=kegg.gs, same.dir =TRUE, compare ="unpaired")

lapply(gobpres, head)
#========================================
#----------------- Define plotting function for applying later
#========================================

plot_pathway = function(pid) pathview(gene.data=foldchanges,gene.idtype="ENTREZID", pathway.id=pid, species="hsa", new.signature=FALSE)
#---------------- plot multiple pathways ups and downs
tmpup = sapply(keggresidsup, function(pid) pathview(gene.data=foldchanges,gene.idtype="ENTREZID", pathway.id=pid, species="hsa"))
tmpdn = sapply(keggresidsdn, function(pid) pathview(gene.data=foldchanges,gene.idtype="ENTREZID", pathway.id=pid, species="hsa"))

#========================================
#--------------PATHWAY HEATMAP STARTS HERE
#========================================
pdf("pathwayheatmap.pdf")
head(rownames(keggres) )
for (gs in rownames(keggres$greater) ){
   outname = gsub(" |:|/", "_", substr(gs, 10, 100))
   geneData(genes = kegg.gs[[gs]], exprs = keggres$greater,
            heatmap= T, Colv = F, Rowv = F, dendrogram = "none",  scatterplot = T)
}
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="SYMBOL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="SYMBOL",
                     multiVals="first")
foldchanges = res$logfc
names(foldchanges) = res$entrez
#---------------Use Kegg and gage to get upregulated and downregulated pathways
datakegg.gs)
keggres = gage(foldchanges, gsets =kegg.gs, same.dir =,True,  compare="unpaired" )

lapply(keggres, head)
keggrespathwaysup = data.frame(id=rownames(keggres$greater), keggres$greater) %>%
   tbl_df() %>%
   filter(row_number()<=5) %>%
   .$id %>%
   as.character()
keggrespathwaysdn = data.frame(id=rownames(keggres$less), keggres$less) %>%
   tbl_df() %>%
   filter(row_number()<=5) %>%
   .$id %>%
   as.character()
write.csv(keggrespathwaysup, "keggspathsup.csv")
write.csv(keggrespathwaysdn, "keggspathsdn.csv")
#-------------------------------------------------------
keggresidsup = substr(keggrespathwaysup, start=1, stop=8)
keggresidsup
keggresidsdn = substr(keggrespathwaysdn, start=1, stop=8)
gobpres = gage(foldchanges, gsets=kegg.gs, same.dir =TRUE, compare ="unpaired")
lapply(gobpres, head)
plot_pathway = function(pid) pathview(gene.data=foldchanges,gene.idtype="ENTREZID", pathway.id=pid,
                                      species="hsa", new.signature=FALSE)

#------------------------------plot multiple pathways ups and downs
tmpup = sapply(keggresidsup, function(pid) pathview(gene.data=foldchanges,gene.idtype="ENTREZID", pathway.id=pid, species="hsa"))
tmpdn = sapply(keggresidsdn, function(pid) pathview(gene.data=foldchanges,gene.idtype="ENTREZID", pathway.id=pid, species="hsa"))
