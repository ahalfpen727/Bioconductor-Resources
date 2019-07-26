```{r}
#devtools::install_github("genomicsclass/ph525x")
library(ph525x)
??ph525x
library(rete)
rete::

source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R")
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/dendroCol.R")

```

```{r import-drosophila-gene-expr,  message=FALSE, warning=FALSE, include=FALSE}}
genome="/home/drew/umb_triley/Reference-Genomes/Drosophila/dmel_b5.41_ncbi/genome.fa"
gtfFile="/home/drew/umb_triley/Reference-Genomes/Drosophila/dmel_b5.41_ncbi/genes.gtf"
cuff.a<- readCufflinks(dir='/home/drew/umb_triley/drosophila/AnalysisOne/cuffdiff_results_ncbi_b5.41_gtf_guided/withVitA-over-woVitA/',genome="/home/drew/umb_triley/Reference-Genomes/Drosophila/dmel_b5.41_ncbi/genome.fa",gtfFile='/home/drew/umb_triley/drosophila/AnalysisOne/cuffcompare_results_ncbi_b5.41_gtf_guided/cuffcmp.combined.gtf', rebuild=F)

cuff.b<- readCufflinks(dir="/home/drew/umb_triley/urine1/cuffdiff_results_hg19_default/LUTS-over-CTRL",genome="/home/drew/umb_triley/ReferenceGenomes/Human/USCS_hg19/genome.fa",gtfFile="/home/drew/umb_triley/ReferenceGenomes/Human/USCS_hg19/genes.gtf", rebuild=F)

cuff.c<- readCufflinks(dir=inDir,genome="/home/drew/umb_triley/ReferenceGenomes/Human/USCS_hg19/genome.fa",gtfFile="/home/drew/umb_triley/urine1/cuffcompare_results_hg19_gtf_guided/cuffcmp.combined.gtf", rebuild=F)

library(gsean)
library(TCGAbiolinks)
library(DESeq2)
library(PPInfer)
# TCGA LUAD
query <- GDCquery(data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  platform = "Illumina HiSeq", 
                  file.type  = "normalized_results",
                  experimental.strategy = "RNA-Seq",
                  legacy = TRUE)

query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  platform = "Illumina HiSeq", 
                  file.type  = "normalized_results",
                  experimental.strategy = "RNA-Seq",
                  legacy = TRUE)
GDCdownload(query, method = "api")
invisible(capture.output(data <- GDCprepare(query)))
exprs.LUAD <- assay(data)
# remove duplicated gene names
exprs.LUAD <- exprs.LUAD[-which(duplicated(rownames(exprs.LUAD))),]
# list of genes
recur.mut.gene <- c("KRAS", "TP53", "STK11", "RBM10", "SPATA31C1", "KRTAP4-11",
                    "DCAF8L2", "AGAP6", "KEAP1", "SETD2", "ZNF679", "FSCB",
                    "BRAF", "ZNF770", "U2AF1", "SMARCA4", "HRNR", "EGFR")

# KEGG_hsa
load(system.file("data", "KEGG_hsa.rda", package = "gsean"))

# GSEA
set.seed(1)
result.GSEA <- gsean(KEGG_hsa, recur.mut.gene, exprs.LUAD, threshold = 0.8)
invisible(capture.output(p <- GSEA.barplot(result.GSEA, category = 'pathway',
                                           score = 'NES', pvalue = 'padj',
                                           sort = 'padj', top = 20)))
p <- p + scale_fill_gradient(low = "red", high = "blue") +
  theme(plot.margin = margin(10, 10, 10, 75))
plotly::ggplotly(p + guides(fill = FALSE))


2.2 GSEA with statistics, based on the degree centrality

Gene expression data are analyzed to identify differentially expressed genes. Consider pasilla RNA-seq count data. By using the Wald statistic in this example, GSEA can be performed with Gene Ontology terms from http://www.go2msig.org/cgi-bin/prebuilt.cgi?taxid=7227. Thus, it is expected that we may find the biological functions related to change in phenotype from the network, rather than individual genes.

library(gsean)
library(pasilla)
library(DESeq2)
library(PPInfer)
# pasilla count data
pasCts <- system.file("extdata", "pasilla_gene_counts.tsv",
                      package = "pasilla", mustWork = TRUE)
cts <- as.matrix(read.csv(pasCts, sep="\t", row.names = "gene_id"))
condition <- factor(c(rep("untreated", 4), rep("treated", 3)))
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData   = data.frame(condition),
  design    = ~ 0 + condition)
# filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# differentially expressed genes
dds <- DESeq(dds)
resultsNames(dds)

## [1] "conditiontreated"   "conditionuntreated"

res <- results(dds, contrast = list("conditiontreated", "conditionuntreated"), listValues = c(1, -1))
statistic <- res$stat
names(statistic) <- rownames(res)
exprs.pasilla <- counts(dds, normalized = TRUE)

# convert gene id
library(org.Dm.eg.db)
gene.id <- AnnotationDbi::select(org.Dm.eg.db, names(statistic), "ENTREZID", "FLYBASE")
names(statistic) <- gene.id[,2]
rownames(exprs.pasilla) <- gene.id[,2]

# GO_dme
load(system.file("data", "GO_dme.rda", package = "gsean"))

# GSEA
set.seed(1)
result.GSEA <- gsean(GO_dme, statistic, exprs.pasilla)

## log2 transformed for correlation

invisible(capture.output(p <- GSEA.barplot(result.GSEA, category = 'pathway',
                                           score = 'NES', top = 50, pvalue = 'padj',
                                           sort = 'padj', decreasing = FALSE,
                                           numChar = 110)))
p <- p + scale_fill_continuous(low = 'red', high = 'green') +
  theme(plot.margin = margin(10, 10, 10, 50))
plotly::ggplotly(p)


# Experiment-Level Gene Expression Diff-Table
The *CummeRbund* Bioconductor package provides functions that facilitate visualizing and annotating the differential expression analysis for all feature levels. A funtion that consolidates the various gene expression feature levels into a single table is provided. The command is computationally intensive and should be run alone without any other processes that represent a computational load.

```{r difftable}
options(mc.cores = parallel::detectCores())
# This command is computationally intensive
all_exp.df<-diffTable(cummeRbund::genes(cuff))
```

 
# GTF-guided Feature Analysis
Tailor employs the *Cufflinks* tool-suite to assemble transcripts and identify potentialy novel features. *Cufflinks* uses a whole genome fasta along with an annotation file that indexes the genome to create a de-bruijn graphs via bipartite matching of the isoforms into the minimum features necessary to explain the paired-end reads. Depending on the type of annotations used initially additional annotation formating may be necessary.


```{r fix-annotation-in-de-novo-analysis}
# gene expr data with all annotated features
gene_exp.df<-diffData(cummeRbund::genes(cuff))
sig_genes.df<-subset(gene_exp.df, gene_exp.df$significant=="yes")
all_cols<-colnames(all_exp.df)
# reformat all_exp
st<-grep(pattern="status", all_cols)
all_exp.df$status<-all_exp.df[,st]
si<-grep(pattern="significant", all_cols)
all_exp.df$significant<-all_exp.df[,si]
p<-grep(pattern="p_value", all_cols)
all_exp.df$p_value<-all_exp.df[,p]
q<-grep(pattern="q_value", all_cols)
all_exp.df$q_value<-all_exp.df[,q]
v1<-grep(pattern="value_1", all_cols)
all_exp.df$value_1<-all_exp.df[,v1]
v2<-grep(pattern="value_2", all_cols)
all_exp.df$value_2<-all_exp.df[,v2]
log<-grep(pattern="log2", all_cols)
all_exp.df$log2_fold_change<-all_exp.df[,log]
te<-grep(pattern="test_stat", all_cols)
all_exp.df$test_stat<-all_exp.df[,te]
o.u<-grep(pattern=over, all_cols)
all_exp.df<-all_exp.df[,-o.u]

all.exp.tib<-all_exp.df %>%
   mutate(length= all_exp.df$width) %>%
   mutate(XLOC.id= all_exp.df$gene_id..28) %>%
   mutate(ref.seq.id = all_exp.df$nearest_ref_id) %>%
   mutate(chrs= all_exp.df$seqnames) %>%
  select(-c(seqnames, score, gene_id..8,nearest_ref_id,
             gene_id..28, coverage, width, class_code..34))

all.exp.df <-all_exp.df %>%
   group_by(gene_id)
   
table(all.exp.df$seqnames)
table(all.exp.df$strand)

all.sig.exp.df <-gene_exp.df %>%
   filter(significant=="yes") %>%
   group_by(gene_id) %>% 
   ungroup()
str(all.sig.exp.df)

iso.rep.matrix<-as.data.frame(repFpkmMatrix(isoforms(cuff)))
head(iso.rep.matrix)
gene.features<-annotation(genes(cuff))
head(gene.features)
table(gene.features$type)

gfeats<-gene.features %>%
   group_by(type) %>%
   summarise(n = n())
gfeats   
gfeats.chr<-grep(pattern="_",gene.features$seqnames, invert=T)
gfeats.p<-gene.features[gfeats.chr,]

```



```{r adj-matrix-non-siggenes, message=FALSE, warning=FALSE, include=FALSE}
over_gene_data<-subset(gene_exp.df, (log2_fold_change > 0)) #& (significant =='yes') 
head(over_gene_data);dim(over_gene_data)
under_gene_data<-subset(gene_exp.df, (log2_fold_change < 0)) #& (significant =='yes') 
head(under_gene_data);dim(under_gene_data)

under.genes<-row.names(g.rep.ma) %in% under_gene_data$gene_id
u.rep.ma<-g.rep.ma[under.genes,under.group]
head(u.rep.ma);dim(u.rep.ma)

over.genes<-row.names(g.rep.ma) %in% over_gene_data$gene_id
o.rep.ma<-g.rep.ma[over.genes,over.group]
head(o.rep.ma);dim(o.rep.ma)

u.rep.sd<-genefilter::rowSds(u.rep.ma)
o.rep.sd<-genefilter::rowSds(o.rep.ma)
summary(u.rep.sd)
summary(o.rep.sd)
u.rep.mu<-rowMeans(u.rep.ma)
o.rep.mu<-rowMeans(o.rep.ma)
summary(u.rep.mu)
summary(o.rep.mu)

o_gene_data<-subset(o.rep.ma, (o.rep.mu > o.rep.sd)) #& (significant =='yes') 
u_gene_data<-subset(u.rep.ma, (u.rep.mu > u.rep.sd)) #& (significant =='yes') 

o.graph <- graph.adjacency(as.matrix(as.dist(cor(t(o_gene_data), method="pearson"))),
                              mode="undirected", weighted=TRUE, diag=FALSE)
u.graph <- graph.adjacency(as.matrix(as.dist(cor(t(u_gene_data), method="pearson"))),
                              mode="undirected", weighted=TRUE, diag=FALSE)

```


```{r tailor-power message=FALSE, warning=FALSE, include=FALSE}
g.f.cnt.sds <- rowSds(as.matrix(g.f.cnt.ma))
effect.size<-sqrt((g.u.cnt.mu - g.o.cnt.mu)^2)/g.f.cnt.sds
#effect.size
g.cnt.fstats = rowFtests(as.matrix(g.f.cnt.ma), as.factor(groups))
g.cnt.s.fstats<-subset(g.cnt.fstats, p.value < 0.05)
ftest.sigs<-row.names(g.cnt.s.fstats)
number.of.tests<-nrow(g.cnt.fstats)
num.sig.pvals<-sum(g.cnt.fstats$p.value < 0.05)
true.sig<-c(num.pval *0.05)
sig.calls<-c((num.pval *0.05)-num.sig.pvals)
prop.non.de<-(number.of.tests-num.sig.pvals)/number.of.tests
typeIIerrors<-true.sig-sig.calls
typeIIerr.rate<-c(true.sig-sig.calls)/number.of.tests
power=1-typeIIerr.rate
library(ssizeRNA)
size1 <- ssizeRNA_single(nGenes = number.of.tests, pi0 = 0.8, m = 200, mu = 10,
disp = 0.1, fc = 2, fdr = 0.05,
power = 0.8, maxN = 20)

number.of.tests
num.sig.pvals
true.sig
prop.non.de
sig.calls
typeIIerrors
typeIIerr.rate
power


# Power Analysis Produced from Cuffdiff Output
Using the difference between the number of significantly differentially expressed genes identified from p-values below $\alpha=0.05$ and the FDR corrected p-values below $\alpha=0.05$ the type II error rate can be calculated and from it the statistical power of the analysis can be extrapolated.
```


```{r correlelogram}
set.seed(955)
vvar <- 1:20 + rnorm(20,sd=3)
wvar <- 1:20 + rnorm(20,sd=5)
xvar <- 20:1 + rnorm(20,sd=3)
yvar <- (1:20)/2 + rnorm(20, sd=10)
zvar <- rnorm(20, sd=6)
# A data frame with multiple variables
data <- data.frame(vvar, wvar, xvar, yvar, zvar)
head(data)
#To make the graph:
library(ellipse)
# Make the correlation table
ctab <- cor(data)
round(ctab, 2)
# Make the graph, with reduced margins
plotcorr(ctab, mar = c(0.1, 0.1, 0.1, 0.1))
# Do the same, but with colors corresponding to value
colorfun <- colorRamp(c("#CC0000","white","#3366CC"), space="Lab")
plotcorr(ctab, col=rgb(colorfun((ctab+1)/2), maxColorValue=255),
         mar = c(0.1, 0.1, 0.1, 0.1))
```

