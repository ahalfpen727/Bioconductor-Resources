BiocManager::install("cummeRbund")
BiocManager::install("limma")
BiocManager::install("DESeq")
BiocManager::install("GeneRfold")
BiocManager::install("GeneRfold")

BiocManager::install("cummeRbund", version = "3.8")

BiocManager::install("limma")

if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("featurecounts", version = "3.5")

if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("edgeR", version = "3.5")

if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db", version = "3.5")

if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("ballgown", version = "3.5")

if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("biomaRt", version = "3.5")


if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("IRanges", version = "3.5")
BiocManager::install("AnnotationHub", version = "3.8")
## ----knitr, echo=FALSE, results="hide"-----------------------------------
library("cummeRbund")
library("limma")
library("knitr")
library("org.Hs.eg.db")
library("markdown")

library("org.Hs.eg.db")
opts_chunk$set(tidy=FALSE,tidy.opts=list(width.cutoff=30),dev="png",fig.show="hide",
               fig.width=4,fig.height=4.5,
               message=FALSE)

## ----style, eval=TRUE, echo=FALSE, results="asis"--------------------------
BiocStyle::latex()
version
## ----options, results="hide", echo=FALSE--------------------------------------
options(digits=3, width=80, prompt=" ", continue=" ")

## ----install_cellTree, eval=FALSE---------------------------------------------
 if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
  BiocManager::install("cellTree")

version
## ----install_missing_bioconductor_packages, eval=FALSE------------------------
#  BiocManager::install(c("HSMMSingleCell", "org.Hs.eg.db", "biomaRt"))

## ----init_sincell, cache=FALSE, eval=TRUE,warning=FALSE-----------------------
library(cellTree)

## ----load_hsmm_data, eval=TRUE------------------------------------------------
# load HSMMSingleCell package and load the data set:
library(HSMMSingleCell)
data(HSMM_expr_matrix)

# Total number of genes * cells:
dim(HSMM_expr_matrix)

## ----compute_lda_maptpx, eval=FALSE-------------------------------------------
#  # Run LDA inference using 'maptpx' method
#  # finding best number of topics k between 3 and 8:
#  lda.results = compute.lda(HSMM_expr_matrix, k.topics=3:8, method="maptpx")

## ----compute_lda_gibbs, eval=FALSE--------------------------------------------
#  # Run LDA inference using 'Gibbs' method for k = 6 topics:
#  lda.results = compute.lda(HSMM_expr_matrix, k.topics=6, method="Gibbs")

## ----compute_lda_with_hgnc, eval=FALSE----------------------------------------
#  HSMM_expr_matrix.hgnc = HSMM_expr_matrix
#
#  library("biomaRt")
#  ensembl.ids = sapply(strsplit(rownames(HSMM_expr_matrix), split=".",fixed=TRUE),
#  					 "[",
#  					 1)
#  ensembl.mart = useMart(host="www.ensembl.org",
#  					   "ENSEMBL_MART_ENSEMBL",
#  					   dataset = "hsapiens_gene_ensembl")
#  gene.map = getBM(attributes = c("ensembl_gene_id", "entrezgene", "hgnc_symbol"),
#  				 filters = "ensembl_gene_id",
#  				 values = ensembl.ids,
#  				 mart = ensembl.mart)
#  idx = match(ensembl.ids, gene.map$ensembl_gene_id)
#  hgnc.ids = gene.map$hgnc_symbol[idx]
#  has.hgnc.ids = !is.na(hgnc.ids)&(hgnc.ids!="")
#  rownames(HSMM_expr_matrix.hgnc)[has.hgnc.ids] = hgnc.ids[has.hgnc.ids]
#
#  HSMM_lda_model = compute.lda(HSMM_expr_matrix.hgnc, k.topics=6)

## ----load_lda_with_hgnc, eval=TRUE--------------------------------------------
# Load pre-computed LDA model for skeletal myoblast RNA-Seq data
# from HSMMSingleCell package:
data(HSMM_lda_model)

# Number of topics of fitted model:
print(HSMM_lda_model$K)

# Model uses HGCN gene names:
head(rownames(HSMM_lda_model$theta))

## ----pairwise_distances, eval=TRUE--------------------------------------------
# Compute pairwise distance between cells
# based on topic distributions in the fitted model:
dists = get.cell.dists(HSMM_lda_model)

print(dists[1:5,1:5])

## ----day_annotation, eval=TRUE------------------------------------------------
# Recover sampling time point for each cell:
library(HSMMSingleCell)
data(HSMM_sample_sheet)
days.factor = HSMM_sample_sheet$Hours
days = as.numeric(levels(days.factor))[days.factor]

# Our grouping annotation (in hours):
print(unique(days))

## ----mst_tree, eval=TRUE------------------------------------------------------
# compute MST from a fitted LDA model:
mst.tree = compute.backbone.tree(HSMM_lda_model, days, only.mst=TRUE)

## ----plot_mst_tree_with_topics, eval=TRUE, echo=TRUE, fig.show="asis", dpi=144, fig.width=5, fig.height=5, out.width="5in", out.height="5in"----
# plot the tree (showing topic distribution for each cell):
mst.tree.with.layout = ct.plot.topics(mst.tree)

## ----plot_mst_tree_with_groups, echo=TRUE, fig.show="asis", dpi=144, fig.width=5, fig.height=5, out.width="5in", out.height="5in"----
# plot the tree (showing time point for each cell):
mst.tree.with.layout = ct.plot.grouping(mst.tree)

## ----plot_btree_with_groups, eval=TRUE, echo=TRUE, fig.show="asis", dpi=144, fig.width=5, fig.height=5, out.width="5in", out.height="5in"----
# compute backbone tree from a fitted LDA model:
b.tree = compute.backbone.tree(HSMM_lda_model, days)

# plot the tree (showing time label for each cell):
b.tree.with.layout = ct.plot.grouping(b.tree)

## ----plot_btree_with_groups_wider, eval=TRUE, echo=TRUE, fig.show="asis", dpi=144, fig.width=5, fig.height=5, out.width="5in", out.height="5in"----
# compute backbone tree from a fitted LDA model:
b.tree = compute.backbone.tree(HSMM_lda_model, days, width.scale.factor=1.5)

# plot the tree (showing time label for each cell):
b.tree.with.layout = ct.plot.grouping(b.tree)

## ----plot_btree_with_topics_wider, eval=TRUE, echo=TRUE, fig.show="asis", dpi=144, fig.width=5, fig.height=5, out.width="5in", out.height="5in"----
# plot the tree (showing topic distribution for each cell):
b.tree.with.layout = ct.plot.topics(b.tree)

## ----go_lib, eval=TRUE--------------------------------------------------------
# Load GO mappings for human:
library(org.Hs.eg.db)

## ----go_terms, eval=TRUE, results="hide"--------------------------------------
# Compute GO enrichment sets (using the Cellular Components category)
# for each topic
go.results = compute.go.enrichment(HSMM_lda_model,
                                   org.Hs.eg.db, ontology.type="CC",
                                   bonferroni.correct=TRUE, p.val.threshold=0.01)

## ----go_terms_print, eval=TRUE------------------------------------------------
# Print ranked table of significantly enriched terms for topic 1
# that do not appear in other topics:
go.results$unique[[1]]

## ----go_terms_dag_files, eval=FALSE-------------------------------------------
#  # Compute GO enrichment sets (using the Biological Process category)
#  # for each topic and saves DAG plots to files:
#  go.results.bp = compute.go.enrichment(HSMM_lda_model,
#                                  org.Hs.eg.db, ontology.type="BP",
#                                  bonferroni.correct=TRUE, p.val.threshold=0.01,
#                                  dag.file.prefix="hsmm_go_")

## ----plot_go_results, eval=TRUE, echo=TRUE, fig.show="asis", dpi=144, fig.width=5, fig.height=5, out.width="5in", out.height="5in"----
# plot GO sub-DAG for topics 1 to 3:
go.dag.subtree = ct.plot.go.dag(go.results,
                                up.generations = 2,
                                only.topics=c(1:3))

## ----cell_ordering_table, eval=TRUE-------------------------------------------
# Generate table summary of cells, ranked by tree position:
cell.table = cell.ordering.table(b.tree)

# Print first 5 cells:
cell.table[1:5,]

## ----cell_ordering_table_latex, eval=FALSE------------------------------------
#  # Generate table summary of cells, ranked by tree position:
#  cell.table = cell.ordering.table(b.tree,
#                                  write.to.tex.file="cell_summary.tex")

## ----session_info, eval=TRUE--------------------------------------------------
sessionInfo()

