library(pathview)
library(gage)
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)
# Note importing BioC pkgs after dplyr requires explicitly using dplyr::select()
library(dplyr)
library(DESeq2)
data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]
res$symbol = mapIds(org.Hs.eg.db,
                    keys=row.names(res),
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")
res$entrez = mapIds(org.Hs.eg.db,
                    keys=row.names(res),
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")
res$name =   mapIds(org.Hs.eg.db,
                    keys=row.names(res),
                    column="GENENAME",
                    keytype="ENSEMBL",
                    multiVals="first")

gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)

lapply(gobpres, head)
# Which data do you want to use? Let's use the sailfish counts.
# browseURL("http://dx.doi.org/10.6084/m9.figshare.1601975")
# countDataURL = "http://files.figshare.com/2439061/GSE37704_featurecounts.csv"
countDataURL = "http://files.figshare.com/2600373/GSE37704_sailfish_genecounts.csv"

# Import countdata
countData = read.csv(countDataURL, row.names=1) %>%
   dplyr::select(-length) %>%
   as.matrix()

# Filter data where you only have 0 or 1 read count across all samples.
countData = countData[rowSums(countData)>1, ]
head(countData)

source("https://bioconductor.org/biocLite.R")
biocLite("tximport")
library(tximport)
library(GenomicFeatures)
library(readr)

Salmon did the quantifiation of the transcript level. We want to see which genes are differentially expressed, so we need to link the transcript names to the gene names. We can use our .gtf annotation for that, and the GenomicFeatures package:

   txdb <- makeTxDbFromGFF("chr22_genes.gtf")
k <- keys(txdb, keytype = "GENEID")
tx2gene <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
head(tx2gene)

Now we can import the salmon quantification.

samples <- read.table("samples.txt", header = TRUE)
files <- file.path("quant", samples$sample, "quant.sf")
names(files) <- paste0(samples$sample)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)

Take a look at the data:

   head(txi.salmon$counts)

Differential expression using DESeq2¶

load DESeq2:

   library(DESeq2)

Instantiate the DESeqDataSet and generate result table. See ?DESeqDataSetFromTximport and ?DESeq for more information about the steps performed by the program.

dds <- DESeqDataSetFromTximport(txi.salmon, samples, ~condition)
dds <- DESeq(dds)
res <- results(dds)

Run the summary command to get an idea of how many genes are up- and downregulated between the two conditions:

   summary(res)

DESeq uses a negative binomial distribution. Such distributions have two parameters: mean and dispersion. The dispersion is a parameter describing how much the variance deviates from the mean.

You can read more about the methods used by DESeq2 in the paper or the vignette

Plot dispersions:

   plotDispEsts(dds, main="Dispersion plot")

For clustering and heatmaps, we need to log transform our data:

   rld <- rlogTransformation(dds)
head(assay(rld))

Then, we create a sample distance heatmap:

   library(RColorBrewer)
library(gplots)

(mycols <- brewer.pal(8, "Dark2")[1:length(unique(samples$condition))])
sampleDists <- as.matrix(dist(t(assay(rld))))
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[samples$condition],
          RowSideColors=mycols[samples$condition],
          margin=c(10, 10), main="Sample Distance Matrix")

We can also plot a PCA:

   DESeq2::plotPCA(rld, intgroup="condition")

It is time to look at some p-values:

   table(res$padj<0.05)
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)

Examine plot of p-values, the MA plot and the Volcano Plot:

   hist(res$pvalue, breaks=50, col="grey")
DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)

# Volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

KEGG pathway analysis¶

As always, load the necessary packages:

   library(AnnotationDbi)
library(org.Hs.eg.db)
library(pathview)
library(gage)
library(gageData)

Let’s use the mapIds function to add more columns to the results. The row.names of our results table has the Ensembl gene ID (our key), so we need to specify keytype=ENSEMBL. The column argument tells the mapIds function which information we want, and the multiVals argument tells the function what to do if there are multiple possible values for a single input value. Here we ask to just give us back the first one that occurs in the database. Let’s get the Entrez IDs, gene symbols, and full gene names.

res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
res$name <- mapIds(org.Hs.eg.db,
                   keys=row.names(res),
                   column="GENENAME",
                   keytype="ENSEMBL",
                   multiVals="first")

head(res)

We’re going to use the gage package for pathway analysis, and the pathview package to draw a pathway diagram.

The gageData package has pre-compiled databases mapping genes to KEGG pathways and GO terms for common organisms:

   data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)

Run the pathway analysis. See help on the gage function with ?gage. Specifically, you might want to try changing the value of same.dir.

foldchanges <- res$log2FoldChange
names(foldchanges) <- res$entrez
keggres <- gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
lapply(keggres, head)

Pull out the top 5 upregulated pathways, then further process that just to get the IDs. We’ll use these KEGG pathway IDs downstream for plotting. The dplyr package is required to use the pipe (%>%) construct.

library(dplyr)

# Get the pathways
keggrespathways <- data.frame(id=rownames(keggres$greater), keggres$greater) %>%
   tbl_df() %>%
   filter(row_number()<=5) %>%
   .$id %>%
   as.character()
keggrespathways

# Get the IDs.
keggresids <- substr(keggrespathways, start=1, stop=8)
keggresids

Finally, the pathview() function in the pathview package makes the plots. Let’s write a function so we can loop through and draw plots for the top 5 pathways we created above.

# Define plotting function for applying later
plot_pathway <- function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)

# Unload dplyr since it conflicts with the next line
detach("package:dplyr", unload=T)

# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp <- sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))

Thanks¶
