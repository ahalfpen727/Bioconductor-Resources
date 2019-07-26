## ---- message=FALSE------------------------------------------------------
source("https://bioconductor.org/biocLite.R")
biocLite("QuaternaryProd")
library(QuaternaryProd)

# Compute the probability mass function
pmf <- QP_Pmf(q_p = 20, q_m = 20, q_z = 20, q_r = 0, n_p = 20, n_m = 20, n_z = 20)

# Plot the mass function
plot(names(pmf), pmf, col="blue", xlab = "scores", ylab = "probabilities")
lines(names(pmf), pmf, col = "blue")

## ---- message=FALSE------------------------------------------------------
# Get the p-value of score 5
pval <- QP_Pvalue(score = 5, q_p = 20, q_m = 20, q_z = 20, q_r = 0,
                                                     n_p = 20, n_m = 20, n_z = 20)
pval

# Compue the p-value only if it is statistically significant otherwise
# return -1
pval <- QP_SigPvalue(score = 5, q_p = 20, q_m = 20, q_z = 20, q_r = 0,
                                                     n_p = 20, n_m = 20, n_z = 20)
pval

## ---- message=FALSE------------------------------------------------------
library(QuaternaryProd)
library(pathview);library(gage);library(gageData)
data(kegg.sets.hs);data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)
# Note importing BioC pkgs after dplyr requires explicitly using dplyr::select()
library(dplyr);library(DESeq2); library(org.Hs.eg.db)
data(go.sets.hs);data(go.subs.hs)

#siggenes.df<-sig_genes_exp.diff[,"gene_id","log2_fold_change", "p_value"]
siggenes.df<-sig_genes_exp.diff[,c(1,7,10)]
head(siggenes.df)
diff.table = file.path(DiffTable)
write.table(isodiff, file = diff.table, sep = "  ", row.names = F , col.names = T,quote = F)


gobpsets = go.sets.hs[go.subs.hs$BP]
columns(org.Hs.eg.db)
siggenes.df$entrezid = mapIds(org.Hs.eg.db,
                              keys=siggenes.df[,1],
                              column="ENTREZID",
                              keytype="SYMBOL",
                              multiVals="first")

length(siggenes.df$entrezid)
length(unique(siggenes.df$entrezid))
sig.genes.df<-subset(siggenes.df, siggenes.df$entrezid  %in%  unique(siggenes.df$entrezid))
SigGenes.df<-na.omit(sig.genes.df)
as.integer(siggenes.df$entrezid)
SigGenes.df<-cbind(entrezid=SigGenes.df[,"entrezid"], SigGenes.df)
SigGenes.df<-SigGenes.df[,-c(2,5)]
head(SigGenes.df); dim(SigGenes.df)
colnames(SigGenes.df)<-c("entrez", "fc", "pvalue")
# Get gene expression data
e2f3 <- system.file("extdata", "e2f3_sig.txt",
                             package = "QuaternaryProd")
e2f3 <- read.table(e2f3, sep = "\t",
                             header = TRUE, stringsAsFactors = FALSE)
myc <- system.file("extdata", "myc_sig.txt",
                             package = "QuaternaryProd")
myc <- read.table(myc, sep = "\t",
                             header = TRUE, stringsAsFactors = FALSE)

# Rename column names appropriately
# and remove duplicated entrez ids in the gene expression data
names(e2f3) <- c("entrez", "pvalue", "fc")
e2f3 <- e2f3[!duplicated(e2f3$entrez),]

names(myc) <- c("entrez", "pvalue", "fc")
myc <- myc[!duplicated(myc$entrez),]

## ---- message=FALSE------------------------------------------------------
# Compute the Quaternary Dot Product Scoring Statistic for only statistically
# significant regulators
quaternary_results <- RunCRE_HSAStringDB(SigGenes.df, method = "Quaternary",
                                   fc.thresh = log2(1.3), pval.thresh = 0.05,
                                   only.significant.pvalues = TRUE,
                                   significance.level = 0.05) #,
#                                   epsilon = 1e-16)
quaternary_results[1:4, c("uid","symbol","regulation","pvalue")]

## ---- message=FALSE------------------------------------------------------
ternary_results <- RunCRE_HSAStringDB(myc, method = "Ternary",
                                      fc.thresh = log2(1.3), pval.thresh = 0.05,
                                      only.significant.pvalues = TRUE,
                                      significance.level = 0.05,
                                      epsilon = 1e-16)
ternary_results[1:4, c("uid","symbol","regulation","pvalue")]

## ---- message=FALSE------------------------------------------------------
enrichment_results <- RunCRE_HSAStringDB(myc, method = "Enrichment",
                                         fc.thresh = log2(1.3), pval.thresh = 0.05,
                                         only.significant.pvalues = TRUE,
                                         significance.level = 0.05,
                                         epsilon = 1e-16)
enrichment_results[1:10, c("uid","symbol","regulation","pvalue")]

