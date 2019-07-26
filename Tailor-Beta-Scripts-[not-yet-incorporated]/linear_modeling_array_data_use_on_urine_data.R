############################################
##          Package Installation          ##
############################################
source("http://bioconductor.org/biocLite.R")
 biocLite("ath1121501probe")
# previously_installed <- row.names(installed.packages())  
# current <- c("affyio", "affycoretools", 
# "GOstats", "Ruuid", "graph", "GO.db", "Category", 
# "plier", "simpleaffy", "ath1121501.db", "ath1121501cdf",
#  "ath1121501probe", "ricecdf", "riceprobe", "biomaRt")
# previously_installed[!previously_installed %in% current]
############################################

#######################################################
## (A) Download a sample set of Affymetrix cel files ##
#######################################################
## Right-click this link 
#(http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/Workshop.zip) 
## and save its content to your computer: Workshop.zip. 
## These sample CEL files are from the GEO data set: GSE5621.
##  After unpacking this file archive, see six *.cel files.

## (A.1) Generate RMA expression data, 
## MAS5 P/M/A calls and export results to Excel

library(affy)
library(ath1121501probe)
mydata <- ReadAffy() 
# Reads cel files in current working directory into affybatch object 'mydata'. 
print.data.frame(ath1121501probe[1:22,]) 
eset_rma <- rma(mydata) 
# Generates RMA expression values and stores them as ExpressionSet.
exprs(eset_rma)[1:4,] 
# Prints first 4 rows in data frame structure.
# exprs(eset_rma)<- 2^(exprs(eset_rma))[1:4,] 
# If need to unlog RMA expression values.
mydf <- 2^exprs(eset_rma)
myList <- tapply(colnames(mydf), c(1,1,2,2,3,3), list)
names(myList) <- sapply(myList, paste, collapse="_")
mymean <- sapply(myList, function(x) rowMeans(mydf[,x]))
# calculating mean values for any sample combination.
eset_pma <- mas5calls(mydata) 
# Generates MAS 5.0 P/M/A calls.
my_frame <- data.frame(exprs(eset_rma), exprs(eset_pma), 
					   assayDataElement(eset_pma, "se.exprs"))
# Combine RMA intensities, P/M/A calls plus their wilcoxon p-values
my_frame <- my_frame[, sort(names(my_frame))] 
# Sorts columns by cel file name.
write.table(my_frame, file="my_file.xls", sep="\t", col.names = NA) 
# Exports data to text file that can be imported into Excel.

## Add annotation information
library("ath1121501.db") 
## Constructs a data framecontaining the gene IDs, gene symbols
##  and descriptions for all probe sets on the chip.
Annot <- data.frame(ACCNUM=sapply(contents(ath1121501ACCNUM),
					paste, collapse=", "), 
					SYMBOL=sapply(contents(ath1121501SYMBOL),
					paste, collapse=", "), 
                    DESC=sapply(contents(ath1121501GENENAME),
                    paste, collapse=", ")) 

all<-merge(Annot, my_frame, by.x=0, by.y=0, all=T) 
# Merges everything with above expression data.
write.table(all, file="my_annot_file.xls", sep="\t", col.names = NA)
# Exports data to text file that can be imported into Excel.

## Visualization and quality control
d <- cor(2^exprs(eset_rma), method="pearson")
plot(hclust(dist(1-d))) 
# Generates a correlation matrix for all-against-all chip comparisons.
library(affyQCReport)
QCReport(mydata, file="ExampleQC.pdf") 
# Generates a comprehensive QC report for the AffyBatch object 'mydata' in PDF format. See affyQCReport for details.
eset <- eset_rma

##  Simple Affy (command summary)
## Save this covdesc.txt to your R working directory.
## Generate expression data with RMA and MAS5.
## Filter each of the three data sets with the following parameters: 2-fold changes, present in all 4 chips and p-score less than 0.001.
## Write the results into separate files.
## Create scatter plots for the filtered data sets and save them to external image files.
## Identify the overlap of the significant changes between the RMA and MAS5 data.
## Perform simpleaffy QC checks: scaling factor, percent present calls, etc.

## (C.1) Identify differtially expressed genes (DEGs) with the limma package
## Commands from the affy Limma Section in manual. 
## cDNA microarray users can save and extract the SWIRL cDNA microarray sample 
## data. For a quick demonstration of the analysis of this data set, one can 
## copy&paste or source the following command-line summary into the R terminal: my_swirl_commands.txt.

#targets file:
# http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/affy_targets.txt
library(limma) 
targets <- readTargets("affy_targets.txt")
# Makes a data frame from file 'affy_targets.txt' 
data <- ReadAffy(filenames=targets$FileName) 
# Reads CEL files (specified in 'targets') into AffyBatch object.
eset <- rma(data)
# Normalizes data with 'rma' function and
#  assigns them to ExpressionSet object
# exprs(eset) <- log2(exprs(eset))
# # If eset contains absolute intensity values like MAS5 results,
#  then they should be transformed to log2 (or loge) values
# RMA/GCRMA generate log2 values and MAS5 produces absolute values.
pData(eset) 
# Lists the analyzed file names.
write.exprs(eset, file="affy_all.xls") 
# Exports all affy expression values to tab delimited text file. 
# The MAS 5.0 P/M/A calls can be retrieved with simpleaffy package 
# or with the affy package like:
#  'eset <- mas5calls(data); write.exprs(eset, file="my_PMA.txt")'
design <- model.matrix(~ -1+factor(c(1,1,2,2,3,3))) 
# Creates appropriate design matrix 
# design matrix can be created in any spreadsht program + imported to R
colnames(design) <- c("group1", "group2", "group3")
fit <- lmFit(eset, design)
# Fits a linear model for each gene 
# based on the given series of arrays.
contrast.matrix<-makeContrasts(group2-group1, group3-group2,
							   group3-group1, levels=design) 
# contrast matrix to perform all pairwise comparisons. 
# a contrast matrix can be created in any spreadsht program 
# and then imported into R.
#  For complex experiments one can also use this function to
#  compute a contrast matrix with all possible pairwise comparisons
#  
fit2 <- contrasts.fit(fit, contrast.matrix) 
# Computes estimated coefficients and standard errors 
# for a given set of contrasts.
fit2 <- eBayes(fit2)
# Computes moderated t-statistics and log-odds of diff-exp
#by empirical Bayes shrinkage of standard errors to common value
topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=10) 
# Generates list of top 10 ('number=10') diff-exp genes
#  sorted by B-values ('sort.by=B') 
#  for each of the three comparison groups ('coef=1') 
# The summary table contains the following information: 
# logFC is the log2-fold change, 
# the AveExpr is the average expression value accross arrays and channels,
# the moderated t-statistic (t) is the logFC to its standard error,
#  the P.Value is the associated p-value, 
#  the adj.P.Value is the p-value adjusted for multiple testing 
#  and the B-value (B) is the log-odds that a gene is diff-expr
#    'limmaUsersGuide()'
write.table(topTable(fit2, coef=1, adjust="fdr", sort.by="B",
					 number=50000), file="limma_complete.xls", 
			         row.names=F, sep="\t") 
# Exports complete limma statistics table for first comparison 
# group ('coef=1') to tab delimited text file.
results <- decideTests(fit2, p.value=0.05); vennDiagram(results) 
# Creates venn diagram of changed genes with p-value less than 0.05.
x <- topTable(fit2, coef=1, adjust="fdr", sort.by="P", number=50000)
y <- x[x$adj.P.Val < 0.05,];
y
print("Number of genes in this list:")
length(y$ID) # Filters candidates with P-values < 0.05 in ('coef=1')
# and provides the number of candidates for each list.
# These numbers should be identical with the sum of the values in each circle of the above venn diagram.
x <- topTable(fit2, coef=1, adjust="fdr", sort.by="P", number=50000); y <- x[x$adj.P.Val < 0.01 & (x$logFC > 1 | x$logFC < -1) & x$AveExpr > 10,]; y; print("Number of genes in this list:"); length(y$ID) # Same as above but with complex filter: P-value < 0.01 AND at least 2-fold change AND expression value A > 10.
results <- decideTests(fit2, p.value=0.000005); heatDiagram(results, fit2$coef, primary=1) # This function plots heat diagram gene expression profiles for genes which are significantly differentially expressed in the primary condition (this is not a cluster analysis heat map). Genes are sorted by differential expression under the primary condition. The argument 'primary=1' selects the first contrast column in the 'results' matrix as primary condition. The plotted genes can be extracted like this 'results[results[,1]==1,]'. More information on this function can be found in the limma manual.


## (D.1) GO Term enrichment analysis with GOstats
library(ath1121501.db); library(GOstats)
affySample <- c("266592_at", "266703_at", "266199_at", "246949_at", "267370_at", "267115_s_at", "266489_at", "259845_at", "266295_at", "262632_at")
geneSample <- as.vector(unlist(mget(affySample, ath1121501ACCNUM, ifnotfound=NA)))
library(ath1121501cdf)
affyUniverse <- ls(ath1121501cdf)
geneUniverse <- as.vector(unlist(mget(affyUniverse, ath1121501ACCNUM, ifnotfound=NA)))
params <- new("GOHyperGParams", geneIds = geneSample, universeGeneIds = geneUniverse, annotation="ath1121501", 
               ontology = "MF", pvalueCutoff = 0.5, conditional = FALSE, testDirection = "over")
hgOver <- hyperGTest(params); summary(hgOver)
htmlReport(hgOver, file = "MyhyperGresult.html") 

## (D.2) GO batch analysis using the GOHyperGAll script
## Download and unzip the annotation objects for Arabidopsis (http://faculty.ucr.edu/%7Etgirke/Documents/R_BioCond/Samples/ArabSampleGOHyperGAll.zip) 
## into your working directory. Then continue with the following commands:
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/GOHyperGAll.txt") # Imports the GOHyperGAll functions.
loadData(); load(file="MF_node_affy_list"); load(file="BP_node_affy_list"); load(file="CC_node_affy_list") # Loads the downloaded annotation objects for Arabidopsis.
GOHyperGAll_result <- GOHyperGAll(gocat="MF", sample=unique(as.vector(GO_MF_DF[1:40,2])), Nannot=2); GOHyperGAll_result[1:10,-8] # Performs the enrichment test for provided set of gene identifiers.
CL_DF <- data.frame(geneID=GO_MF_DF[1:400,2], ClusterID=sort(rep(1:4,100)), ClusterSize=rep(100,400)) # Create sample data set for batch processing.
BatchResult <- GOCluster_Report(CL_DF=CL_DF, method="all", id_type="gene", CLSZ=10, cutoff=0.001, gocats=c("MF", "BP", "CC"), recordSpecGO=c("GO:0003674", "GO:0008150", "GO:0005575")); BatchResult[1:4,-10] # Performs all three GO analyses for many sample set at once. When the method argument is set to "slim" then the goSlim method is used.
write.table(BatchResult, file="GO_Batch.xls", sep="\t", col.names = NA) # Exports batch result to Excel file.

## (E.1) Clustering of differentially expressed genes
my_fct <- function(x) hclust(x, method="complete") # Creates function to perform hierarchical clustering (complete linkage).
heatmap(as.matrix(2^exprs(eset)[1:40,]), col = cm.colors(256), hclustfun=my_fct) # Plots heatmap with dendrograms.
## Replace in last step 'exprs(eset)[1:40,]' by matrix of differentially expressed genes from limma analysis.

## (F.1) Clean up this R script and execute it with the source() function
source("array.R") # Executes all of the above commands and generates the corresponding output files in the current working directory.



