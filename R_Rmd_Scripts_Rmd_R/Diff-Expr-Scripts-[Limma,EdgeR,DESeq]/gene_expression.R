# expressionAnalysisSampleCode.R
#
# Purpose: Sample processing of gene expression studies with RNA seq and
#            microarray platforms
# Version: 1.0
# Date:    2018 01
# Author:  Boris Steipe <boris.steipe@utoronto.ca>
#
# Input:
# Output:
# Dependencies:
#
# Version history:
#
# ToDo:
# Notes:  R Code adapted from RPR-GEO2R  ABC learning unit.
#
#         Quantile normalization needs to cite:
#           Bolstad, B. M., Irizarry R. A., Astrand, M, and Speed, T. P. (2003)
#           A Comparison of Normalization Methods for High Density
#           Oligonucleotide Array Data Based on Bias and Variance.
#           Bioinformatics 19(2) ,pp 185-193.
#
#         Should we calculate DE instead and normalize that? There is otherwise
#         a danger of ML fitting the noise more than the otherwise sparse
#         signal. We could use (test - control) * pDE and set all insignificant
#         changes to 0.?
#
# ==============================================================================


# ====  PARAMETERS  ============================================================

# The microarray demo set is GSE35330
# The RNAseq demo set is GSE63807


# ====  PACKAGES  ==============================================================

if (! require(Biobase, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("Biobase")
  library(Biobase)
}

if (! require(GEOquery, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("GEOquery")
  library(GEOquery)
}

# for quantile normalization ...
if (! require(preprocessCore, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("preprocessCore")
  library(preprocessCore)
}



# ====  FUNCTIONS  =============================================================
# ...

# ====  PROCESS  ===============================================================

# ==============================================================================
#
# ==== Microarray data =========
#
# ==============================================================================


# GSE35330 <- getGEO("GSE35330", GSEMatrix =TRUE, AnnotGPL=TRUE)
# Loads the data set as well as the platform information (GPL570 in this case).
#
# GSE35330 <- GSE35330[[1]]
# save(GSE35330, file="GSE35330.RData")
load("GSE35330.RData")

# ==== What data do we have?
nrow(exprs(GSE35330))      # 54675
ncol(exprs(GSE35330))      # 24 slides: 6 conditions with 4 replicates each
                           # I want to use control and test for the 16h
                           # samples.
colnames(exprs(GSE35330))  # sample names ...

# Assess distributions
#    define a color scheme: controls: greens, test: purples
c1 <- colorRampPalette(c("#C27E7E", "#816EBA"))(3)   # test
c2 <- colorRampPalette(c("#758AC9", "#82C9B6"))(3)   # ctrl

# a color vector for the 24 samples ...
myArrayCols <- c(rep(c2[1], 4),   # ctrl  4h
                 rep(c2[2], 4),   # ctrl 16h  ***
                 rep(c2[3], 4),   # ctrl 64h
                 rep(c1[1], 4),   # test  4h
                 rep(c1[2], 4),   # test 16h  ***
                 rep(c1[3], 4))   # test 64h

# samples are stored as test1 - ctrl1 - test2 - ctrl2 - test3 - ctrl3.
# I want test1 - test2 - test3 - ctrl1 - ctrl2 - ctrl3 instead. Reorder
# the columns
iReorder <- c(5:8, 13:16, 21:24, 1:4, 9:12, 17:20)

boxplot(log(exprs(GSE35330)[ , iReorder]),
        boxwex = 0.6,
        notch = TRUE,
        main = "GSE35330",
        outline = FALSE,
        col = myArrayCols)

# ==== extract columnd in new order

myEx <- exprs(GSE35330)[ , iReorder]
colnames(myEx) <- c("c04a", "c04b", "c04c", "c04d",
                    "c16a", "c16b", "c16c", "c16d",
                    "c64a", "c64b", "c64c", "c64d",
                    "t04a", "t04b", "t04c", "t04d",
                    "t16a", "t16b", "t16c", "t16d",
                    "t64a", "t64b", "t64c", "t64d")

boxplot(log(myEx),
        boxwex = 0.6,
        notch = TRUE,
        main = "GSE35330",
        outline = FALSE,
        col = myArrayCols)

# How to normalize? One way is to quantile normalize all replicate sets
# individually, but keep the trend differences between them.
myEx[ ,   1:4] <- normalize.quantiles(myEx[ ,   1:4], copy = FALSE)
myEx[ ,   5:8] <- normalize.quantiles(myEx[ ,   5:8], copy = FALSE)
myEx[ ,  9:12] <- normalize.quantiles(myEx[ ,  9:12], copy = FALSE)
myEx[ , 13:16] <- normalize.quantiles(myEx[ , 13:16], copy = FALSE)
myEx[ , 17:20] <- normalize.quantiles(myEx[ , 17:20], copy = FALSE)
myEx[ , 21:24] <- normalize.quantiles(myEx[ , 21:24], copy = FALSE)

boxplot(log(myEx),
        boxwex = 0.6,
        notch = TRUE,
        main = "GSE35330",
        outline = FALSE,
        col = myArrayCols)

# ==== Prepare annotation data
str(GSE35330@featureData@data)   # Have annotations been properly loaded ?
myAnnot <- GSE35330@featureData@data[ , c("ID", "Gene symbol")]
str(myAnnot)                                        # confirm
colnames(myAnnot) <- c("probeIDs", "symbols")       # rename
myAnnot$probeIDs <- as.character(myAnnot$probeIDs)  # convert to character
myAnnot$symbols  <- as.character(myAnnot$symbols)

any(is.na(myAnnot$symbols))       # FALSE
sum(myAnnot$symbols == "")        # 9557 probes are not annotated with a
                                  # HUGO symbol.We will ignore these.

sum(grepl("/", myAnnot$symbols))  # 2217 probes are annotated with more than
                                  # one symbol. We will ignore these too.

sum(duplicated(myAnnot$symbols[myAnnot$symbols != ""])) # 22929 duplicated
                                                        # symbols (not counting
                                                        # the un-annotated
                                                        # rows).

# How many unique symbols do these contain?
x <- unique(myAnnot$symbols[duplicated(myAnnot$symbols)]) # 11404 ...
# i.e. more than half of symbols are measured more than once on the array.

# Let's look at some examples:
head(which(duplicated(myAnnot$symbol[myAnnot$symbol != ""])))

# first example: seven measurements for MAPK1
myAnnot[myAnnot$symbol == myAnnot$symbol[19], ]

sel <- which(myAnnot$symbol == myAnnot$symbol[19])  # MAPK1
sel <- which(myAnnot$symbol == myAnnot$symbol[23])  # PRR22
sel <- which(myAnnot$symbol == myAnnot$symbol[25])  # PXK
sel <- which(myAnnot$symbol == myAnnot$symbol[29])  # SLC46A1
sel <- which(myAnnot$symbol == myAnnot$symbol[36])  # CILP2
sel <- which(myAnnot$symbol == myAnnot$symbol[44])  # TMEM106A

# plot expression values
plot(1:24, 1:24,
     ylim = c(min(log10(myEx[sel, ])) * 0.95, max(log10(myEx[sel, ])) * 1.05),
     xlab = ("# sample"),
     ylab = ("log10(expression)"),
     type = "n")
for (idx in sel) {
  points(1:24,
         log10(myEx[idx, ]),
         pch = 16,
         col = myArrayCols,
         type = "b")
}

# For MAPK1 we have expression data sets that span more than one order of
# magnitude. Remember: all lines are different probes OF THE SAME mRNA POOL,
# and all identical colours are replicates that should have the same value.

# simple average of all profiles
points(1:24,
     apply(log10(myEx[sel, ]), 2, mean),
     col = "seagreen",
     lwd = 2,
     type = "l")

# But are all these experiments equally good? We could conjecture that
# experiments with lower absolute values have greater noise. If this were the
# case the relative standard deviation should be smaller if the mean is
# larger.

# Look at relative standard deviations (coefficient of variation) as a function
# of means.
# Randomly choose 5000 sets
set.seed(11235)
N <- 5000
SDM <- matrix(numeric(2 * N), nrow = N) # matrix to hold sd() and mean() values

Nrow <- nrow(myEx)
col1 <- c(5, 17)  # first column of the replicate sets we are interested in
for (i in 1:N) {
  r <- sample(1:Nrow, 1) # random row
  c <- sample(col1, 1)   # randomly either of the replicate sets
  x <- myEx[r , c:(c+3)] # fetch the four expression values
  SDM[i, 2] <- mean(x)             # calculate mean
  SDM[i, 1] <- sd(x) / SDM[i, 2]   # calculate relative SD
}

# plot the distribution of points
plot(SDM[ , 2],
     SDM[ , 1],
     pch = 16,
     col = "#3366FF22",
     log = "xy",
     xlab = "mean",
     ylab = "relative SD")
abline(h =  1, lwd = 1.5, col = "#33FF6644") # at this level rSD == mean
abline(v = 50, lwd = 1.5, col = "#33FF6644") # below this, the relationship
                                             # is no longer linear

# Not only do we see that the rSDs fall off with intensity, and that this
# is linear on the log/log plot (i.e. a polynomial relationship) we also see
# that one order of magnitude of rSD is gained for two orders of magnitude
# of intensity increase. This is exactly what we would expect from random
# population sampling, where the relative error drops with the square root
# of the sample size. There are no other systematic effects (above a threshold
# of ~ 50).

# Therefore, in order to estimate the expression values for each sample, we
# should take a weighted column-wise average of those samples that have been
# measured more than once. Then we take the averages within the replicates.

# But BEFORE that - while we still have individual measurements, should we
# remove outliers? And if we do, how should we impute the actual values?


# ==================================================

# === More considerations of the symbol annotations. How do the existing
# symbols compare to our traget vector of gene symbols?

load("./inst/extdata/HUGOsymbols.RData")

# How many target symbols do we have covered?
sum(HUGOsymbols %in% unique(myAnnot$symbol)) # 16861: 83 %

# make a vector of missing Symbols
missingSymbols <- HUGOsymbols[!(HUGOsymbols %in% unique(myAnnot$symbol))]

# What annotation symbols do we have that are NOT in our symbol list?
x <- unique(myAnnot$symbol)
extraSymbols <- x[!(x %in% HUGOsymbols)]
head(extraSymbols, 50)

# [45] "ATRIP///TREX1"  - some of these are probes that identify
# more than one gene.
"ATRIP" %in% HUGOsymbols # TRUE
"TREX1" %in% HUGOsymbols # TRUE

# Since we can't attribute these uniquely to a gene, we should ignore them.

# What are the others? Look them up at https://www.genenames.org/
# [5] "AFG3L1P"        - pseudogene
# [9] "C10orf25"
# [11] "FAM167A-AS1"   - antisense
# [13] "LINC00161"
# [14] "DSCR10"        - "Down syndrome critical region"
# [15] "LACE1"         - NEW NAME: AFG1L
# [16] "LINC01547"
# [17] ""
# [18] "NUDT9P1"       - pseudogene
# [19] "ALMS1P1"       - pseudogene
# [20] "ABCC13"        - pseudogene
# [21] "LINC00308"
# [24] "RFWD2"         - NEW NAME: COP1
# [29] "JMJD1C-AS1"
# [30] "TUBA3FP"       - pseudogene
# [44] "ZCCHC5"        - NEW NAME: RTL3
# [50] "ZNF645"        - NEW NAME: CBLL2

# Estimate 8% (~ 400) recoverable symbols ( 10% of missing, raising coverage
# to 85%).

# I have prepared a symbolMap dataframe so that we can use match() to replace
# many symbols that are not in the HUGO list with some that are ...

# But how does match() work? Sample code:
tmp <- data.frame(sym = c("A", "oldB", "C", "oldD", "unkE", "F"),
                  val = 16:11,
                  stringsAsFactors = FALSE)

aliasMap <- data.frame(old = c("oldE", "oldD", "oldBi", "oldA", "oldX", "oldB"),
                       new = c(   "E",    "D",     "B",    "A",    "X",    "B"),
                       stringsAsFactors = FALSE)

(m <- match(tmp$sym, aliases$old))  # this vector  describes the matches
(sel <- ! is.na(m))                 # this selects the elements we can replace
tmp$sym[sel] <- aliases$new[m[sel]] # this replaces old with new
tmp

# ==============================================================================
#
# ==== RNAseq data =========
#
# ==============================================================================

# Work with RNAseq seems much easier by comparison: we already get unique gene
# symbols.

# download the normalized, annotated dataset from GEO
#
# URL <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63807/suppl/"
# fName <- "GSE63807_NormalizedData_withAnnotationsKE.txt.gz"
# library(readr)
# GSE63807 <- read_tsv(paste0(URL, fName))
# str(GSE63807)
# GSE63807 <- data.frame(symbols = GSE63807$Symbol,
#                        ctrl1   = GSE63807$Neg_co_1,
#                        ctrl2   = GSE63807$Neg_co_2,
#                        ctrl3   = GSE63807$Neg_co_3,
#                        test1   = GSE63807$Stab1siRNA_1,
#                        test2   = GSE63807$Stab1siRNA_2,
#                        test3   = GSE63807$Stab1siRNA_2,
#                        stringsAsFactors = FALSE)
# str(GSE63807)
# save(GSE63807, file="GSE63807.RData")
load("GSE63807.RData")

# what does this data look like
str(GSE63807)

nrow(GSE63807) # 22479 rows
any(is.na(GSE63807$symbols))       # FALSE
any(GSE63807$symbols == "")        # FALSE
any(grepl("/", GSE63807$symbols))  # FALSE
any(duplicated(GSE63807$symbols))  # FALSE

sum(HUGOsymbols %in% GSE63807$symbols) # 17593: 86 %

missingRNAseqSymbols <- HUGOsymbols[!(HUGOsymbols %in% GSE63807$symbols)]

# how many of these can be recovered, because they are aliases?

extraRNAseqSymbols <- GSE63807$symbols[!(GSE63807$symbols %in% HUGOsymbols)]
# what are these?
cat(sprintf("# [%d] %s\n", 1:20, head(extraRNAseqSymbols, 20)))
# [1] A1BG-AS1       -  antisense
# [2] A2MP1          -  pseudogene
# [3] AA06           -  ???
# [4] AAA1           -  synonym for NPSR1-AS1
# [5] AACSP1
# [6] AATK-AS1
# [7] ABCA11P
# [8] ABCA17P
# [9] ABCC13         - pseudogene
# [10] ABCC6P1
# [11] ABCC6P2
# [12] ABHD14A-ACY1  - readthrough (NMD candidate)
# [13] ABP1          - NEW SYMBOL: AOC1
# [14] ACN9          - NEW SYMBOL: SDHAF3
# [15] ACPL2         - NEW SYMBOL: PXYLP1
# [16] ACPT          - NEW SYMBOL: ACP4
# [17] ACRC          - NEW SYMBOL: GCNA
# [18] ACTR3BP2
# [19] ACTR3BP5
# [20] ADAM21P1

load("./inst/extdata/synMap.RData")
str(synMap)
# mappable:
x <- extraRNAseqSymbols[extraRNAseqSymbols %in% synMap$synonyms] # 772

# recovering these 772 symbols will raise coverage
# to 17593 + 772 = 18365 : i.e. 90%

# value distributions:

boxplot(GSE63807[ , 2:7],
        boxwex = 0.6,
        notch = TRUE,
        main = "GSE63807",
        outline = FALSE,
        col = c(rep(myArrayCols[5], 3), rep(myArrayCols[17], 3)))

# etc ...


# [END]
