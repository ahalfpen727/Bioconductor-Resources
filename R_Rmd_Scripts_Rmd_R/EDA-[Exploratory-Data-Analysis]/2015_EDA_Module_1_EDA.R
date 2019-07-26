# ================================================== #
# Canadian Bioinformatics Workshops Series           #
# Toronto, May 21 and 22 2015                        #
# Exploratory Analysis of Biological Data using R    #
#                                                    #
# Faculty: Boris Steipe <boris.steipe@utoronto.ca>   #
#                                                    #
# Some prior contributions by:                          #
#   Raphael Gottardo, FHCRC                          #
#   Sohrab Shah, UBC                                 #
#                                                    #
#                                                    #
# Module 1: EDA (Exploratory Data Analysis)          #
#                                                    #
# ================================================== #


# ==================================================
# Setup a project and work-environment
# ==================================================

# Task:
# 1 - Create a project directory on your computer,
#     call it "CBW".
# 2 - Open R studio, and use the menu to create
#     a new project in this directory.
# 3 - Open a new R script file.
# 4 - Open this file on the workshop wiki. Copy its
#     contents and paste it into your script file.
# 5 - Save the script under a name of 01-EDA.R
#     The code and comments you write during this
#     workshop will all be kept together in this file.
# 6 - Also load PlottingReference.R
# 7 - Open a third script file, call this "tmp.R";
#     use this script as a scratchpad.
# 8 - Set the current "working directory" to your CBW
#     workshop directory.
# 9 - Confirm:
getwd()   # Confirm the correct directory
dir()     # Confirm that the right files are present.

# Optional task:
# If you are waiting for others to finish, try the
# following:
# 1 - Install/load the "RUnit" package.
# 2 - Explore the vignettes it contains.
# 3 - Explore the functions it contains.


# ==================================================
# Utilities
# ==================================================

# It's useful to keep some utility functions loaded
# for your project.

# In yesterday's course we have used the following 
# function a lot:

typeInfo <- function(x) {
    print(x, digits=22)  
    cat("str:    ")                
    str(x)  
    cat("mode:   ", mode(x), "\n")
    cat("typeof: ", typeof(x), "\n")
    cat("class:  ", class(x), "\n")
    # if there are attributes, print them too
    if (! is.null(attributes(x))) {
        cat("attributes:\n")
        print(attributes(x))
    }
}

# Task:
# 1 - If you haven't done so already, put typeInfo() into
#     a file "utilities.R" in your working directory.
# 2 - Then:
source("utilities.R")

# ==================================================
# Load Data
# ==================================================

# Task:
# In yesterday's workshop we have worked with a
# a supplementary datafile from a 2014 publication
# on single-cell RNAseq for the automatic definition
# of tissue types. The data consists of gene-names,
# Expression values of clustered cells in the
# presence and absence of LPS stimulation, and
# cluster assignment labels. The corresponding
# paper is on th Wiki. ()

# Task:
# 1 - If you don't have the file yet, download the
#     .csv version of the file from the Workshop
#     Wiki.
# 2 - Read it into R with the following command:

sup3 <- read.csv("table_S3.csv",
                  skip = 6,
                  header = FALSE,
                  colClasses = c("character", rep("numeric", 15)),
                  col.names = c("genes",
                                  "B.ctrl", "B.LPS",
                                  "MF.ctrl", "MF.LPS",
                                  "NK.ctrl", "NK.LPS",
                                  "Mo.ctrl", "Mo.LPS",
                                  "pDC.ctrl", "pDC.LPS",
                                  "DC1.ctrl", "DC1.LPS",
                                  "DC2.ctrl", "DC2.LPS",
                                  "cluster"),
                  stringsAsFactors = FALSE)
str(sup3)
head(sup3)


# ==================================================
# Descriptive statistics
# ==================================================

set.seed(100)
x <- rnorm(100, mean=0, sd=1)
mean(x)
median(x)
IQR(x)
var(x)
sd(x)
summary(x)

# Task:
# 1 - What would be interesting characterizations
#     of sup3? Explore them.

# =============================================
# Quantiles: what is the threshold that has a 
# given fraction of values above/below it?

x <- seq(-4, 4, 0.1)
f <- dnorm(x, mean=0, sd=1)
q90 <- qnorm(0.90, mean = 0, sd = 1)
plot(x, f, xlab="x", ylab="density", type="l", lwd=5)
abline(v=q90, col=2, lwd=5)

# =============================================
# empirical quantiles

set.seed(100)
x <- rnorm(100, mean=0, sd=1)
quantile(x)
quantile(x, probs=c(.1, .2, .9))

# =============================================

set.seed(100)
x <- rnorm(1000, mean=0, sd=1)
boxplot(x)

# Careful - a boxplot per se can obscure
# important structure in your data. Consider
# for example this bimodal distribution
x <- rnorm(100, mean=-2, sd=1)
x <- c(x, rnorm(100, mean=2, sd=1))
hist(x)
boxplot(x)

# Violin plot
if (!require(ggplot2, quietly=TRUE)) {
  install.packages("ggplot2")
  library(ggplot2)
}
X <- as.data.frame(x)

p <- ggplot(X, aes(1,x))
p + geom_violin()
# See ggplot2 introductory tutorial at 
#     http://www.r-bloggers.com/basic-introduction-to-ggplot2/
# and violin plot documentation at
#     http://docs.ggplot2.org/current/geom_violin.html


# Plotting more than one column with a boxplot
# places the plots side by side.
colnames(sup3)
boxplot(sup3[ ,c("MF.ctrl", "MF.LPS")])
boxplot(sup3[ ,2:14])


# Task:
# 1 - What would be interesting quantiles
#     and boxplots in sup3? Explore this.
#     Interpret the result.



# ==================================================
# Explore plot types
# (Section 1 of PlottingReference.R)
# ==================================================

# ==================================================
# Probability distributions
# ==================================================

# The binomial distribution
?dbinom
x <- 0:1
f <- dbinom(x, size=1, 0.1)
plot(x, f, xlab="x", ylab="density", type="h", lwd=5)

set.seed(100)
x <- rbinom(100, size=1, prob=.1)
hist(x)


# ==================================================
# The Normal distribution
# ==================================================
?dnorm
x <- seq(-4, 4, 0.1)
f <- dnorm(x, mean=0, sd=1)
plot(x, f, xlab="x", ylab="density", lwd=5, type="l")


# ==================================================
# Explore lines (Section 3 of PlottingReference.R)
# ==================================================

x <- seq(-4, 4, 0.1)
f <- dnorm(x, mean=0, sd=1)
plot(x, f, xlab="x", ylab="density", lwd=5, lty=3, type="l")

# =============================================
# Overlay a histogram with a line
set.seed(100)
x <- rnorm(100, mean=0, sd=1)
hist(x)
lines(seq(-3,3,0.1),50*dnorm(seq(-3,3,0.1)), col="red")


# ==== QQ PLOTS ====================================

set.seed(100)
x <- rnorm(100, mean=0, sd=1)
qqnorm(x)
qqline(x, col=2)

# =============================================
# Plotting lines and legends
# Example: compare the normal distribution with
# the t-distribution
x <- seq(-4, 4, 0.1)
f1 <- dnorm(x, mean=0, sd=1)
f2 <- dt(x, df=2)
plot(x, f1, xlab="x", ylab="density", lwd=5, type="l")
lines(x, f2, lwd=5, col=2)
legend(-4, .4, c("Normal", "t2"), col=1:2, lwd=5)

# =============================================
# use qqplot to compare normally distributed samples with
# t-distributed samples

set.seed(100)
t <- rt(100, df=2)
qqnorm(t)
qqline(t, col=2)


# =============================================
# QQ- plot: sample against sample
set.seed(100)
x <- rnorm(100, mean=0, sd=1)
t <- rt(100, df=2)
qqplot(x, t)

# Task:
# 1 - What columns of sup3 could be compared with
#     a qqplot? Explore this. Interpret the result.

# ==================================================
# Flow cytometry data
# ==================================================

# Task:
# 1 - Download the data file GvHD.txt and save it to
#     your project directory.


# GvHD flow cytometry data
gvhd <- read.table("GvHD.txt", header=TRUE)
# Only extract the CD3 positive cells
gvhdCD3p <- as.data.frame(gvhd[gvhd[, 5]>280, 3:6])
plot(gvhdCD3p[, 1:2])

# ==================================================
# Explore scatter plots
# Topics:
# Section 6 - Plotting symbols and characters
#         2 - Colors
#         4 - Coordinates
# ... of PlottingReference.R)
# ==================================================
# Scatter plots are extremely useful, but learning
# a bit about R's plotting capabilities will help
# tremendously to create INFORMATIVE plots.


# Some special packages
# The overlap in the GvHD data can obscure
# data relationships. Here are some alternatives
# for dense plots:

# load "hexbin" package from CRAN
if (!require(hexbin, quietly=TRUE)) {
  install.packages("hexbin")
  library(hexbin)
}

# variant one: hexbin
hb <- hexbin(gvhdCD3p[, 1:2], xbins = 20)
plot(hb, colramp = colorRampPalette(c("#FFFFDD",
                                      "#77EE99",
                                      "#3377AA",
                                      "#0033BB")))

# load "prada" package from BioConductor
if (!require(prada, quietly=TRUE)) {
    source("http://www.bioconductor.org/biocLite.R")
    biocLite("prada")
}

# variant two: smoothScatter
smoothScatter(gvhdCD3p[, 1:2], nrpoints=50, pch=20, cex=0.5, col="#6633BB55")

# variant three: colors vary by density
plot (gvhdCD3p[, c(1,2)], col=densCols(gvhdCD3p[, 1], gvhdCD3p[, 2]), pch=16)




# Task:
# 1 - Create a scatterplot from differences in
#     LPS activation between a pair of cell types in sup3.
# 2 - Add a line that shows the situation if all
#     activations were equal.
# 3 - Redo the plot with density Coloring

plot(sup3[ ,"Mo.ctrl"] - sup3[ ,"Mo.LPS"], 
     sup3[ ,"MF.ctrl"] - sup3[ ,"MF.LPS"])
abline(0,1,col = "red")

# =============================================
# Trellis plots: all against all

plot(gvhdCD3p, pch=".")

# =============================================

boxplot(gvhdCD3p)

# =============================================

oPar <- par(mfrow=c(2, 2))
hist(gvhdCD3p[, 1], 50, main=names(gvhdCD3p)[1], xlab="fluorescent intensity", ylim=c(0, 120))
hist(gvhdCD3p[, 2], 50, main=names(gvhdCD3p)[2], xlab="fluorescent intensity", ylim=c(0, 120))
hist(gvhdCD3p[, 3], 50, main=names(gvhdCD3p)[3], xlab="fluorescent intensity", ylim=c(0, 120))
hist(gvhdCD3p[, 4], 50, main=names(gvhdCD3p)[4], xlab="fluorescent intensity", ylim=c(0, 120))
par(oPar)


# =============================================

# ==================================================
# Explore 3D plots
# 
# Section 8 - 8 - Plots of X-Y-Z coordinates
# ... of PlottingReference.R)
# ==================================================


# [End]
