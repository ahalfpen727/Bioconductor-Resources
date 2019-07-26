# ==============================================================================
# Exploratory Analysis of Biological Data using R, 2015
# Integrated Assignment Part 1 Answers
#
# @Authors:  David JH Shih <djh.shih@gmail.com>
#            Catalina Anghel <catalina.anghel@oicr.on.ca>
# @License:  GNU General Public License v3 
# @Created:  2015-05-11
#
# @Input:    Data from the Cancer Cell Line Encyclopedia (CCLE)
#            Downloaded from: http://www.broadinstitute.org/ccle
#            Associated publication:
#
#            Barretina et al. The Cancer Cell Line Encylopedia enables
#            predictive modelling of anticancer drug sensitivity.
#            Nature 483, 603-607 (2012).
#
# @Output:   Statistics, Plots


# ==============================================================================
# SECTION 1  Install the required R package
#

# You should have built and installed the CLLE package.


# ==============================================================================
# SECTION 2  Import the data into R
#

# 2a  Load the "CCLE" R library by using the library() function.
library(CCLE);

# 2b  Import the data into the workspace
data(ccleCgc);


# ==============================================================================
# SECTION 3  Examine the environment and objects
# 

# 3a  Get the list of objects in the workspace using ls().
ls();

# 3b  What is the difference between ls() and list.files()?
list.files();
# Answer:  ls() lists objects in the R environment and list.files() lists the
#          files stored on the disk.

# For now, we will only work with `expr`, `cn`, and `pheno`.

# 3c  What are the classes of `expr`, `cn`, and `pheno`?
#     e.g. class(expr) outputs the class of `expr`
class(expr);
class(cn);
class(pheno);

# 3d  How many columns and rows does `expr` have? (Hint: nrow(), ncol())
nrow(expr);
ncol(expr);
# Alternatively
dim(expr);

# 3e  How many columns and rows do the data frames `cn`, `expr`, `pheno` have?
dim(cn);
dim(expr);
dim(pheno);

# 3f  What are the names of the rows of `expr`? The names of the columns?
rownames(expr);
colnames(expr);

# 3g  What are the names of the rows of `pheno`? The names of the columns?
rownames(pheno);
colnames(pheno);


# ==============================================================================
# SECTION 4  Convert and rearrange the data
# 

# 4a  Restructure the data into an easy-to-use format.
# Already done.

# 4b  Examine the description of the data by typing `?ccleCgc`.
#     Type `q` to get back to the command line, after reading the description.
#     For all the data frames (e.g. `cn`), what do the rows represent?
# Answer:  Each row represents a cancer cell line.

# We should make sure that the data frames describe the same cell lines
# and that the data are arranged in the same order.
# We can do so by ensuring that the row names are the same.
# The code below does this check for for the data frames `cn` and `expr`.
identical(rownames(cn), rownames(expr));

# 4c  Now repeat this comparison between `cn` and `pheno`.
identical(rownames(cn), rownames(pheno));


# ==============================================================================
# SECTION 5  Examine the data distributions
# 

# 5a  Examine the usage instructions for hist(), using ? hist.
# Note:  hist(x) does not work if x is a data.frame.
#        x must be a vector. Matrices are automatically converted to vectors.
#        as.matrix(x) converts x into a matrix.
#     Plot a histogram of `expr`.
hist(as.matrix(expr));

# 5b  Plot a prettier histogram of `expr`.
#     Use 50 breaks and set the colour of the bars to "skyblue".
#     Additionally, set the pararmeter freq to FALSE.
# Hint:  In the usage instructions of hist(), what do the parameters col, 
#        breaks, and freq do?
# Note:  Don't forget to convert `expr` into a matrix.
hist(as.matrix(expr), breaks=50, col="skyblue", freq=FALSE);

# If hist() have been called with freq=FALSE, you can superimpose a 
# density curve by the following command:
lines(density(as.matrix(expr)), col="red");

# 5c  Plot a histogram of `cn` with 100 breaks and orange bars.
#     Superimpose a density plot.
hist(as.matrix(cn), breaks=100, freq=FALSE, col="orange");
lines(density(as.matrix(cn)), col="red");

# 5d  Create a quantile-quantile plot of `cn` against the normal 
#     distribution. Draw a reference line through the data.
# Hint:  Use qqnorm() to compare against a normal distribution.
qqnorm(as.matrix(cn), pch='.');
qqline(as.matrix(cn), col="red")

# You can also superimpose a normal distribution on top of a histogram:
z <- as.numeric(as.matrix(cn));
hist(z, breaks=100, freq=FALSE, col="orange");
curve(dnorm(x, mean=mean(z), sd=sd(z)), col="blue", add=TRUE);


# ==============================================================================
# SECTION 6  Explore the gene expression data
# 

# 6a  Optional: Check if the tumor suppressor gene TP53 is in the dataset. 
# Hint:  There are two ways: one using `%in%` and the other using `which()`.
#        Recall whether the genes are along rows or columns of the data.
"TP53" %in% colnames(expr);
which(colnames(expr)=="TP53");

# 6b  Print the expression values of TP53 across all cell lines.
expr[,"TP53"];

# 6c  What are the main statistics?
#     (e.g. min, max, median, mean, standard deviation)
summary(expr[,"TP53"]);
# looking at them one-by-one:
min(expr[,"TP53"]);
max(expr[,"TP53"]);
median(expr[,"TP53"]);
mean(expr[,"TP53"]);
sd(expr[,"TP53"]);

# 6d  Look at the histogram and QQ-plot of TP53 expression values.
#     Are the values normally distributed?
hist(expr[,"TP53"]);
qqnorm(expr[,"TP53"]);
qqline(expr[,"TP53"], col="red");

# 6e  From the histogram, we can see that some cell lines have low TP53 
#     expression, while others have high TP53 expression.
#     Let's consider expression values of 5.5 or lower as low expression.
#     First, check which entries of the TP53 expression values are < 5.5.
expr[,"TP53"] < 5.5
#    How many cell lines with low TP53 expression are there? 
#    Use the answer above in place of the ???? in the below partial code:
nrow(expr[expr[,"TP53"] < 5.5,])    # nrow(expr[????,])

# 6f  What are the names of these cell lines?
rownames(expr[expr[,"TP53"] < 5.5,]);

# 6g  Create a scatterplot of TP53 copy-number and expression.
#     What do you notice about the trend in the values?  TP53 encodes a protein
#     involved in DNA damage response; TP53 is a tumour suppressor gene because
#     it induces programmed cell death or cell cycle arrest upon DNA damage.
plot(cn[,"TP53"], expr[,"TP53"]);

# Let's make a boxplot of TP53 expression value across different cancer sites.

# First, make a data frame of the site and T53 expression value 
tp53.df <- data.frame(
  site = pheno$site,
  expr = expr[,"TP53"]
);

# 6h  Print out the first 10 rows of this data frame.
#     Notice that some of the cell lines come from the same site (e.g. skin).
tp53.df[1:10,]

# Then, adjust the margins of the plot and make the labels horizontal
old.par <- par(mar = c(16, 4, 1, 2), las=2);

# Finally, create box plot, grouping by site
boxplot(tp53.df$expr ~ tp53.df$site);

# Reset plotting parameters back to old values
par(old.par);


# ==============================================================================
# SECTION 7  Principal component analysis
# 

# 7a  Apply PCA on the expression data.
# Note:  prcomp() expects the samples to be along the rows
expr.pr <- prcomp(expr);

# 7b  Plot the first two principal components.
plot(expr.pr$x[,1], expr.pr$x[,2]);

# 7c  What are the groups? Do they represent different cancer types?
# Hint: You can plot data points in different colours by specifiying the `col`
#       parameter using a vector of numbers.
#       The phenotype information for the cell lines are in `pheno`.
plot(expr.pr$x[,1], expr.pr$x[,2], col=as.numeric(pheno$site),
	xlim=c(-10, 90));

# 7d  Add a legend to the plot.
sites <- levels(pheno$site);
legend("topright", legend=sites, fill=1:length(sites), bty="n");

# 7e  There are too many different cancer types.
#     Let's subset a few: breast, prostate, ovary, lung, skin, bone, 
#     haematopoietic_and_lymphoid_tissue, central_nervous_system.
#     Create a vector of cell lines correpsonding to these cancer types.
cell.lines <- rownames(pheno)[pheno$site %in%
	c("breast", "prostate", "ovary", "lung", "skin", "bone",
	"haematopoietic_and_lymphoid_tissue", "central_nervous_system")];

# 7f  Create a subset of `pheno` and call it `pheno.sub`.
pheno.sub <- pheno[cell.lines, ];

# Note:  You need to re-create factor variables to remove missing factor levels.
pheno.sub$site <- factor(pheno.sub$site);

# 7g  Create a subset of `expr` and apply PCA.
expr.sub <- expr[cell.lines, ];
expr.sub.pr <- prcomp(expr.sub);

# 7h  Plot the first two principal components of the data subset.
plot(expr.sub.pr$x[,1], expr.sub.pr$x[,2], col=as.numeric(pheno.sub$site),
	xlim=c(-30, 15), ylim=c(-20, 15));
sites <- levels(pheno.sub$site);
legend("bottomleft", legend=sites, fill=1:length(sites), bty="n");
