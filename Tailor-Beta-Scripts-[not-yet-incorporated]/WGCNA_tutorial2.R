#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================
nstall.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("GO.db", "preprocessCore", "impute"))
install.packages("path/to/file", repos = NULL, lib = "path/to/library")
library
library
library
instead of the usual
library("GO.db")
library("preprocessCore")
library("impute")
library("matrixStats")
library("Hmisc")
library("splines")
library("foreach")
library("doParallel")
library("fastcluster")
library("dynamicTreeCut")
library("survival") 
library("KEGG.db")
library("topGO")
library("hgu133a.db")
library("hgu95av2.db")
library("annotate")
library("hgu133plus2.db")
library("SNPlocs.Hsapiens.dbSNP.20100427")
library("minet")
library("OrderedList")
library(WGCNA);
library("survival") 

orgCodes = c("Hs", "Mm", "Rn", "Pf", "Sc", "Dm", "Bt", "Ce", "Cf", "Dr", "Gg");
orgExtensions = c(rep(".eg", 4), ".sgd", rep(".eg", 6));
packageNames = paste("org.", orgCodes, orgExtensions, ".db", sep="");

#biocLite(c("KEGG.db", "topGO GO.db", packageNames,
		   "hgu133a.db", "hgu95av2.db", "annotate",
"hgu133plus2.db", "SNPlocs.Hsapiens.dbSNP.20100427", 
"minet", "OrderedList"))


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the package
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
femData = read.csv("LiverFemale3600.csv");
# Read in the male liver data set
maleData = read.csv("LiverMale3600.csv");
# Take a quick look at what is in the data sets (caution, longish output):
dim(femData)
names(femData)
dim(maleData)
names(maleData)


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# We work with two sets:
nSets = 2;
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("Female liver", "Male liver")
shortLabels = c("Female", "Male")
# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)

multiExpr[[1]] = list(data = as.data.frame(t(femData[-c(1:8)])));
names(multiExpr[[1]]$data) = femData$substanceBXH;
rownames(multiExpr[[1]]$data) = names(femData)[-c(1:8)];
multiExpr[[2]] = list(data = as.data.frame(t(maleData[-c(1:8)])));
names(multiExpr[[2]]$data) = maleData$substanceBXH;
rownames(multiExpr[[2]]$data) = names(maleData)[-c(1:8)];
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


pdf(file = "Plots/SampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
dev.off();


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


# Choose the "base" cut height for the female data set
baseHeight = 16
# Adjust the cut height for the male data set for the number of samples
cutHeights = c(16, 16*exprSize$nSamples[2]/exprSize$nSamples[1]);
# Re-plot the dendrograms including the cut lines
pdf(file = "Plots/SampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
{
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
  abline(h=cutHeights[set], col = "red");
}
dev.off();


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


for (set in 1:nSets)
{
  # Find clusters cut by the line
  labels = cutreeStatic(sampleTrees[[set]], cutHeight = cutHeights[set])
  # Keep the largest one (labeled by the number 1)
  keep = (labels==1)
  multiExpr[[set]]$data = multiExpr[[set]]$data[keep, ]
}
collectGarbage();
# Check the size of the leftover data
exprSize = checkSets(multiExpr)
exprSize


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


traitData = read.csv("ClinicalTraits.csv");
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData[, -c(31, 16)];
allTraits = allTraits[, c(2, 11:36) ];
# See how big the traits are and what are the trait and sample names
dim(allTraits)
names(allTraits)
allTraits$Mice
# Form a multi-set structure that will hold the clinical traits.
Traits = vector(mode="list", length = nSets);
for (set in 1:nSets)
{
  setSamples = rownames(multiExpr[[set]]$data);
  traitRows = match(setSamples, allTraits$Mice);
  Traits[[set]] = list(data = allTraits[traitRows, -1]);
  rownames(Traits[[set]]$data) = allTraits[traitRows, 1];
}
collectGarbage();
# Define data set dimensions
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


save(multiExpr, Traits, nGenes, nSamples, setLabels, shortLabels, exprSize, 
     file = "Consensus-dataInput.RData");

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. 
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "Consensus-dataInput.RData");
# The variable lnames contains the names of loaded variables.
lnames
# Get the number of sets in the multiExpr structure.
nSets = checkSets(multiExpr)$nSets


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]]);
collectGarbage();
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
"Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
pdf(file = "Plots/scaleFreeAnalysis.pdf", wi = 8, he = 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off();


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


net = blockwiseConsensusModules(
        multiExpr, power = 6, minModuleSize = 30, deepSplit = 2,
        pamRespectsDendro = FALSE, 
        mergeCutHeight = 0.25, numericLabels = TRUE,
        minKMEtoStay = 0,
        saveTOMs = TRUE, verbose = 5)


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


consMEs = net$multiMEs;
moduleLabels = net$colors;
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
consTree = net$dendrograms[[1]]; 


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


sizeGrWindow(8,6);
pdf(file = "Plots/ConsensusDendrogram-auto.pdf", wi = 8, he = 6)
plotDendroAndColors(consTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")

dev.off()


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


save(consMEs, moduleLabels, moduleColors, consTree, file = "Consensus-NetworkConstruction-auto.RData")

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. 
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "Consensus-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Get the number of sets in the multiExpr structure.
nSets = checkSets(multiExpr)$nSets


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]]);
collectGarbage();
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
"Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


softPower = 6;
# Initialize an appropriate array to hold the adjacencies
adjacencies = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate adjacencies in each individual data set
for (set in 1:nSets)
  adjacencies[set, , ] = abs(cor(multiExpr[[set]]$data, use = "p"))^softPower;


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Initialize an appropriate array to hold the TOMs
TOM = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate TOMs in each individual data set
for (set in 1:nSets)
  TOM[set, , ] = TOMsimilarity(adjacencies[set, , ]);


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


# Define the reference percentile
scaleP = 0.95
# Set RNG seed for reproducibility of sampling
set.seed(12345)
# Sample sufficiently large number of TOM entries
nSamples = as.integer(1/(1-scaleP) * 1000);
# Choose the sampled TOM entries
scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
TOMScalingSamples = list();
# These are TOM values at reference percentile
scaleQuant = rep(1, nSets)
# Scaling powers to equalize reference TOM values
scalePowers = rep(1, nSets)
# Loop over sets
for (set in 1:nSets)
{
  # Select the sampled TOM entries
  TOMScalingSamples[[set]] = as.dist(TOM[set, , ])[scaleSample]
  # Calculate the 95th percentile
  scaleQuant[set] = quantile(TOMScalingSamples[[set]],
                             probs = scaleP, type = 8);
  # Scale the male TOM
  if (set>1)
  {
    scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
    TOM[set, ,] = TOM[set, ,]^scalePowers[set];
  }
}


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# For plotting, also scale the sampled TOM entries
scaledTOMSamples = list();
for (set in 1:nSets)
  scaledTOMSamples[[set]] = TOMScalingSamples[[set]]^scalePowers[set]
# Open a suitably sized graphics window
sizeGrWindow(6,6)
pdf(file = "Plots/TOMScaling-QQPlot.pdf", wi = 6, he = 6);
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[2]], plot.it = TRUE, cex = 0.6,
                    xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[2]),
                    main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[2]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20);
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off();


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


consensusTOM = pmin(TOM[1, , ], TOM[2, , ]);


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


# Clustering
consTree = hclust(as.dist(1-consensusTOM), method = "average");
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
unmergedLabels = cutreeDynamic(dendro = consTree, distM = 1-consensusTOM,
                deepSplit = 2, cutHeight = 0.995,
                minClusterSize = minModuleSize,
                pamRespectsDendro = FALSE );
unmergedColors = labels2colors(unmergedLabels)


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


sizeGrWindow(8,6)
plotDendroAndColors(consTree, unmergedColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


# Calculate module eigengenes
unmergedMEs = multiSetMEs(multiExpr, colors = NULL, universalColors = unmergedColors)
# Calculate consensus dissimilarity of consensus module eigengenes
consMEDiss = consensusMEDissimilarity(unmergedMEs);
# Cluster consensus modules
consMETree = hclust(as.dist(consMEDiss), method = "average");
# Plot the result
sizeGrWindow(7,6)
par(mfrow = c(1,1))
plot(consMETree, main = "Consensus clustering of consensus module eigengenes",
     xlab = "", sub = "")
abline(h=0.25, col = "red")


#=====================================================================================
#
#  Code chunk 11
#
#=====================================================================================


merge = mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.25, verbose = 3)


#=====================================================================================
#
#  Code chunk 12
#
#=====================================================================================


# Numeric module labels
moduleLabels = merge$colors;
# Convert labels to colors
moduleColors = labels2colors(moduleLabels)
# Eigengenes of the new merged modules:
consMEs = merge$newMEs;


#=====================================================================================
#
#  Code chunk 13
#
#=====================================================================================


sizeGrWindow(9,6)
plotDendroAndColors(consTree, cbind(unmergedColors, moduleColors),
                    c("Unmerged", "Merged"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#=====================================================================================
#
#  Code chunk 14
#
#=====================================================================================


save(consMEs, moduleColors, moduleLabels, consTree, file = "Consensus-NetworkConstruction-man.RData")

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. 
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "Consensus-dataInput.RData");
# The variable lnames contains the names of loaded variables.
lnames
# Get the number of sets in the multiExpr structure.
nSets = checkSets(multiExpr)$nSets


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]]);
collectGarbage();
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
"Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
pdf(file = "Plots/scaleFreeAnalysis.pdf", wi = 8, he = 6);
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off();


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


bnet = blockwiseConsensusModules(
    multiExpr, maxBlockSize = 2000, power = 6, minModuleSize = 30,
    deepSplit = 2, 
    pamRespectsDendro = FALSE, 
    mergeCutHeight = 0.25, numericLabels = TRUE,
    minKMEtoStay = 0,
    saveTOMs = TRUE, verbose = 5)


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


load(file = "Consensus-NetworkConstruction-auto.RData")
bwLabels = matchLabels(bnet$colors, moduleLabels, pThreshold = 1e-7);
bwColors = labels2colors(bwLabels)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


# Here we show a more flexible way of plotting several trees and colors on one page
sizeGrWindow(12,6)
pdf(file = "Plots/BlockwiseGeneDendrosAndColors.pdf", wi = 12, he = 6);
# Use the layout function for more involved screen sectioning
layout(matrix(c(1:4), 2, 2), heights = c(0.8, 0.2), widths = c(1,1))
#layout.show(4);
nBlocks = length(bnet$dendrograms)
# Plot the dendrogram and the module colors underneath for each block
for (block in 1:nBlocks)
  plotDendroAndColors(bnet$dendrograms[[block]], moduleColors[bnet$blockGenes[[block]]],
                      "Module colors", 
                      main = paste("Gene dendrogram and module colors in block", block), 
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      setLayout = FALSE)
dev.off();


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


sizeGrWindow(12,9);
pdf(file="Plots/SingleDendro-BWColors.pdf", wi = 12, he = 9);
plotDendroAndColors(consTree,
                    cbind(moduleColors, bwColors),
                    c("Single block", "Blockwise"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Single block consensus gene dendrogram and module colors")
dev.off();

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the data saved in the first part
lnames = load(file = "Consensus-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load the results of network analysis, tutorial part 2.a
lnames = load(file = "Consensus-NetworkConstruction-auto.RData");
lnames


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


lnames = load("../Mouse-Female/FemaleLiver-02-networkConstruction-auto.RData")
lnames
# Rename variables to avoid conflicts
femaleLabels = moduleLabels;
femaleColors = moduleColors;
femaleTree = geneTree;
femaleMEs = orderMEs(MEs, greyName = "ME0");


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


lnames = load("Consensus-NetworkConstruction-auto.RData")
lnames


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Isolate the module labels in the order they appear in ordered module eigengenes
femModuleLabels = substring(names(femaleMEs), 3)
consModuleLabels = substring(names(consMEs[[1]]$data), 3)
# Convert the numeric module labels to color labels
femModules = labels2colors(as.numeric(femModuleLabels))
consModules = labels2colors(as.numeric(consModuleLabels))
# Numbers of female and consensus modules
nFemMods = length(femModules)
nConsMods = length(consModules)
# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = nFemMods, ncol = nConsMods);
CountTbl = matrix(0, nrow = nFemMods, ncol = nConsMods);
# Execute all pairwaise comparisons
for (fmod in 1:nFemMods)
  for (cmod in 1:nConsMods)
  {
    femMembers = (femaleColors == femModules[fmod]);
    consMembers = (moduleColors == consModules[cmod]);
    pTable[fmod, cmod] = -log10(fisher.test(femMembers, consMembers, alternative = "greater")$p.value);
    CountTbl[fmod, cmod] = sum(femaleColors == femModules[fmod] & moduleColors ==
                      consModules[cmod])
  }


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


# Truncate p values smaller than 10^{-50} to 10^{-50} 
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50 ;
# Marginal counts (really module sizes)
femModTotals = apply(CountTbl, 1, sum)
consModTotals = apply(CountTbl, 2, sum)
# Actual plotting
sizeGrWindow(10,7 );
pdf(file = "Plots/ConsensusVsFemaleModules.pdf", wi = 10, he = 7);
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", consModules),
               yLabels = paste(" ", femModules),
               colorLabels = TRUE,
               xSymbols = paste("Cons ", consModules, ": ", consModTotals, sep=""),
               ySymbols = paste("Fem ", femModules, ": ", femModTotals, sep=""),
               textMatrix = CountTbl,
               colors = greenWhiteRed(100)[50:100],
               main = "Correspondence of Female set-specific and Female-Male consensus modules",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE);
dev.off();

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the data saved in the first part
lnames = load(file = "Consensus-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Also load results of network analysis
lnames = load(file = "Consensus-NetworkConstruction-auto.RData");
lnames
exprSize = checkSets(multiExpr);
nSets = exprSize$nSets;


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Set up variables to contain the module-trait correlations
moduleTraitCor = list();
moduleTraitPvalue = list();
# Calculate the correlations
for (set in 1:nSets)
{
  moduleTraitCor[[set]] = cor(consMEs[[set]]$data, Traits[[set]]$data, use = "p");
  moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set]);
}


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


# Convert numerical lables to colors for labeling of modules in the plot
MEColors = labels2colors(as.numeric(substring(names(consMEs[[1]]$data), 3)));
MEColorNames = paste("ME", MEColors, sep="");
# Open a suitably sized window (the user should change the window size if necessary)
sizeGrWindow(10,7)
pdf(file = "Plots/ModuleTraitRelationships-female.pdf", wi = 10, he = 7);
# Plot the module-trait relationship table for set number 1
set = 1
textMatrix =  paste(signif(moduleTraitCor[[set]], 2), "\n(",
                           signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module--trait relationships in", setLabels[set]))
dev.off();
# Plot the module-trait relationship table for set number 2
set = 2
textMatrix =  paste(signif(moduleTraitCor[[set]], 2), "\n(",
                           signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
sizeGrWindow(10,7)
pdf(file = "Plots/ModuleTraitRelationships-male.pdf", wi = 10, he = 7);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module--trait relationships in", setLabels[set]))
dev.off();


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Initialize matrices to hold the consensus correlation and p-value
consensusCor = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));
consensusPvalue = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));
# Find consensus negative correlations
negative = moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0;
consensusCor[negative] = pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative]);
consensusPvalue[negative] = pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative]);
# Find consensus positive correlations
positive = moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0;
consensusCor[positive] = pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive]);
consensusPvalue[positive] = pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive]);


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


textMatrix =  paste(signif(consensusCor, 2), "\n(",
                           signif(consensusPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
sizeGrWindow(10,7)
pdf(file = "Plots/ModuleTraitRelationships-consensus.pdf", wi = 10, he = 7);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = consensusCor,
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Consensus module--trait relationships across\n",
                            paste(setLabels, collapse = " and ")))


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


file = gzfile(description = "GeneAnnotation.csv.gz");
annot = read.csv(file = file);
# Match probes in the data set to the probe IDs in the annotation file 
probes = names(multiExpr[[1]]$data)
probes2annot = match(probes, annot$substanceBXH)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


consMEs.unord = multiSetMEs(multiExpr, universalColors = moduleLabels, excludeGrey = TRUE)
GS = list();
kME = list();
for (set in 1:nSets)
{
  GS[[set]] = corAndPvalue(multiExpr[[set]]$data, Traits[[set]]$data);
  kME[[set]] = corAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data);
}


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


GS.metaZ = (GS[[1]]$Z + GS[[2]]$Z)/sqrt(2);
kME.metaZ = (kME[[1]]$Z + kME[[2]]$Z)/sqrt(2);
GS.metaP = 2*pnorm(abs(GS.metaZ), lower.tail = FALSE);
kME.metaP = 2*pnorm(abs(kME.metaZ), lower.tail = FALSE);


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


GSmat = rbind(GS[[1]]$cor, GS[[2]]$cor, GS[[1]]$p, GS[[2]]$p, GS.metaZ, GS.metaP);
nTraits = checkSets(Traits)$nGenes
traitNames = colnames(Traits[[1]]$data)
dim(GSmat) = c(nGenes, 6*nTraits)
rownames(GSmat) = probes;
colnames(GSmat) = spaste(
    c("GS.set1.", "GS.set2.", "p.GS.set1.", "p.GS.set2.", "Z.GS.meta.", "p.GS.meta"),
    rep(traitNames, rep(6, nTraits)))
# Same code for kME:
kMEmat = rbind(kME[[1]]$cor, kME[[2]]$cor, kME[[1]]$p, kME[[2]]$p, kME.metaZ, kME.metaP);
MEnames = colnames(consMEs.unord[[1]]$data);
nMEs = checkSets(consMEs.unord)$nGenes
dim(kMEmat) = c(nGenes, 6*nMEs)
rownames(kMEmat) = probes;
colnames(kMEmat) = spaste(
    c("kME.set1.", "kME.set2.", "p.kME.set1.", "p.kME.set2.", "Z.kME.meta.", "p.kME.meta"),
    rep(MEnames, rep(6, nMEs)))


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


info = data.frame(Probe = probes, GeneSymbol = annot$gene_symbol[probes2annot],
             EntrezID = annot$LocusLinkID[probes2annot],
             ModuleLabel = moduleLabels,
             ModuleColor = labels2colors(moduleLabels),
             GSmat,
             kMEmat);
write.csv(info, file = "consensusAnalysis-CombinedNetworkResults.csv",
          row.names = FALSE, quote = FALSE);


#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Basic settings: we work with two data sets
nSets = 2;
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("Female liver", "Male liver")
shortLabels = c("Female", "Male")
# Load the data saved in the first part
lnames = load(file = "Consensus-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load the results of network analysis, tutorial part 2.a
lnames = load(file = "Consensus-NetworkConstruction-auto.RData");
lnames


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Create a variable weight that will hold just the body weight of mice in both sets
weight = vector(mode = "list", length = nSets);
for (set in 1:nSets)
{
  weight[[set]] = list(data = as.data.frame(Traits[[set]]$data$weight_g));
  names(weight[[set]]$data) = "weight"
}
# Recalculate consMEs to give them color names
consMEsC = multiSetMEs(multiExpr, universalColors = moduleColors);
# We add the weight trait to the eigengenes and order them by consesus hierarchical clustering:
MET = consensusOrderMEs(addTraitToMEs(consMEsC, weight));


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


sizeGrWindow(8,10);
pdf(file = "Plots/EigengeneNetworks.pdf", width= 8, height = 10);
par(cex = 0.9)
plotEigengeneNetworks(MET, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
                      zlimPreservation = c(0.5, 1), xLabelsAngle = 90)
# Version: 03/14/2008
#   Change:  ModuleEnrichment1: The version is more robust: it avoids and error when only grey modules are present (i.e. no modules have been found).
#
# Version: 03/13/2008
#   Change:  PlotClusterTreeSamples: a new function that produces a cluster tree of the microarray samples possibly colored by a sample trait.
#
# Version: 02/18/2008
#   Change:  signedKME: new option for customizing the output column names.
#            The user can either choose the old labels "kME" (or the new labels module membership "MM."). See the option outputColumnName="kME"
#
# Version: 02/13/2008
#   Change:  signedKME: now it can handle missing data and constant gene expressions.
#
# Version: 02/04/2008
#   Change:  ModulePrinComps1: a comment on "FiveComponents=T" was removed.
#            SoftConnectivity: a criterion for checking "no.genes < no.samples" was removed. 
#
# Version: 01/07/2008
#   Change:  ModulePrinComps1: allows that the user does not input a color vector.
#
# Version: 11/16/2007
#   Change: \! changed to \n in modulecolor2 to get rid of a warning message.
#   No functional changes.

# Version: 09/21/2007
#   Change: ModulePrinComps1: code checks the structure of value returned by
#   impute.knn to make it compatible both with the CRAN version 1.0-5 and
#   Bioconductor version 1.8.

# This document contains functions that we find useful for gene co-expression network analysis
# We recommend to read the tutorial first.
# The functions were written by Steve Horvath, Jun Dong, Peter Langfelder, Andy Yip, Bin Zhang
# To cite this code or the statistical methods please use the following 2 references:

#1) Horvath S, Zhang B, Carlson M, Lu KV, Zhu S, Felciano RM, Laurance MF, Zhao W, Shu, Q, Lee Y, 
# Scheck AC, Liau LM, Wu H, Geschwind DH, Febbo PG, Kornblum HI, Cloughesy TF, Nelson SF, Mischel PS (2006) 
#"Analysis of Oncogenic Signaling Networks in Glioblastoma Identifies ASPM as a Novel Molecular Target", 
# PNAS | November 14, 2006 | vol. 103 | no. 46 | 17402-17407

#2) Zhang B, Horvath S (2005) "A General Framework for Weighted Gene Co-Expression Network Analysis", 
#Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17 


# Technical Report and software code at: www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork.

# CONTENTS  
# This document contains function for carrying out the following tasks
# A) Assessing scale free topology and choosing the parameters of the adjacency function
#    using the scale free topology criterion (Zhang and Horvath 05)
# B) Computing the topological overlap matrix 
# C) Defining gene modules using clustering procedures
# D) Summing up modules by their first principal component (first eigengene)
# E) Relating a measure of gene significance to the modules 
# F) Carrying out a within module analysis (computing intramodular connectivity) 
#    and relating intramodular connectivity to gene significance.
# G) Miscellaneous other functions, e.g. for computing the cluster coefficient.
# H) Network functions from Bin Zhang and Peter Langfelder for dynamic tree cutting 
     #of a hierarchical clustering tree (dendrogram)
# I) Peter Langfelder's functions for merging modules, e.g. when their module eigengenes have a high correlation
# J) General statistical functions, e.g. for scatterplots

# The functions below are organized according into categories A-J.

#####################################################################################################
################################################################################################################################
# A) Assessing scale free topology and choosing the parameters of the adjacency function
#    using the scale free topology criterion (Zhang and Horvath 05)
#################################################################################################################


# ===================================================
#For hard thresholding, we use the signum (step) function
if(exists("signum1") ) rm(signum1); 
signum1=function(corhelp,tau1)  {
adjmat1= as.matrix(abs(corhelp)>=tau1)
dimnames(adjmat1) <- dimnames(corhelp)
diag(adjmat1) <- 0
adjmat1}

# ===================================================
# For soft thresholding, one can use the sigmoid function 
# But we have focused on the power adjacency function in the tutorial...
if (exists("sigmoid1") ) rm(sigmoid1); sigmoid1=function(ss,mu1=0.8,alpha1=20) {
1/(1+exp(-alpha1*(ss-mu1)))}





#This function is useful for speeding up the connectivity calculation.
#The idea is to partition the adjacency matrix into consecutive baches of a #given size.
#In principle, the larger the batch size the faster is the calculation. But #smaller batchsizes require #less memory...
# Input: gene expression data set where *rows* correspond to microarray samples #and columns correspond to genes. 
# If fewer than MinimumNoSamples contain gene expression information for a given
# gene, then its connectivity is set to missing. 
if(exists("SoftConnectivity")) rm(SoftConnectivity);
SoftConnectivity=function(datE, power=6,batchsize=1500,MinimumNoSamples=10) {
no.genes=dim(datE)[[2]]
no.samples=dim(datE)[[1]]
if (no.genes<10 | no.samples<5 ) {print("Error: Something seems to be wrong. Make sure that the input data frame has genes as rows and array samples as columns. Alternatively, there could be fewer than 10 genes or fewer than 5 samples. ") } else {
if (no.genes<no.samples ) print("Warning: There are fewer genes than samples in the function SoftConnectivity. Maybe you should transpose the data?")
sum1=function(x) sum(x,na.rm=T)
k=rep(NA,no.genes)
no.batches=as.integer(no.genes/ batchsize)
if (no.batches>0) {
for (i in 1:no.batches) {
print(paste("batch number = ", i))
index1=c(1:batchsize)+(i-1)* batchsize
ad1=abs(cor(datE[,index1], datE,use="p"))^power
ad1[is.na(ad1)]=0
k[index1]=apply(ad1,1,sum1)
# If fewer than MinimumNoSamples contain gene expression information for a given
# gene, then we set its connectivity to 0.
NoSamplesAvailable=apply(!is.na(datE[,index1]),2,sum)
k[index1][NoSamplesAvailable< MinimumNoSamples]=NA
} # end of for (i in 1:no.batches
} # end of if (no.batches>0)...
if (no.genes-no.batches*batchsize>0 ) {
restindex=c((no.batches*batchsize+1):no.genes)
ad1=abs(cor(datE[,restindex], datE,use="p"))^power
ad1[is.na(ad1)]=0
k[restindex]=apply(ad1,1,sum1)
NoSamplesAvailable=apply(!is.na(datE[,restindex]),2,sum)
k[restindex][NoSamplesAvailable< MinimumNoSamples]=NA
} # end of if
} # end of else statement
k
} # end of function





# ===================================================
# The function PickHardThreshold can help one to estimate the cut-off value 
# when using the signum (step) function.
# The first column lists the threshold ("cut"),
# the second column lists the corresponding p-value based on the Fisher Transform 
# of the correlation. 
# The third column reports the resulting scale free topology fitting index R^2.
# The fourth column reports the slope of the fitting line, it shoud be negative for 
# biologically meaningul networks.
# The fifth column reports the fitting index for the truncated exponential model. 
# Usually we ignore it.
# The remaining columns list the mean, median and maximum resulting connectivity.
# To pick a hard threshold (cut) with the scale free topology criterion:
# aim for high scale free R^2 (column 3), high connectivity (col 6) and negative slope 
# (around -1, col 4).
# The output is a list with 2 components. The first component lists a sugggested cut-off
# while the second component contains the whole table.
# The removeFirst option removes the first point (k=0, P(k=0)) from the regression fit.
# no.breaks specifies how many intervals used to estimate the frequency p(k) i.e. the no. of points in the 
# scale free topology plot.
if (exists("PickHardThreshold")) rm(PickHardThreshold);
PickHardThreshold=function(datExpr1,RsquaredCut=0.85, cutvector=seq(0.1,0.9,by=0.05) ,removeFirst=FALSE,no.breaks=10) {
no.genes   <- dim(datExpr1)[[2]]
no.genes <- dim(datExpr1)[[2]]
no.samples= dim(datExpr1)[[1]]
colname1=c("Cut","p-value", "scale law R^2", "slope="  ,"truncated R^2","mean(k)","median(k)","max(k)")
datout=data.frame(matrix(666,nrow=length(cutvector),ncol=length(colname1) ))
names(datout)=colname1
datout[,1]=cutvector
for (i in c(1:length(cutvector) ) ){
cut1=cutvector[i]
datout[i,2]=2*(1-pt(sqrt(no.samples-1)*cut1/sqrt(1-cut1^2),no.samples-1))}
if(exists("fun1")) rm(fun1)
fun1=function(x) {
corx=abs(cor(x,datExpr1,use="p"))
out1=rep(NA, length(cutvector) )
for (j in c(1:length(cutvector))) {out1[j]=sum(corx>cutvector[j])}
out1
} # end of fun1
datk=t(apply(datExpr1,2,fun1))
for (i in c(1:length(cutvector) ) ){
nolinkshelp <- datk[,i]-1
cut2=cut(nolinkshelp,no.breaks)
binned.k=tapply(nolinkshelp,cut2,mean)
freq1=as.vector(tapply(nolinkshelp,cut2,length)/length(nolinkshelp))
# The following code corrects for missing values etc
breaks1=seq(from=min(nolinkshelp),to=max(nolinkshelp),length=no.breaks+1)
hist1=hist(nolinkshelp,breaks=breaks1,equidist=F,plot=FALSE,right=TRUE)
binned.k2=hist1$mids
binned.k=ifelse(is.na(binned.k),binned.k2,binned.k)
binned.k=ifelse(binned.k==0,binned.k2,binned.k)
freq1=ifelse(is.na(freq1),0,freq1)
xx= as.vector(log10(binned.k))
if(removeFirst) {freq1=freq1[-1]; xx=xx[-1]}
plot(xx,log10(freq1+.000000001),xlab="log10(k)",ylab="log10(p(k))" )
lm1= lm(as.numeric(log10(freq1+.000000001))~ xx )
lm2=lm(as.numeric(log10(freq1+.000000001))~ xx+I(10^xx) )
datout[i,3]=summary(lm1)$adj.r.squared 
datout[i,4]=summary(lm1)$coefficients[2,1]  
datout[i,5]=summary(lm2)$adj.r.squared
datout[i,6]=mean(nolinkshelp)
datout[i,7]= median(nolinkshelp)
datout[i,8]= max(nolinkshelp) 
} 
datout=signif(datout,3) 
print(data.frame(datout));
# the cut-off is chosen as smallest cut with R^2>RsquaredCut 
ind1=datout[,3]>RsquaredCut
indcut=NA
indcut=ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)
# this estimates the cut-off value that should be used. 
# Don't trust it. You need to consider slope and mean connectivity as well!
cut.estimate=cutvector[indcut][[1]]
list(cut.estimate, data.frame(datout));
} # end of function











# ===========================================================
# The function PickSoftThreshold allows one to estimate the power parameter when using
# a soft thresholding approach with the use of the power function AF(s)=s^Power
# The function PickSoftThreshold allows one to estimate the power parameter when using
# a soft thresholding approach with the use of the power function AF(s)=s^Power
# The removeFirst option removes the first point (k=1, P(k=1)) from the regression fit.
if (exists("PickSoftThreshold")) rm(PickSoftThreshold);
PickSoftThreshold=function(datExpr1,RsquaredCut=0.85, powervector=c(seq(1,10,by=1),seq(12,20,by=2)),
removeFirst=FALSE,no.breaks=10) {
no.genes   <- dim(datExpr1)[[2]]
no.genes <- dim(datExpr1)[[2]]
no.samples= dim(datExpr1)[[1]]
colname1=c("Power", "scale law R^2" ,"slope", "truncated R^2","mean(k)","median(k)","max(k)")
datout=data.frame(matrix(666,nrow=length(powervector),ncol=length(colname1) ))
names(datout)=colname1
datout[,1]=powervector
if(exists("fun1")) rm(fun1)
fun1=function(x) {
corx=abs(cor(x,datExpr1,use="p"))
out1=rep(NA, length(powervector) )
for (j in c(1:length(powervector))) {out1[j]=sum(corx^powervector[j])}
out1
} # end of fun1
datk=t(apply(datExpr1,2,fun1))
for (i in c(1:length(powervector) ) ){
nolinkshelp <- datk[,i]-1
cut2=cut(nolinkshelp,no.breaks)
binned.k=tapply(nolinkshelp,cut2,mean)
freq1=as.vector(tapply(nolinkshelp,cut2,length)/length(nolinkshelp))
# The following code corrects for missing values etc
breaks1=seq(from=min(nolinkshelp),to=max(nolinkshelp),length=no.breaks+1)
hist1=hist(nolinkshelp,breaks=breaks1,equidist=F,plot=FALSE,right=TRUE)
binned.k2=hist1$mids
binned.k=ifelse(is.na(binned.k),binned.k2,binned.k)
binned.k=ifelse(binned.k==0,binned.k2,binned.k)
freq1=ifelse(is.na(freq1),0,freq1)

xx= as.vector(log10(binned.k))
if(removeFirst) {freq1=freq1[-1]; xx=xx[-1]}
plot(xx,log10(freq1+.000000001),xlab="log10(k)",ylab="log10(p(k))" )
lm1= lm(as.numeric(log10(freq1+.000000001))~ xx )
lm2=lm(as.numeric(log10(freq1+.000000001))~ xx+I(10^xx) )
datout[i,2]=summary(lm1)$adj.r.squared 
datout[i,3]=summary(lm1)$coefficients[2,1]  
datout[i,4]=summary(lm2)$adj.r.squared
datout[i,5]=mean(nolinkshelp)
datout[i,6]= median(nolinkshelp)
datout[i,7]= max(nolinkshelp) 
} 
datout=signif(datout,3) 
print(data.frame(datout));
# the cut-off is chosen as smallest cut with R^2>RsquaredCut 
ind1=datout[,2]>RsquaredCut
indcut=NA
indcut=ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)
# this estimates the power value that should be used. 
# Don't trust it. You need to consider slope and mean connectivity as well!
power.estimate=powervector[indcut][[1]]
list(power.estimate, data.frame(datout));
}






# ===================================================
# The function ScaleFreePlot1 creates a plot for checking scale free topology
# when truncated1=T is specificed, it provides the R^2 measures for the following
# degree distributions: a) scale free topology, b) log-log R^2 and c) truncated exponential R^2

# The function ScaleFreePlot1 creates a plot for checking scale free topology
if(exists("ScaleFreePlot1")) rm(ScaleFreePlot1) ; 
ScaleFreePlot1=function(kk,no.breaks=10,AF1="" ,truncated1=FALSE, removeFirst=FALSE,cex.lab1=1){
cut1=cut(kk,no.breaks)
binned.k=tapply(kk,cut1,mean)
freq1=tapply(kk,cut1,length)/length(kk)
# The following code corrects for missing values etc
breaks1=seq(from=min(kk),to=max(kk),length=no.breaks+1)
hist1=hist(kk,breaks=breaks1,equidist=F,plot=FALSE,right=TRUE)
binned.k2=hist1$mids
binned.k=ifelse(is.na(binned.k),binned.k2,binned.k)
binned.k=ifelse(binned.k==0,binned.k2,binned.k)
freq1=ifelse(is.na(freq1),0,freq1)
plot(log10(binned.k),log10(freq1+.000000001),xlab="log10(k)",ylab="log10(p(k))",cex.lab=cex.lab1 )
xx= as.vector(log10(binned.k))
if(removeFirst) {freq1=freq1[-1]; xx=xx[-1]}
lm1=lm(as.numeric(log10(freq1+.000000001))~ xx )
lines(xx,predict(lm1),col=1)
OUTPUT=data.frame(ScaleFreeRsquared=round(summary(lm1)$adj.r.squared,2),Slope=round(lm1$coefficients[[2]],2))
if (truncated1==TRUE) { 
lm2=lm(as.numeric(log10(freq1+.000000001))~ xx+I(10^xx) );
OUTPUT=data.frame(ScaleFreeRsquared=round(summary(lm1)$adj.r.squared,2),Slope=round(lm1$coefficients[[2]],2),
TruncatedRsquared=round(summary(lm2)$adj.r.squared,2))
print("the red line corresponds to the truncated exponential fit")
lines(xx,predict(lm2),col=2);
title(paste(AF1, 
", scale free R^2=",as.character(round(summary(lm1)$adj.r.squared,2)),
", slope=", round(lm1$coefficients[[2]],2),
", trunc.R^2=",as.character(round(summary(lm2)$adj.r.squared,2))))} else { 
title(paste(AF1, ", scale R^2=",as.character(round(summary(lm1)$adj.r.squared,2)) , 
", slope=", round(lm1$coefficients[[2]],2)))
}
OUTPUT
} # end of function 









#################################################################################################################
################################################################################################################################
# B) Computing the topological overlap matrix 
#################################################################################################################
#################################################################################################################



# ===================================================
#The function TOMdist1 computes a dissimilarity 
# based on the topological overlap matrix (Ravasz et al)
# Input: an Adjacency matrix with entries in [0,1]
if(exists("TOMdist1")) rm(TOMdist1);
TOMdist1=function(adjmat1, maxADJ=FALSE) {
diag(adjmat1)=0;
adjmat1[is.na(adjmat1)]=0;
maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) ); 
if (maxh1>1 | minh1 < 0 ) {print(paste("ERROR: the adjacency matrix contains entries that are larger than 1 or smaller than 0!!!, max=",maxh1,", min=",minh1)) } else { 
if (  max(c(as.dist(abs(adjmat1-t(adjmat1)))))>10^(-12)   ) {print("ERROR: non-symmetric adjacency matrix!!!") } else { 
adjmat1= (adjmat1+ t(adjmat1) )/2
kk=apply(adjmat1,2,sum)
maxADJconst=1
if (maxADJ==TRUE) maxADJconst=max(c(as.dist(adjmat1 ))) 
Dhelp1=matrix(kk,ncol=length(kk),nrow=length(kk))
denomTOM= pmin(as.dist(Dhelp1),as.dist(t(Dhelp1)))   +as.dist(maxADJconst-adjmat1); 
gc();gc();
numTOM=as.dist(adjmat1 %*% adjmat1 +adjmat1);
#TOMmatrix=numTOM/denomTOM
# this turns the TOM matrix into a dissimilarity 
out1=1-as.matrix(numTOM/denomTOM) 
diag(out1)=1 
# setting the diagonal to 1 is unconventional (it should be 0)
# but it leads to nicer looking TOM plots... 
out1
}}}





# ===================================================
# This function computes a TOMk dissimilarity
# which generalizes the topological overlap matrix (Ravasz et al)
# Input: an Adjacency matrix with entries in [0,1]
# WARNING:  ONLY FOR UNWEIGHTED NETWORKS, i.e. the adjacency matrix contains binary entries...
# This function is explained in Yip and Horvath (2005)
# http://www.genetics.ucla.edu/labs/horvath/GTOM/
if(exists("TOMkdist1")) rm(TOMkdist1);
TOMkdist1 = function(adjmat1,k=1){
    maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) );
    if (k!=round(abs(k))) {
        stop("k must be a positive integer!!!", call.=TRUE);}
    if (maxh1>1 | minh1 < 0 ){
        print(paste("ERROR: entries of the adjacency matrix must be between inclusively 0 and 1!!!, max=",maxh1,", min=",minh1))}
    else {
if (  max(c(as.dist(abs(adjmat1-t(adjmat1)))))>0   ) {print("ERROR: non-symmetric adjacency matrix!!!") } else { 

        B <- adjmat1;
        if (k>=2) {
            for (i in 2:k) {
                diag(B) <- diag(B) + 1;
                B = B %*% adjmat1;}}   # this gives the number of paths with length at most k connecting a pair
        B <- (B>0);   # this gives the k-step reachability from a node to another
        diag(B) <- 0;   # exclude each node being its own neighbor
        B <- B %*% B   # this gives the number of common k-step-neighbor that a pair of nodes share

        Nk <- diag(B);
        B <- B +adjmat1;   # numerator
        diag(B) <- 1;
        denomTOM=outer(Nk,Nk,FUN="pmin")+1-adjmat1;
        diag(denomTOM) <- 1;
        1 - B/denomTOM   # this turns the TOM matrix into a dissimilarity
}}
}


# IGNORE THIS function...
# The function TOMdistROW computes the TOM distance of a gene (node)
# with that of all other genes in the network.
# WhichRow is an integer that specifies which row of the adjacency matrix
# corresponds to the gene of interest.
# Output=vector of TOM distances.
if (exists("TOMdistROW") ) rm(TOMdistROW) 
TOMdistROW=function(WhichRow=1, adjmat1, maxADJ=FALSE) {
diag(adjmat1)=0;
maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) ); 
if (maxh1>1 | minh1 < 0 ) {print(paste("ERROR: the adjacency matrix contains entries that are larger than 1 or smaller than 0!!!, max=",maxh1,", min=",minh1)) } else { 
kk=apply(adjmat1,2,sum)
numTOM=adjmat1[WhichRow,] %*% adjmat1 +adjmat1[WhichRow,]; 
numTOM[WhichRow]=1
maxADJconst=1
if (maxADJ==TRUE) maxADJconst=max(c(as.dist(adjmat1 ))) 
denomTOM=pmin(kk[WhichRow],kk)+maxADJconst-adjmat1[WhichRow,]; denomTOM[WhichRow]=1
#TOMmatrix=numTOM/denomTOM
# this turns the TOM matrix into a dissimilarity 
1-numTOM/denomTOM 
}
}


#####################################################################################################
################################################################################################################################
# C) Defining gene modules using clustering procedures
#####################################################################################################
################################################################################################################################


# ===================================================
#The function modulecolor2 function assigns colors to the observations 
# in the branches of a dendrogram. Now enhanced to handle all possible R colors. 
# Needs GlobalStandardColors.

# we use it to define modules....

if (exists("modulecolor2")) rm(modulecolor2);
modulecolor2=function(hier1, h1=0.9,minsize1=50) {
# here we define modules by using a height cut-off for the branches
labelpred= cutree(hier1,h=h1)
sort1=-sort(-table(labelpred))
modulename= as.numeric(names(sort1))
modulebranch= sort1>minsize1
no.modules=sum(modulebranch)
# now we assume that there are fewer than a certain number of colors
colorcode=GlobalStandardColors
# "grey" means not in any module;
colorhelp=rep("grey",length(labelpred))
if ( no.modules==0 | no.modules >length(colorcode)){ print(paste("The number of modules is problematic! \n Number of modules = ", as.character(no.modules)))} else { for (i in c(1:no.modules)) {colorhelp=ifelse(labelpred==modulename[i],colorcode[i],colorhelp)};
colorhelp=factor(colorhelp,levels=c(colorcode[1:no.modules],"grey"))
}
factor(colorhelp, levels=unique(colorhelp[hier1$order] ))
}





#------------------------------------------------------------------------------------------------
# Create a barplot of colors given in the couleur parameter ordered according
# to the order of the hierarchical tree hier1. 
# Parameters:
#    hier1:		Hierarchical tree (such as returned by hclust).
#    couleur:		Either a vector of colors for each element in the tree, or a matrix in which every row 
#  			is a set of colors for the corresponding object in the tree. If a matrix, the colors
#			will be plotted starting at the	bottom. 
#    RowLabels:		If given, each color row will be labeled by the corresponding entry in RowLabels; 
#			otherwise the number of the row will be used as a label.
#    cex.RowLabels:	Scale factor for the font size used for labeling rows.
#
# Return value:	None.
#    
# hclustplot1.old implements the case where couleur is a vector, hclustplotn
# the case of a matrix; hclustplot1 is a common wrapper.






#------------------------------------------------------------------------------------------------
#
# hclustplot1, hclustplot1.old, hclustplotn
# The funtion hclusplotn was created by Peter Langfelder.
#------------------------------------------------------------------------------------------------
# Create a barplot of colors given in the couleur parameter ordered according
# to the order of the hierarchical tree hier1. 
# Parameters:
#    hier1:		Hierarchical tree (such as returned by hclust).
#    couleur:		Either a vector of colors for each element in the tree, or a matrix in which every row 
#  			is a set of colors for the corresponding object in the tree. If a matrix, the colors
#			will be plotted starting at the	bottom. 
#    RowLabels:		If given, each color row will be labeled by the corresponding entry in RowLabels; 
#			otherwise the number of the row will be used as a label.
#    cex.RowLabels:	Scale factor for the font size used for labeling rows.
#
# Return value:	None.
#    
# hclustplot1.old implements the case where couleur is a vector, hclustplotn
# the case of a matrix; hclustplot1 is a common wrapper.

if (exists("hclustplot1")) rm(hclustplot1);
hclustplot1 = function(hier1, couleur, title1="Colors sorted by hierarchical clustering", 
                     RowLabels=NULL, cex.RowLabels = 0.9, ...)
{
  if (length(dim(couleur))>1)
  {
    hclustplotn(hier1, couleur, RowLabels, cex.RowLabels, main = title1, ...);
  } else
  { 
    hclustplot1.old(hier1,couleur,title1);
  }
}

if (exists("hclustplot1.old")) rm(hclustplot1.old);
hclustplot1.old=function(hier1,couleur,title1="Colors sorted by hierarchical clustering") {
if (length(hier1$order) != length(couleur) ) {print("ERROR: length of color vector not compatible with no. of objects in the hierarchical tree")};
if (length(hier1$order) == length(couleur) ) {
barplot(height=rep(1, length(couleur)), col= as.character(couleur[hier1$order]),border=F, main=title1,space=0, axes=F)}
}
 
 
# This function extends hclustplot1 by allows Color to be a data frame of different colors
if (exists("hclustplotn") ) rm(hclustplotn);
hclustplotn=function(hier1, Color, RowLabels=NULL, cex.RowLabels = 0.9, ...) 
{
if ( is.null(RowLabels) ) RowLabels=names(data.frame(Color))
  options(stringsAsFactors=FALSE);
  if (length(hier1$order) != dim(as.matrix(Color))[[1]] ) 
  { 
    stop("ERROR: length of color vector not compatible with no. of objects in the hierarchical tree");
  } else {
       No.Sets = dim(as.matrix(Color))[[2]];
       C = as.matrix(Color)[hier1$order, ]; 
       step = 1/dim(as.matrix(Color))[[1]];
       ystep = 1/No.Sets;
       barplot(height=1, col = "white", border=F,space=0, axes=F, ...)
       for (j in 1:No.Sets)
       {
         for (i in 1:dim(as.matrix(C))[[1]])
         { 
           #lines(x=rep((i*step), times=2), y=c(ystep*(j-1),ystep*j),  col = as.character(as.matrix(C)[i,j]));
           rect(xleft=((i-1)*step), ybottom=ystep*(j-1), xright = (i) * step, ytop = ystep*j,  
                    col = as.character(as.matrix(C)[i,j]), border = as.character(as.matrix(C)[i,j]));
         } 
         if (is.null(RowLabels))
         {
             text(as.character(j), pos=2, x=0, y=ystep*(j-0.5), cex=cex.RowLabels, xpd = TRUE);
         } else {
             text(RowLabels[j], pos=2, x=0, y=ystep*(j-0.5), cex=cex.RowLabels, xpd = TRUE);
         }
       }
       for (j in 1:No.Sets) lines(x=c(0,1), y=c(ystep*j,ystep*j));
  }
}

 




# ===================================================
# The function TOMplot1 creates a TOM plot
# Inputs:  distance measure, hierarchical (hclust) object, color label=couleur
if (exists("TOMplot1")) rm(TOMplot1);
TOMplot1=function(disttom,hier1, couleur,terrainColors=FALSE) {
no.nodes=length(couleur)
if (no.nodes != dim(disttom)[[1]] ) {print("ERROR: number of color labels does not equal number of nodes in disttom")} else {
labeltree=as.character(couleur)
labelrow  = labeltree
labelrow[hier1$order[length(labeltree):1]]=labelrow[hier1$order]
options(expressions = 10000)
if (terrainColors) heatmap(as.matrix(disttom),Rowv=as.dendrogram(hier1),Colv= as.dendrogram(hier1), scale="none",revC=T, ColSideColors=as.character(labeltree),RowSideColors=as.character(labelrow),
labRow=F, labCol=F, col = terrain.colors(1000)) else heatmap(as.matrix(disttom),Rowv=as.dendrogram(hier1),Colv= as.dendrogram(hier1), scale="none",revC=T, ColSideColors=as.character(labeltree),RowSideColors=as.character(labelrow),
labRow=F, labCol=F)
}
} #end of function


# ===================================================
# The function TOMplot2 creates a TOM plot where the top and left color bars can be different
# Inputs:  distance measure, hierarchical (hclust) object, color label=couleurTop, couleurLeft
if (exists("TOMplot2")) rm(TOMplot2);
TOMplot2=function(disttom,hier1, couleurTop, couleurLeft) {
no.nodes=length(couleurTop)
if (no.nodes != length(couleurLeft)) {stop("ERROR: number of top color labels does not equal number of left color labels")}
if (no.nodes != dim(disttom)[[1]] ) {stop("ERROR: number of color labels does not equal number of nodes in disttom")} else {
labeltree = as.character(couleurTop)
labelrow  = as.character(couleurLeft)
labelrow[hier1$order[length(labeltree):1]]=labelrow[hier1$order]
options(expressions = 10000)
heatmap(as.matrix(disttom),Rowv=as.dendrogram(hier1),Colv= as.dendrogram(hier1), scale="none", revC=T, ColSideColors=as.character(labeltree),RowSideColors=as.character(labelrow),
labRow=F, labCol=F)
}
} #end of function



# IGNORE THIS FUNCTION...
# The function "BestHeightCut" allows one to find the best height cut-off
# for a hierarchical clustering tree when external gene information is available
# It computes a Kruskal Wallis-test p-value for each height cut-off
# based on determining whether gene significance differs across branch membership.
if(exists("BestHeightCut")) rm(BestHeightCut);
BestHeightCut=function(hier1, GeneSignif, hcut=seq(0.1,.95,by=0.01) ) {
pvalues=rep(NA, length(hcut))
for (i in c(1:length(hcut))) {
colorhelp=modulecolor2(hier1,hcut[i])
if (length(unique(colorhelp))>1 ) {pvalues[i]=kruskal.test(GeneSignif, colorhelp)$p.value}
data.frame(hcut,pvalues)
}}




#####################################################################################################
################################################################################################################################
# D) Summing up modules using their first principal components (first eigengene)
#####################################################################################################
################################################################################################################################

# ===================================================
#The function ModulePrinComps1 finds the first principal component (eigengene) in each 
# module defined by the colors of the input vector "couleur" (Pardon my French).
# It also reports the variances explained by the first 5 principal components.
# This requires the R library impute
# The output is a list with 2 components: 
# 1) a data frame of module eigengenes (PCs), 
# 2) a data frame that lists the percent variance explained by the first 5 PCs of a module
# Options: if removeGrey=T, then no output is generated for the grey module.
# Recall that grey often denotes genes outside proper modules. 
if (exists("ModulePrinComps1" ) ) rm(ModulePrinComps1);
ModulePrinComps1=function(datexpr,couleur=NULL,removeGrey=F, FiveComponents=T) {
if (is.null(couleur) )  {couleur=rep("turquoise", dim(data.frame(datexpr))[[2]])}
modlevels=levels(factor(couleur))
if ( removeGrey ) modlevels=setdiff(modlevels, c("grey") ); 
PrinComps=data.frame(matrix(666,nrow=dim(datexpr)[[1]],ncol= length(modlevels))) 
varexplained= data.frame(matrix(666,nrow= 5,ncol= length(modlevels)))
names(PrinComps)=paste("PC",modlevels,sep="")
for(i in c(1:length(modlevels)) ){
#print(i)
modulename    = modlevels[i]
restrict1= as.character(couleur)== modulename
# in the following, rows are genes and columns are samples
    datModule=t(datexpr[, restrict1])
    datModule=impute.knn(as.matrix(datModule))
datModule=t(scale(t(datModule)))
if (FiveComponents ) { svd1 =svd(datModule, nu = 5, nv = 5) } else {svd1=svd(datModule)}
    mtitle=paste("PCs of ", modulename," module", sep="")
varexplained[,i]= (svd1$d[1:5])^2/sum(svd1$d^2)
# this is the first principal component
    pc1=svd1$v[,1]
signh1=sign(sum(cor(pc1,  t(datModule))))
if (signh1 != 0)  pc1=signh1* pc1
PrinComps[,i]= pc1
}
list(PrinComps=PrinComps, varexplained=varexplained)
}






#####################################################################################################
################################################################################################################################
# E) Relating a measure of gene significance to the modules 
#####################################################################################################
################################################################################################################################

# ===================================================
# The function ModuleEnrichment1 creates a bar plot that shows whether modules are enriched with
# significant genes.
# More specifically, it reports the mean gene significance for each module.
# The gene significance can be a binary variable or a quantitative variable. 
# It also plots the 95% confidence interval of the mean (CI=mean +/- 1.96* standard error).
# It also reports a Kruskal Wallis P-value.
if( exists("ModuleEnrichment1") ) rm(ModuleEnrichment1);
ModuleEnrichment1=function(genesignif1,couleur,title1="gene significance across modules",labely="Gene Significance",boxplot=F) {
if (length(genesignif1) != length(couleur) ) print("Error: vectors don't have the same lengths") 
else {
if (boxplot != TRUE) {
mean1=function(x) mean(x,na.rm=T) 
means1=as.vector(tapply(genesignif1,couleur,mean1));
se1= as.vector(tapply(genesignif1,couleur,stderr1))
par(mfrow=c(1,1))
barplot(means1, names.arg=names(table(couleur) ),col= names(table(couleur) ) ,ylab=labely)
err.bp(as.vector(means1), as.vector(1.96*se1), two.side=T)}
 else {
boxplot(split(genesignif1,couleur),notch=T,varwidth=T, col= names(table(couleur) ),ylab=labely)}
no.colors=length(names(table(couleur) ))
if (no.colors==1) pp=NA
if (no.colors>1) {
pp=try(kruskal.test(genesignif1,factor(couleur))$p.value)
 if (class(pp)=='try-error') pp=NA
}
title(paste(title1,", p-value=", signif(pp,2) ))
}
} # end of function




# IGNORE THIS...
# ===================================================
#The function fisherPvector allows one to compute Fisher exact p-values
# Thus it allows one to carry out an EASE analysis
# Output: a table of Fisher��s exact p-values
# Input: annotation1 is a vector of gene annotations
# Input: couleur (French for color) denotes the module color of each gene
# Only those gene functions (levels of annotation1) that occur a certain mininum number of times
# (parameter= minNumberAnnotation) in the data set will be considered.  
if (exists("fisherPvector" ) ) rm(fisherPvector);
fisherPvector=function(couleur,annotation1,minNumberAnnotation=50) {
levelsannotation1=levels(annotation1)
levelscouleur=levels(factor(couleur))
no.couleur=length(levelscouleur)
restannotation1=table(annotation1)>minNumberAnnotation
no.annotation=sum( restannotation1)
datoutP=data.frame(matrix(666,nrow=no.annotation,ncol=no.couleur) )
#datoutProp=data.frame(matrix(666,nrow=no.annotation,ncol=2*no.couleur) )
#names(datoutProp)=paste("Prop",paste( rep(levelscouleur ,rep(2, length(levelscouleur))) ) , c("Y","N")  ,sep=".")
datoutProp=data.frame(matrix(666,nrow=no.annotation,ncol=no.couleur) )
names(datoutProp)=paste("Perc",levelscouleur , sep=".")
names(datoutP)=paste("p",levelscouleur,sep=".")
restlevelsannotation1= levelsannotation1[restannotation1]
row.names(datoutP)= restlevelsannotation1
for (i in c(1:no.annotation) ) {
for (j in c(1:no.couleur) ){
tab1=table( annotation1 !=restlevelsannotation1[i], couleur !=levelscouleur[j])
datoutP[i,j]=signif(fisher.test(tab1)$p.value,2) 
#datoutProp[i,2*j-1]=signif(tab1[1,1]/sum(tab1[,1] ),2)
#datoutProp[i,2*j]= signif(tab1[1,2]/sum(tab1[,2]) ,2)
} 
table2=table(annotation1 !=restlevelsannotation1[i], couleur)
datoutProp[i,]= signif(table2[1,]/apply(table2,2,sum),2)
}
data.frame(datoutP,datoutProp)
} # end of function fisherPvector



#####################################################################################################
################################################################################################################################
# F) Carrying out a within module analysis (computing intramodular connectivity etc) 
#####################################################################################################
################################################################################################################################

# ===================================================
#The function DegreeInOut computes for each gene 
#a) the total number of connections, 
#b) the number of connections with genes within its module, 
#c) the number of connections with genes outside its module
# When scaledToOne=TRUE, the within module connectivities are scaled to 1, i.e. the max(K.Within)=1 for each module
if (exists("DegreeInOut")) rm(DegreeInOut); DegreeInOut =function(adj1, couleur,scaledToOne=FALSE) {
no.nodes=length(couleur)
couleurlevels=levels(factor(couleur))
no.levels=length(couleurlevels)
kWithin=rep(-666,no.nodes )
diag(adj1)=0
for (i in c(1:no.levels) ) {
rest1=couleur==couleurlevels[i];
if (sum(rest1) <3 ) { kWithin[rest1]=0 } else {
kWithin[rest1]=apply(adj1[rest1,rest1],2,sum)
if (scaledToOne) kWithin[rest1]=kWithin[rest1]/max(kWithin[rest1])}
}
kTotal= apply(adj1,2,sum) 
kOut=kTotal-kWithin
if (scaledToOne) kOut=NA
kDiff=kWithin-kOut
data.frame(kTotal,kWithin,kOut,kDiff)
}


# =======================================================================
# The function WithinModuleCindex1 relates the node measures (e.g. connectivities) 
# to "external" node significance information within each  module,
# i.e. it  carries out a by module analysis.
# Output: first column reports the spearman correlation p-value between the network variable and the 
# node significance. The next columns contain the Spearman correlations between the variables.
if (exists("WithinModuleAnalysis1")) rm(WithinModuleAnalysis1);
WithinModuleAnalysis1=function(datnetwork,nodesignif, couleur) {
cortesthelp=function( x ) {
len1=dim(x)[[2]]-1
out1=rep(666, len1);
for (i in c(1:len1) ) {out1[i]= signif( cor.test(x[,i+1], x[,1], method="s",use="p" )$p.value ,2) }
data.frame( variable=names(x)[-1] , NS.CorPval=out1, NS.cor=t(signif(cor (x[,1], x[,-1],use="p",method="s"),2)), signif(cor(x[,-1],use="p",method="s"),2) )
} #end of function cortesthelp
print("IGNORE  the warnings...");
by( data.frame(nodesignif, datnetwork),couleur,cortesthelp);
} #end of function WithinModuleAnalysis




# =======================================================================
# The function WithinModuleCindex1 relates the node measures (e.g. connectivities) 
# to "external" node significance information within each  module, 
# i.e. it  carries out a by module analysis.
# BUT it focuses on the C-index also known as area under the ROC curve
# This measure is related to Kendall's Tau statistic and Somer's D, 
# see F. Harrel (Regression Modeling Strategies). Springer. 
# It requires the following library
library(Hmisc)
# Output: the first column reports the C-index and the second, p-value 
if (exists("WithinModuleCindex1")) rm(WithinModuleCindex1);
WithinModuleCindex1=function(datnetwork,nodesignif, couleur) {
CindexFun=function( x ) {
len1=dim(x)[[2]]-1
outC=rep(666, len1);
outP=rep(666, len1);
for (i in c(1:len1) ) {rcor1=rcorr.cens(x[,i+1], x[,1])
outC[i]=rcor1[[1]] 
outP[i]=1- pnorm(abs(rcor1[[2]]/rcor1[[3]]))
}
data.frame( variable=names(x)[-1] , C.index=outC, p.value=outP)
} #end of function CindexFun
#print("IGNORE  the warnings...");
by( data.frame(nodesignif, datnetwork),couleur,CindexFun);
} #end of function WithinModuleAnalysis


# The following function allows on to plot a gene (node) significance measure versus
# connectivity.
if(exists("plotConnectivityGeneSignif1") ) rm( plotConnectivityGeneSignif1);
plotConnectivityGeneSignif1=function(degree1,genesignif1,color1="black", 
title1="Gene Significance vs Connectivity" , xlab1="Connectivity", ylab1="GeneSignificance") {
lm1=lm(genesignif1~degree1 ,na.action="na.omit")
plot(degree1, genesignif1, col=color1,ylab=ylab1,xlab=xlab1,main=paste(title1, ", cor=",  
signif(cor( genesignif1,degree1, method="s",use="p" )   ,2) ))
abline(lm1)
}




 var1=function(x) var(x, na.rm=T)
 
 
if (exists("no.present")) rm(no.present);
no.present=function(x) sum(!is.na(x))
 




# The following function computes the module eigengene based connectivity.
# Input: datExpr= a possibly very large gene expression data set where the rows
# correspond to samples and the columns represent genes
# datPC=data frame of module eigengenes (columns correspond to module eigengenes or PCs)
# A module eigengene based connectivity KME value will be computed if the gene has 
# a non missing expression value in at least MinimumNoSamples arrays.
# Output a data frame where columns are the KME values 
# corresponding to different modules.
# By splitting the expression data into different batches, the function can deal with 
# tens of thousands of gene expression data. 
# If there are many eigengenes (say hundreds) consider decreasing the batch size.
if(exists("signedKME") ) rm(signedKME) 
signedKME=function(datExpr, datPC, outputColumnName="kME") {
datExpr=data.frame(datExpr)
datPC=data.frame(datPC)
output=list()
if (dim(as.matrix(datPC))[[1]] != dim(as.matrix(datExpr))[[1]] ) stop("ERROR in signedKME function: dim(datPC)[[1]] != dim(datExpr)[[1]] ")
varianceZeroIndicatordatExpr=as.vector(apply(as.matrix(datExpr),2,var1))==0
varianceZeroIndicatordatPC =as.vector(apply(as.matrix(datPC),2,var1))==0
if (sum(varianceZeroIndicatordatExpr,na.rm=T)>0 ) warning("warning in signedKME: some genes are constant. Hint: consider removing constant columns from datExpr." )
if (sum(varianceZeroIndicatordatPC,na.rm=T)>0 ) warning("warning in signedKME: some module eigengenes are constant, which is very weird. Module eigengenes should have mean zero and variance 1. Hint: consider removing constant columns from datPC." )
 no.presentdatExpr=as.vector(apply(!is.na(as.matrix(datExpr)),2, sum) )
if (min(no.presentdatExpr)<4 ) warning("warning in signedKME: some gene expressions have fewer than 4 observations. Hint: consider removing genes with too many missing values or collect more arrays.")
output=data.frame(cor(datExpr, datPC, use="p"))
output[no.presentdatExpr<4,]=NA
names(output)=paste(outputColumnName, substring(names(datPC), first=3, last=100), sep="")  
dimnames(output)[[1]] = names(datExpr) 
output
} # end of function signedKME
 








#####################################################################################################
################################################################################################################################
# G) Miscellaneous other functions, e.g. for computing the cluster coefficient.
#####################################################################################################
################################################################################################################################



# ===================================================
# The function ClusterCoef.fun computes the cluster coefficients.
# Input is an adjacency matrix 
if(exists("ClusterCoef.fun")) rm(ClusterCoef.fun) ; ClusterCoef.fun=function(adjmat1) {
diag(adjmat1)=0
no.nodes=dim(adjmat1)[[1]]
computeLinksInNeighbors <- function(x, imatrix){x %*% imatrix %*% x}
nolinksNeighbors <- c(rep(-666,no.nodes))
total.edge <- c(rep(-666,no.nodes))
maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) ); 
if (maxh1>1 | minh1 < 0 ) {print(paste("ERROR: the adjacency matrix contains entries that are larger than 1 or smaller than 0!!!, max=",maxh1,", min=",minh1)) } else { 
nolinksNeighbors <- apply(adjmat1, 1, computeLinksInNeighbors, imatrix=adjmat1)
plainsum  <- apply(adjmat1, 1, sum)
squaresum <- apply(adjmat1^2, 1, sum)
total.edge = plainsum^2 - squaresum
CChelp=rep(-666, no.nodes)
CChelp=ifelse(total.edge==0,0, nolinksNeighbors/total.edge)
CChelp}
} # end of function



# ===================================================
# The function err.bp  is used to create error bars in a barplot
# usage: err.bp(as.vector(means), as.vector(stderrs), two.side=F)
err.bp<-function(daten,error,two.side=F){
 if(!is.numeric(daten)) {
      stop("All arguments must be numeric")}
 if(is.vector(daten)){ 
    xval<-(cumsum(c(0.7,rep(1.2,length(daten)-1)))) 
 }else{
    if (is.matrix(daten)){
      xval<-cumsum(array(c(1,rep(0,dim(daten)[1]-1)),
dim=c(1,length(daten))))+0:(length(daten)-1)+.5
    }else{
      stop("First argument must either be a vector or a matrix") }
 }
 MW<-0.25*(max(xval)/length(xval)) 
 ERR1<-daten+error 
 ERR2<-daten-error
 for(i in 1:length(daten)){
    segments(xval[i],daten[i],xval[i],ERR1[i])
    segments(xval[i]-MW,ERR1[i],xval[i]+MW,ERR1[i])
    if(two.side){
      segments(xval[i],daten[i],xval[i],ERR2[i])
      segments(xval[i]-MW,ERR2[i],xval[i]+MW,ERR2[i])
    } 
 } 
} 

# ===================================================
# this function computes the standard error
if (exists("stderr1")) rm(stderr1)
stderr1 <- function(x){ sqrt( var(x,na.rm=T)/sum(!is.na(x))   ) }



# ===================================================
# The following two functions are for displaying the pair-wise correlation in a panel when using the command "pairs()"
# Typically, we use "pairs(DATA, upper.panel=panel.smooth, lower.panel=panel.cor, diag.panel=panel.hist)" to
# put the correlation coefficients on the lower panel.
panel.cor <- function(x, y, digits=2, prefix="", cex.cor){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex * r)
}
panel.hist <- function(x, ...){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}



# ===================================================
# this function computes the standard error
if (exists("stderr1")) rm(stderr1);
stderr1 <- function(x){ sqrt( var(x,na.rm=T)/sum(!is.na(x))   ) }




# ===================================================
# This function collects garbage
if (exists("collect_garbage")) rm(collect_garbage);
collect_garbage=function(){while (gc()[2,4] != gc()[2,4] | gc()[1,4] != gc()[1,4]){}}
collect_garbage()


# this function is used for computing the Rand index below...
# ===================================================
if (exists("choosenew") ) rm(choosenew)
choosenew <- function(n,k){
  n <- c(n)
  out1 <- rep(0,length(n))
  for (i in c(1:length(n)) ){
    if (n[i]<k) {out1[i] <- 0}
    else {out1[i] <- choose(n[i], k)}}
  out1	
}


# ===================================================
# the following function computes the Rand index between 2 clusterings
# assesses how similar two clusterings are
if (exists("Rand1") ) rm(Rand1)
Rand2 <- function(tab,adjust=T) {
  a <- 0; b <- 0; c <- 0; d <- 0; nn <- 0
  m <- nrow(tab);
  n <- ncol(tab);
  for (i in 1:m) {
    c<-0
    for(j in 1:n) {
      a <- a+choosenew(tab[i,j],2)
      nj <- sum(tab[,j])
      c <- c+choosenew(nj,2)
    }
    ni <- sum(tab[i,])
    b <- b+choosenew(ni,2)
    nn <- nn+ni
  }
  if(adjust==T) {
    d <- choosenew(nn,2)
    adrand <- (a-(b*c)/d)/(0.5*(b+c)-(b*c)/d)
    adrand
  } else {
    b <- b-a
    c <- c-a
    d <- choosenew(nn,2)-a-b-c
    rand <- (a+d)/(a+b+c+d)
    rand
  }
}

# ===================================================
# This function is used in "pairs()" function. The problem of the original  panel.cor is that 
# when the correlation coefficient is very small, the lower panel will have a large font 
# instead of a mini-font in a saved .ps file. This new function uses a format for corr=0.2 
# when corr<0.2, but it still reports the original value of corr, with a minimum format.

panel.cor1=function(x, y, digits=2, prefix="", cex.cor){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    txt1=txt
    r1=r
    if (r<0.2) {
        r1=0.2
        txt1 <- format(c(r1, 0.123456789), digits=digits)[1]
        txt1 <- paste(prefix, txt1, sep="")
        }
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt1)
    cex = cex * r1
    r <- round(r, digits)
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    text(0.5, 0.5, txt, cex=cex)
}

 







# H) THE FOLLOWING FUNCTIONS were originally created by Bin Zhang and later enhanced by Peter Langfelder
# These functions allow one to cut small branches of a hierarchical tree.
# Instead of a fixed (static) height cut-off, the cutreeDynamic adaptively
# finds cut-off heights for each branch separately. 
# See Ghazalpour et al 2006.

if( exists("cutreeDynamic") ) rm(cutreeDynamic)

cutreeDynamic = function(hierclust, maxTreeHeight=1, deepSplit=TRUE, minModuleSize=50, minAttachModuleSize=100, nameallmodules=FALSE, useblackwhite=FALSE)
{
    if(maxTreeHeight >=1){
      staticCutCluster = cutTreeStatic(hiercluster=hierclust, heightcutoff=0.99, minsize1=minModuleSize)
    }else{
      staticCutCluster = cutTreeStatic(hiercluster=hierclust, heightcutoff=maxTreeHeight, minsize1=minModuleSize)
    }      

    #get tree height for every singleton
    #node_index   tree_height
    demdroHeiAll= rbind( cbind(hierclust$merge[,1], hierclust$height), cbind(hierclust$merge[,2], hierclust$height) )

    #singletons will stand at the front of the list
    myorder = order(demdroHeiAll[,1])

    #get # of singletons
    no.singletons = length(hierclust$order)
    demdroHeiAll.sort = demdroHeiAll[myorder, ]
    demdroHei.sort = demdroHeiAll.sort[c(1:no.singletons), ]
    demdroHei = demdroHei.sort[seq(no.singletons, 1, by=-1), ]
    demdroHei[,1]     = -demdroHei[,1]

    # combine with prelimilary cluster-cutoff results
    demdroHei         = cbind(demdroHei, as.integer(staticCutCluster))

    # re-order the order based on the dendrogram order hierclust$order
    demdroHei.order = demdroHei[hierclust$order, ]

    static.clupos = locateCluster(demdroHei.order[, 3])

    if (is.null(static.clupos) ){
        module.assign   = rep(0, no.singletons)
        colcode.reduced = assignModuleColor(module.assign, minsize1=minModuleSize,  anameallmodules=nameallmodules, auseblackwhite=useblackwhite)
        return ( colcode.reduced )
    }

    static.no     = dim(static.clupos)[1]
    static.clupos2 =     static.clupos
    static.no2     =     static.no

    #split individual cluster if there are sub clusters embedded
    mcycle=1
    while(1==1){
        clupos = NULL
        for (i in c(1:static.no)){
           mydemdroHei.order = demdroHei.order[ c(static.clupos[i,1]:static.clupos[i,2]), ] #index to [1, clusterSize]
           mydemdroHei.order[, 1] = mydemdroHei.order[, 1] - static.clupos[i, 1] + 1

           #cat("Cycle ", as.character(mcycle), "cluster (", static.clupos[i,1], static.clupos[i,2], ")\n")
           #cat("i=", as.character(i), "\n")

           iclupos = processInvididualCluster(mydemdroHei.order, 
                                              cminModuleSize       = minModuleSize, 
                                              cminAttachModuleSize = minAttachModuleSize)

           iclupos[,1] = iclupos[,1] + static.clupos[i, 1] -1 #recover the original index
           iclupos[,2] = iclupos[,2] + static.clupos[i, 1] -1

           clupos  = rbind(clupos, iclupos) #put in the final output buffer
        }

        if(deepSplit==FALSE){
          break
        }

        if(dim(clupos)[1] != static.no) {
           static.clupos = clupos
           static.no     = dim(static.clupos)[1]
        }else{
          break
        }
       mcycle = mcycle + 1
       #static.clupos  
    }
    final.cnt = dim(clupos)[1]
    #assign colors for modules
    module.assign = rep(0, no.singletons)
    module.cnt=1
    for (i in c(1:final.cnt ))
    {
       sdx = clupos[i, 1] #module start point
       edx = clupos[i, 2] #module end point
       module.size = edx - sdx +1 
       if(module.size <minModuleSize){
         next
       }
       #assign module lable
       module.assign[sdx:edx] = rep(module.cnt, module.size)
       #update module label for the next module
       module.cnt = module.cnt + 1
    }
    colcode.reduced.order = assignModuleColor(module.assign, minsize1=minModuleSize,  anameallmodules=nameallmodules, auseblackwhite=useblackwhite)
   recov.order = order( demdroHei.order[,1])
    colcode.reduced = colcode.reduced.order[recov.order]
    if(1==2){
        aveheight = averageSequence(demdroHei.order[,2], 2)
        procHei = demdroHei.order[,2]-mean(demdroHei.order[,2])
        par(mfrow=c(3,1), mar=c(0,0,0,0) )
        plot(h1row, labels=F, xlab="",ylab="",main="",sub="",axes = F)
           barplot(procHei, 
                col= "black", space=0,
                border=F,main="", axes = F, axisnames = F)
        barplot(height=rep(1, length(colcode.reduced)), 
                col= as.character(colcode.reduced[h1row$order]), space=0,
                border=F,main="", axes = F, axisnames = F)
        par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
    }

    colcode.reduced
}




#use height cutoff to remove
cutTreeStatic = function(hiercluster,heightcutoff=0.99, minsize1=50) {

    # here we define modules by using a height cut-off for the branches
    labelpred= cutree(hiercluster,h=heightcutoff)
    sort1=-sort(-table(labelpred))
    sort1
    modulename= as.numeric(names(sort1))
    modulebranch= sort1 >= minsize1
    no.modules=sum(modulebranch)

    colorhelp = rep(-1, length(labelpred) )
    if ( no.modules==0){
        print("No module detected")
    }
    else{
        for (i in c(1:no.modules)) {
            colorhelp=ifelse(labelpred==modulename[i],i ,colorhelp)
        }
    }
    colorhelp
}

#leftOrright >0 : running length (with same sign) to right, otherwise to the left
#mysign = -1: negative value, mysign = -1: positive value
runlengthSign = function(mysequence, leftOrright=-1, mysign=-1){
    seqlen = length(mysequence)   
    if(leftOrright<0){
        pseq = rev(mysequence)
    }else{
        pseq = mysequence
    }

    if(mysign<0){ #see where the first POSITIVE number occurs
       nonezero.bool = (pseq > 0)
    }else{ #see where the first NEGATIVE number occur
       nonezero.bool = (pseq < 0)
    }
    if( sum(nonezero.bool) > 0){
      runlength = min( c(1:seqlen)[nonezero.bool] ) - 1
    }else{
      runlength = 0
    }
}




#"0" is for grey module
assignModuleColor = function(labelpred, minsize1=50, anameallmodules=FALSE, auseblackwhite=FALSE) {
    # here we define modules by using a height cut-off for the branches
    #labelpred= cutree(hiercluster,h=heightcutoff)
    #cat(labelpred)

    #"0", grey module doesn't participate color assignment, directly assigned as "grey"
    labelpredNoZero = labelpred[ labelpred >0 ]
    sort1=-sort(-table(labelpredNoZero))
    sort1
    modulename= as.numeric(names(sort1))
    modulebranch= sort1 >= minsize1
    no.modules=sum(modulebranch)

    # now we assume that there are fewer than 10 modules
    colorcode=GlobalStandardColors;

    #"grey" means not in any module;
    colorhelp=rep("grey",length(labelpred))
    if ( no.modules==0){
        print("No module detected\n")
    }
    else{
        if ( no.modules > length(colorcode)  ){
            print( paste("Too many modules (more than colors):", as.character(no.modules)) )
        }

        if ( (anameallmodules==FALSE) || (no.modules <=length(colorcode)) ){
            labeledModules = min(no.modules, length(colorcode) )
            for (i in c(1:labeledModules)) {
               colorhelp=ifelse(labelpred==modulename[i],colorcode[i],colorhelp)
            }
            colorhelp=factor(colorhelp,levels=c(colorcode[1:labeledModules],"grey"))
        }else{#nameallmodules==TRUE and no.modules >length(colorcode)
            maxcolors=length(colorcode)
            labeledModules = no.modules
            extracolors=NULL
            blackwhite=c("red", "black")
            for(i in c((maxcolors+1):no.modules)){
              if(auseblackwhite==FALSE){
                  icolor=paste("module", as.character(i), sep="")
              }else{#use balck white alternatively represent extra colors, for display only
                  #here we use the ordered label to avoid put the same color for two neighboring clusters
                  icolor=blackwhite[1+(as.integer(modulename[i])%%2) ]
              }
              extracolors=c(extracolors, icolor)
            }

            #combine the true-color code and the extra colorcode into a uniform colorcode for 
            #color assignment
            allcolorcode=c(colorcode, extracolors)

            for (i in c(1:labeledModules)) {
               colorhelp=ifelse(labelpred==modulename[i],allcolorcode[i],colorhelp)
            }
            colorhelp=factor(colorhelp,levels=c(allcolorcode[1:labeledModules],"grey"))
        }
    }


    colorhelp
}





#locate the start/end positions of each cluster in the ordered cluster label sequence 
#where "-1" indicating no cluster
#3-1 -1 1 1 1 1 2 2 2
#3 3 -1-1 1 1 1 1 2 2 2  (shift)
#---------------------------------
#0-4  0 2 0 0 0 1 0 0 0   (difference)
#       *     * @
locateCluster = function(clusterlabels)
{
 no.nodes = length(clusterlabels)
 clusterlabels.shift = c(clusterlabels[1], c(clusterlabels[1:(no.nodes-1)]) )
 
 #a non-zero point is the start point of a cluster and it previous point is the end point of the previous cluster
 label.diff = abs(clusterlabels - clusterlabels.shift)

 #process the first and last positions as start/end points if they belong to a cluster instead of no cluster "-1"
 if(clusterlabels[1]       >0) {label.diff[1]=1} 
 if(clusterlabels[no.nodes]>0) {label.diff[no.nodes]=1} 

 flagpoints.bool = label.diff > 0
 if( sum(flagpoints.bool) ==0){
     return(NULL)
 }

 flagpoints = c(1:no.nodes)[flagpoints.bool]
 no.points  = length(flagpoints)

 myclupos=NULL
 for(i in c(1:(no.points-1)) ){
   idx = flagpoints[i]
   if(clusterlabels[idx]>0){
      myclupos = rbind(myclupos, c(idx, flagpoints[i+1]-1) )
   }
 }
 myclupos
}


#input is the cluster demdrogram of an individual cluster, we want to find its embbeded subclusters
#execution order: mean-height ==> (mean+max)/2 ==> (mean+min)/2
#useMean: =0 ~ use mean-height   as calibation line
#         =1 ~ use (mean+max)/2  as calibation line to detect relatively a small cluster sitting on the head of a bigger one,
#                      so mean-height is too low to detect the two modules.
#         =-1~ use (mean+min)/2  as calibation line to detect relatively a small cluster sitting on the tail of a bigger one,
#                      so mean-height & (mean+max)/2 are too high to detect the two modules

processInvididualCluster = function(clusterDemdroHei, cminModuleSize=50, cminAttachModuleSize=100, minTailRunlength=12, useMean=0){
    #for debug: use all genes
    #clusterDemdroHei =demdroHei.order

    no.cnodes = dim(clusterDemdroHei)[1]
    
    cmaxhei   = max(clusterDemdroHei[, 2])
    cminhei   = min(clusterDemdroHei[, 2])
    
    cmeanhei  = mean(clusterDemdroHei[, 2])
    cmidhei = (cmeanhei + cmaxhei)/2.0
    cdwnhei = (cmeanhei + cminhei)/2.0

    if (useMean==1){
        comphei = cmidhei
    }else if (useMean==-1){
        comphei = cdwnhei
    }else{ #normal case
        comphei = cmeanhei
    }
        
    # compute height diffrence with mean height
    heidiff       = clusterDemdroHei[,2] - comphei
    heidiff.shift = shiftSequence(heidiff, -1)

    # get cut positions
    # detect the end point of a cluster, whose height should be less than meanhei 
    #  and the node behind it is the start point of the next cluster which has a height above meanhei
    cuts.bool = (heidiff<0) & (heidiff.shift > 0)
    cuts.bool[1]         = TRUE
    cuts.bool[no.cnodes] = TRUE
    
    if(sum(cuts.bool)==2){
          if (useMean==0){
             new.clupos=processInvididualCluster(clusterDemdroHei=clusterDemdroHei, cminModuleSize=cminModuleSize, 
                                         cminAttachModuleSize=cminAttachModuleSize,
                                         useMean=1)
          }else if(useMean==1){
             new.clupos=processInvididualCluster(clusterDemdroHei=clusterDemdroHei, cminModuleSize=cminModuleSize, 
                                         cminAttachModuleSize=cminAttachModuleSize,
                                         useMean=-1)          
          }else{
             new.clupos = rbind(c(1, no.cnodes))
          }
          return (new.clupos)
    }

    #a good candidate cluster-end point should have significant # of ahead nodes with head < meanHei
    cutindex =c(1:no.cnodes)[cuts.bool]
    no.cutps = length(cutindex)
    runlens  = rep(999, no.cutps)
    cuts.bool2 = cuts.bool
    for(i in c(2:(no.cutps-1)) ){
       seq = c( (cutindex[i-1]+1):cutindex[i] )
       runlens[i] = runlengthSign(heidiff[seq], leftOrright=-1, mysign=-1)

       if(runlens[i] < minTailRunlength){
          #cat("run length=", runlens[i], "\n")
          cuts.bool2[ cutindex[i] ] = FALSE
       }
    }

    #attach SMALL cluster to the left-side BIG cluster if the small one has smaller mean height
    cuts.bool3=cuts.bool2
    if(sum(cuts.bool2) > 3) {
       curj = 2
       while (1==1){
           cutindex2 =c(1:no.cnodes)[cuts.bool2]
           no.clus = length(cutindex2) -1
           if (curj>no.clus){
              break
           }
           pre.sdx = cutindex2[ curj-1 ]+1 #previous module start point
           pre.edx = cutindex2[ curj ] #previous module end   point
           pre.module.size = pre.edx - pre.sdx +1 
           pre.module.hei  = mean(clusterDemdroHei[c(pre.sdx:pre.edx) , 2])
         
           cur.sdx = cutindex2[ curj ]+1 #previous module start point
           cur.edx = cutindex2[ curj+1 ] #previous module end   point
           cur.module.size = cur.edx - cur.sdx +1 
           cur.module.hei  = mean(clusterDemdroHei[c(cur.sdx:cur.edx) , 2])

           #merge to the leftside major module, don't change the current index "curj"
           #if( (pre.module.size >minAttachModuleSize)&(cur.module.hei<pre.module.hei)&(cur.module.size<minAttachModuleSize) ){
           if( (cur.module.hei<pre.module.hei)&(cur.module.size<cminAttachModuleSize) ){
                cuts.bool2[ cutindex2[curj] ] = FALSE
           }else{ #consider next cluster
                curj = curj + 1
           }
       }#while
   }#if 

   cutindex2 =c(1:no.cnodes)[cuts.bool2]
   no.cutps = length(cutindex2)

    #we don't want to lose the small cluster at the tail, attch it to the previous big cluster
    #cat("Lclu= ", cutindex2[no.cutps]-cutindex2[no.cutps-1]+1, "\n")
    if(no.cutps > 2){
      if( (cutindex2[no.cutps] - cutindex2[no.cutps-1]+1) < cminModuleSize ){
        cuts.bool2[ cutindex2[no.cutps-1] ] =FALSE  
      }
    }

    if(1==2){
        myseqnce = c(2300:3000)
        cutdisp = ifelse(cuts.bool2==T, "red","grey" )
        #re-order to the normal one with sequential signleton index
        par(mfrow=c(3,1), mar=c(0,0,0,0) )
        plot(h1row, labels=F, xlab="",ylab="",main="",sub="",axes = F)
        barplot(heidiff[myseqnce],
                col= "black", space=0,
                border=F,main="", axes = F, axisnames = F)
        barplot(height=rep(1, length(cutdisp[myseqnce])), 
                col= as.character(cutdisp[myseqnce]), space=0,
                border=F,main="", axes = F, axisnames = F)
        par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
    }

   cutindex2  = c(1:no.cnodes)[cuts.bool2]
   cutindex2[1]=cutindex2[1]-1 #the first 
   no.cutps2  = length(cutindex2)

   if(no.cutps2 > 2){
     new.clupos = cbind( cutindex2[c(1:(no.cutps2-1))]+1, cutindex2[c(2:no.cutps2)] )
   }else{
     new.clupos = cbind( 1, no.cnodes)
   }

   if ( dim(new.clupos)[1] == 1 ){   
          if (useMean==0){
             new.clupos=processInvididualCluster(clusterDemdroHei=clusterDemdroHei, cminModuleSize=cminModuleSize, 
                                         cminAttachModuleSize=cminAttachModuleSize,
                                         useMean=1)
          }else if(useMean==1){
             new.clupos=processInvididualCluster(clusterDemdroHei=clusterDemdroHei, cminModuleSize=cminModuleSize, 
                                         cminAttachModuleSize=cminAttachModuleSize,
                                         useMean=-1)          
          }   
   }
   new.clupos
}






findClustersSignificant=function(mysequence, modulecolor)
{
  
   modnames= names( table(modulecolor) )
   mysize     = length(modulecolor)
   validseq     = rep(TRUE, mysize)
   for (each in modnames ){
      mybool = (modulecolor==each)
      mymodulesig = mysequence[mybool]
      mydiff      = abs(mymodulesig - mymodulesig[1])
      if(sum(mydiff)==0){
         validseq = ifelse(mybool==TRUE, FALSE, validseq)
      }
   }
  validseq
}


#leftOrright >0 : running length (with same sign) to right, otherwise to the left
#mysign = -1: negative value, mysign = -1: positive value
runlengthSign = function(mysequence, leftOrright=-1, mysign=-1){
    seqlen = length(mysequence)   
    if(leftOrright<0){
        pseq = rev(mysequence)
    }else{
        pseq = mysequence
    }

    if(mysign<0){ #see where the first POSITIVE number occurs
       nonezero.bool = (pseq > 0)
    }else{ #see where the first NEGATIVE number occur
       nonezero.bool = (pseq < 0)
    }
    if( sum(nonezero.bool) > 0){
      runlength = min( c(1:seqlen)[nonezero.bool] ) - 1
    }else{
      runlength = 0
    }
}


#delta >0 : shift to right, otherwise to the left
shiftSequence = function(mysequence, delta){
    seqlen = length(mysequence)
    if(delta>0){
        finalseq=c(mysequence[1:delta], mysequence[1:(seqlen-delta)])
    }else{
        posdelta = -delta 
        finalseq=c(mysequence[(posdelta+1):seqlen], mysequence[(seqlen-posdelta+1):seqlen])
    }
    finalseq
}


#no of neighbors behind and before the point used for average
averageSequence=function(mysequence, noneighbors){
   sumseq = mysequence
   for(i in c(1:noneighbors)){
        iseq = shiftSequence(mysequence, i)
        sumseq =    sumseq + iseq
        iseq = shiftSequence(mysequence, -i)
        sumseq =    sumseq + iseq
   }
   sumseq = sumseq/(1+2*noneighbors)
}

#delta >0 : shift to right, otherwise to the left
shiftSequence = function(mysequence, delta){
    seqlen = length(mysequence)
    if(delta>0){
        finalseq=c(mysequence[1:delta], mysequence[1:(seqlen-delta)])
    }else{
        posdelta = -delta 
        finalseq=c(mysequence[(posdelta+1):seqlen], mysequence[(seqlen-posdelta+1):seqlen])
    }
    finalseq
}

#find the middle of each cluster and label the middle position with the corrsponding color
getDisplayColorSequence=function(colordered){
 mylen = length(colordered)
 colordered2 = c(colordered[1], colordered[1:(mylen-1)] )
 colordiff   = (colordered != colordered2)
 colordiff[1] = TRUE
 colordiff[mylen] = TRUE
 #mydispcolor = ifelse(colordiff==TRUE, colordered, "")
 mydispcolor = rep("", mylen)
 mytrueseq = c(1:mylen)[colordiff]
 for (i in c(1:(length(mytrueseq)-1)) ){
    midi =  (mytrueseq[i] + mytrueseq[i+1])/2
    mydispcolor[midi] = colordered[midi]
 }
 fdispcolor = ifelse(mydispcolor=="grey", "", mydispcolor)
 fdispcolor
}

#use height cutoff to remove
cutTreeStatic = function(hiercluster,heightcutoff=0.99, minsize1=50) {

    # here we define modules by using a height cut-off for the branches
    labelpred= cutree(hiercluster,h=heightcutoff)
    sort1=-sort(-table(labelpred))
    sort1
    modulename= as.numeric(names(sort1))
    modulebranch= sort1 >= minsize1
    no.modules=sum(modulebranch)

    colorhelp = rep(-1, length(labelpred) )
    if ( no.modules==0){
        print("No module detected")
    }
    else{
        for (i in c(1:no.modules)) {
            colorhelp=ifelse(labelpred==modulename[i],i ,colorhelp)
        }
    }
    colorhelp
}




#merge the minor cluster into the major cluster
merge2Clusters = function(mycolorcode, mainclusterColor, minorclusterColor){
  mycolorcode2 = ifelse(as.character(mycolorcode)==minorclusterColor, mainclusterColor, as.character(mycolorcode) )
  fcolorcode   =factor(mycolorcode2)
  fcolorcode
}



# A set of global variables and functions that should help handling color names for some 400+ modules.
# A vector called GlobalStandardColors is defined that holds color names with first few entries 
# being the well-known and -loved colors. The rest is randomly chosen from the color names of R,
# excluding grey colors.
#---------------------------------------------------------------------------------------------------------
## GlobalStandardColors 
#
#---------------------------------------------------------------------------------------------------------
# This code forms a vector of color names in which the first entries are given by BaseColors and the rest
# is "randomly" chosen from the rest of R color names that do not contain "grey" nor "gray".

BaseColors = c("turquoise","blue","brown","yellow","green","red","black","pink","magenta",
"purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan","grey60", "lightgreen", 
"lightyellow", "royalblue", "darkred", "darkgreen", "darkturquoise", "darkgrey", "orange", "darkorange",
 "white", "skyblue", "saddlebrown", "steelblue", "paleturquoise", "violet", "darkolivegreen",
 "darkmagenta", "white" );
RColors = colors()[-grep("grey", colors())];
RColors = RColors[-grep("gray", RColors)];
InBase = match(BaseColors, RColors);
ExtraColors = RColors[-c(InBase[!is.na(InBase)])];
No.Extras = length(ExtraColors);
# Here is the vector of colors that should be used by all functions:
GlobalStandardColors = c(BaseColors, ExtraColors[rank(sin(13*c(1:No.Extras) +sin(13*c(1:No.Extras))) )] );
rm(BaseColors, RColors, ExtraColors, No.Extras);
#---------------------------------------------------------------------------------------------------------
#
# ColorsFromLabels
#
#---------------------------------------------------------------------------------------------------------
# This function converts integer numerical labels Labels into color names in the order either given by
# StandardColors,
# or (if StandardColors==NULL) by GlobalStandardColors. If GreyIsZero == TRUE, labels 0 will be assigned
# the color grey; otherwise presence of labels below 1 will trigger an error.
ColorsFromLabels = function(Labels, ZeroIsGrey = TRUE, StandardColors = NULL)
{
  if (is.null(StandardColors)) StandardColors = GlobalStandardColors;
  if (ZeroIsGrey) MinLabel = 0 else MinLabel = 1
  if (sum( (Labels>=MinLabel) & (Labels <= length(StandardColors)) )!= length(Labels))
    stop(paste("Input error: something's wrong with Labels. Either they are not a numeric vector,",
               "or some values are below", MinLabel,  "or some values are above the maximum number of colors", length(StandardColors)));
  Colors = rep("grey", length(Labels));
  Colors[Labels!=0] = StandardColors[Labels[Labels!=0]];
  Colors;}
#---------------------------------------------------------------------------------------------------------
#
# DisplayColors
#
#---------------------------------------------------------------------------------------------------------
# This function plots a barplot with all colors given. If Colors are not given, GlobalStandardColors are
# used, i.e. if you want to see the GlobalStandardColors, just call this function without parameters.
DisplayColors = function(Colors = NULL)
{
  if (is.null(Colors)) Colors = GlobalStandardColors;
  barplot(rep(1, length(Colors)), col = Colors, border = Colors);}


#----------------------------------------------------------## ColorsFromLabels##------------------------------



###############################################################################
# I) Functions for merging modules based on a high correlation of the module eigengenes
###############################################################################
# The following network functions were created by Peter Langfelder.


# This function is an alternative to the print command that enforces immediate output.

if (exists("print.flush")) { remove(print.flush); collect_garbage(); }
print.flush = function(...)
{
  x = print(...)
  if (exists("flush.console")) x=flush.console();
}

if (exists("PrintSpaces")) { remove(PrintSpaces); collect_garbage(); }
PrintSpaces = function(print.level)
{
  if (print.level>0) 
    {
      spaces = paste(rep("  ", times=print.level), collapse="");
    } else
    {
      spaces = "";
    }
  spaces;
}




# This function forms the average gene expression for each color in the color vector.
# Input: gene expression data (rows are samples, columns are genes)
# In many applications, we recommend to input normalized gene expression values, i.e.
# NormExprData=scale(datExpr). colors= is a vector that encodes the module color of each gene.

AverageExprMatrix = function(NormExprData, colors) {
  no.genes = dim(NormExprData)[2]
  no.samples = dim(NormExprData)[1]
  colorsf = as.factor(colors)
  AverageExpr = matrix(ncol=nlevels(colorsf), nrow = no.samples)
  ExprDataMtrx = as.matrix(NormExprData)
  for (i in (1:no.samples)) AverageExpr[i,] = tapply(ExprDataMtrx[i,], colorsf, mean)
  AverageExpr
}

#---------------------------------------------------------------------------------
#
# ModuleNumber. Function from Peter Langfelder
#
#---------------------------------------------------------------------------------
# Similar to modulecolor2 above, but returns numbers instead of colors, which is oftentimes more useful.
# 0 means unassigned.
# This is particularly useful when dealing with a lot of branches in the tree.
# Or when integer values instead of strings are preferred to label modules.
# Return value is a simple vector, not a factor.
# Caution: the module numbers are neither sorted nor sequential, the only guarranteed fact is that grey
# probes are labeled by 0 and all probes belonging to the same module have the same number.


ModuleNumber = function(HierTree, CutHeight = 0.9, MinSize = 50)
{
  Branches = cutree(HierTree, h = CutHeight);
  NOnBranches = table(Branches);
  #NOnBranches[i][[1]] somehow gives the number of elements on branch i.
  TrueBranch = NOnBranches >= MinSize;
  Branches[!TrueBranch[Branches]] = 0;
  
  #NewLabels = levels(factor(Branches));
  #for (lab in 1:length(NewLabels)) if (NewLabels[lab]!=0)
  #  Branches[Branches==NewLabels[lab]] = lab;

  Branches;

}

#-------------------------------------------------------------------------------------
#
#  ModulePrincipalComponents
#
#-------------------------------------------------------------------------------------
# This function is an alternative to ModulePrincomps1. Probably it produces the same output.
# This function computes the principal components of modules of a given network.
# The first PC of each module is the module eigengene.
# Input: Data: expression data, module colors. 
# In general, the first PC is only defined up to a sign.
# This function optionally orients the first PC so that it has a positive correlation
# with the average normalized gene expression profile, see the function AverageExpression
# AlignPCs can take the values "", "along average".
# output : a dataframe of module eigengenes (PCs)

ModulePrincipalComponents = function(Data, ModuleColors, AlignPCs = "along average", verbose = 1,
  print.level=0) 
{
  spaces = PrintSpaces(print.level);

  AlignPCsRecognizedValues =  c("", "along average");
  if (!is.element(AlignPCs, AlignPCsRecognizedValues)) {
    print.flush(paste("ModulePrincipalComponents: Error:",
                "parameter AlignPCs has an unrecognised value:", 
                AlignPCs, "; Recognized values are ", AlignPCsRecognizedValues));
    stop()
  }

  if (verbose>0) print.flush(paste(spaces, "ModulePrincipalComponents: Calculating PCs"));
  FullPCs = ModulePrinComps1(Data, ModuleColors);
  PCs = FullPCs$PrinComps;

  if (AlignPCs == "") AlignedPCs = PCs;
  if (AlignPCs == "along average") 
  {
    if (verbose>0) print.flush(paste(spaces,"ModulePrincipalComponents:", 
                     "Aligning PCs with average expression for each module."))
    if (verbose>1) print.flush(paste(spaces,"  ++ Calculating averages..."));
    NormData = scale(Data);
    AverageModuleExpr = data.frame(AverageExprMatrix(NormData, ModuleColors));
    if (verbose>1) print.flush(paste(spaces,"  ++ Aligning principal components..."));
    AverageAndPCCor = diag(cor(PCs, AverageModuleExpr, use = "p"));
    sign.matrix = matrix(0, nrow = dim(PCs)[2], ncol = dim(PCs)[2]);
    diag(sign.matrix) = sign(AverageAndPCCor);
    AlignedPCs = as.data.frame(as.matrix(PCs) %*% sign.matrix);
    names(AlignedPCs) = names(PCs);
    rownames(AlignedPCs) = rownames(PCs);
    names(AverageModuleExpr) = names(PCs);
    rownames(AverageModuleExpr) = rownames(PCs);
    if (verbose>1) print.flush(paste(spaces,"  ++ done."));
  }
  RetPCs = list(data = AlignedPCs, VarExplained = FullPCs$varexplained, 
                ModuleConformity = FullPCs$ModuleConformity, AverageExpr = AverageModuleExpr);
  RetPCs;
}


#--------------------------------------------------------------------------------------
# 
# ModulePCs
# This helper function is used for the function MergeCloseModules.
#--------------------------------------------------------------------------------------
# Input: Data are gene expression data where columns correspond to genes 
# and rows correspond to samples that may correspond to different subsets (e.g. study 1 and study 2)
# Subsets is a vector of labels "1", "2", "3", etc (i.e. positive integers).
# Output: a vector of lists where each list contains the following
# data=module eigengenes, 
# AverageExpr = average normalized module expression
# VarExplained= variance explained by the eigengenes 
# Options: AlignPCs= chooses the sign of each module eigengene by
# enforcing a positive correlation with the average normalized expression.


ModulePCs = function(Data, Subsets, ModuleColors, OnlySet = NULL,
                     AlignPCs="along average", verbose=1, print.level=0)
{
  spaces = PrintSpaces(print.level)
  SubsetFactor = factor(Subsets);
  No.Subsets = nlevels(SubsetFactor);
  if (verbose>0) print.flush(paste(spaces,"ModulePCs: Looking for module PCs."));
  PCs = vector(mode="list", length=No.Subsets);
  if (is.null(OnlySet))
  {
      CalculatedSubsets = c(1:No.Subsets);
  } else {
      CalculatedSubsets = c(OnlySet);
  }
  for (subs.ind in CalculatedSubsets) {
    subset = levels(SubsetFactor)[subs.ind];
    if (verbose>0) print.flush(paste(spaces,"  Working on subset", as.character(subset), "...")); 
    SubsetNetworkData = Data[SubsetFactor == subset, ];
    SubsetColors = ModuleColors; 
    SubsetPCs = ModulePrincipalComponents(Data = SubsetNetworkData, 
                          ModuleColors = SubsetColors, AlignPCs = AlignPCs, verbose = verbose-1,
                          print.level = print.level+1);
    PCs[[subs.ind]] = list(data = SubsetPCs$data, AverageExpr = SubsetPCs$AverageExpr, 
                           ModuleConformity = SubsetPCs$ModuleConformity, 
                           VarExplained = SubsetPCs$VarExplained);
    rm(SubsetNetworkData); rm(SubsetColors); rm(SubsetPCs); collect_garbage();
  }
  PCs;
}


#---------------------------------------------------------------------------------------------
#
# MergeCloseModules
#
#---------------------------------------------------------------------------------------------
# This function merges modules whose PCs fall on one branch of a hierarchical clustering tree
# Method: 
# First, the module eigenes are computed corresponding to each color (but see option PCs)
# Second, we use average linkage hierarchical clustering of the module eigengenes to arrive at a dendrogram
# Third, branches are cut-off the dendrogram using a given height (CutHeight)
# Fourth, modules whose PCs fall on one branch are merged.
# If the option Relabel=True, then the colors of the resulting merged modules 
# are assigned according to size, e.g. turquoise encodes the largest module
# If Relabel=F, then the color of the merged module is chosen according to that of the first module.
# Options: the parameter PCs can be used to input a vector of lists 
# where each list must have a component named data that contains the data frame of module eigengenes.
# UseAbs specifies whether the absolute value of the correlation should be used for constructing the clustering tree of the PCs.
# if CalculateNewPCs = TRUE is specified then the PCs are newly calculated based on the new merged colors.
# The options Subsets=NULL, OnlySet = NULL are useful when the data are comprised of multiple subsets (e.g. study 1 and 2).

# StandardColors allows one to specify your own order of colors to be assigned according to module size.
# Thus, specifying StandardColors only makes sense if Relabel=True.


if(exists("MergeCloseModules")) rm(MergeCloseModules);
MergeCloseModules = function(datE, OriginalColors, CutHeight=0.20, 
                             PCs = NULL, UseAbs = F, 
                             Relabel = FALSE, StandardColors = NULL, CalculateNewPCs = TRUE,
                             Subsets = NULL, OnlySet = NULL, 
                             verbose = 1, print.level=0)
{

  PCsInSingleFrame = FALSE;
  spaces = PrintSpaces(print.level);

  if (verbose>0) print.flush(paste(spaces, 
            "MergeCloseModules: Merging modules whose distance is less than", CutHeight));

  if (is.null(Subsets)) 
  {
    Subsets = rep(1, dim(datE)[1]);
    if (!is.null(PCs))
    {
      if (is.null(PCs[[1]]$data))
      {
        if (verbose>1) print.flush(paste(spaces, "  PCs appear to be a single data frame - will work on",
                                         "that assumption."));
        if (length(dim(PCs))!=2)
          stop(paste("PCs must be given either as a vector of lists, one list for each subset, with",
                     "element 'data' containing the frame or matrix of eigengenes",
                     " OR for a single data set PCs can be given as a dataframe or matrix of correct",
                     "dimensions."));
        PCsX = vector(mode="list", length = 1);
        PCsX[[1]] = list(data = PCs);
        rm(PCs); PCs = PCsX; rm(PCsX);
        PCsInSingleFrame = TRUE;
      } 
      if (dim(PCs[[1]]$data)[1]!=dim(datE)[1])
        stop(paste("Number of samples in PCs is incompatible with number of samples in datE."));
    }
  } else
  {
    SubsetsX = Subsets;
    Subsets = as.numeric(factor(Subsets));
  }

  if (!is.null(PCs))
  {
    SubsetsX = Subsets;
    No.Sets = max(Subsets);
    if (length(PCs)!=No.Sets)
      stop(paste("Number of sets given by Subsets is incompatible with the number of eigengene sets",
                 "given in PCs."));
    SamplesInSet = table(Subsets);
    for (i in 1:No.Sets)
    {
      if (dim(PCs[[i]]$data)[1]!=SamplesInSet[i])
          stop(paste("Number of samples in PCs is incompatible with subset length for subset",
                     levels(factor(SubsetsX))[i], " (first occurence)."));
    }
  }

  if (dim(datE)[1]!=length(Subsets))
    stop("Number of genes in datE is different from the length of given subset vector. They must equal.");

  if (dim(datE)[2]!=length(OriginalColors))
    stop("Number of genes in datE is different from the length of original colors. They must equal.");

  if ((CutHeight <0) | (CutHeight>(1+as.integer(UseAbs)))) 
    stop(paste("Given CutHeight is out of sensible range between 0 and", 1+as.integer(UseAbs) ));

  # If ordered PCs were not given, calculate them

  if (is.null(PCs)) 
  {
    PCs = ModulePCs(datE, Subsets, ModuleColors = OriginalColors,
                    OnlySet = OnlySet, 
                    verbose = verbose-1, print.level = print.level+1);
    collect_garbage();
  } else if (nlevels(as.factor(OriginalColors))!=dim(PCs[[1]]$data)[2])
  {
    if (verbose>0) print.flush(paste(spaces, "  Number of given module colors", 
              "does not match number of given MEs => recalculating the MEs."))
    PCs = ModulePCs(datE, Subsets, ModuleColors = OriginalColors,
                    OnlySet = OnlySet, 
                    verbose = verbose-1, print.level = print.level+1);
    collect_garbage();
  }

  # Cluster the found module eigengenes and merge ones that are too close to one another _in both sets_.

  No.Sets = nlevels(as.factor(Subsets));
  
  MEDiss = vector(mode="list", length = No.Sets);
  if (is.null(OnlySet))
  {
      CalculatedSubsets = c(1:No.Sets);
  } else {
      CalculatedSubsets = c(OnlySet);
  }
  for (set in CalculatedSubsets)
  {
    IndexRange = c(1:(nlevels(as.factor(OriginalColors))-1));
    if (UseAbs)
    {
        diss = 1-abs(cor(PCs[[set]]$data[, IndexRange], use = "p"));
    } else {
        diss = 1-cor(PCs[[set]]$data[, IndexRange], use = "p");
    }
    MEDiss[[set]] = list(Diss = diss);
  }
  
  if (is.null(OnlySet))
  {
    ConsDiss = (MEDiss[[1]]$Diss)
    if (No.Sets>1) for (set in 2:No.Sets)
      ConsDiss = pmax(ConsDiss, MEDiss[[set]]$Diss);
  } else {
    ConsDiss = MEDiss[[OnlySet]]$Diss;
  }
  
  METree = hclust(as.dist(ConsDiss), method = "average");
  METreeBranches = as.factor(ModuleNumber(HierTree = METree, CutHeight = CutHeight, MinSize = 1));

  # Analyze the branches: look for the ones that contain more than one original module
  
  MEUniqueBranches = levels(METreeBranches);
  MENo.Branches = nlevels(METreeBranches)
  MENumberOnBranch = rep(0, times = MENo.Branches);
  for (branch in 1:MENo.Branches)
  {
    MENumberOnBranch[branch] = sum(METreeBranches==MEUniqueBranches[branch]);
  }
  
  MergedColors = OriginalColors;
  
  # Merge modules on the same branch

  for (branch in 1:MENo.Branches) if (MENumberOnBranch[branch]>1)
  {
    ModulesOnThisBranch = names(METreeBranches)[METreeBranches==MEUniqueBranches[branch]];
    ColorsOnThisBranch = substring(ModulesOnThisBranch, 3);
    if (verbose>3) print.flush(paste(spaces, "  Merging original colors", paste(ColorsOnThisBranch, 
                                         collapse=", ")));
    for (color in 2:length(ColorsOnThisBranch))
      MergedColors[MergedColors==ColorsOnThisBranch[color]] = ColorsOnThisBranch[1];
  }
  
  No.Mods = nlevels(as.factor(MergedColors));
  RawModuleColors = levels(as.factor(MergedColors));

  # print(paste("No. of new modules: ", No.Mods));
  # print(paste("Merged module colors:"));
  # print(table(as.factor(MergedColors)));
 
  if (Relabel) 
  {
     # Relabel the merged colors to the usual order based on the number of genes in each module
 
     if (is.null(StandardColors))
     {
       StandardColors = c("turquoise","blue","brown","yellow","green","red","black","pink","magenta",
                   "purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan","grey60", 
                   "lightgreen", "lightyellow", "royalblue", "darkred", "darkgreen", "darkturquoise", 
                   "darkgrey", "orange", "darkorange", "white" );
   
     }
     
     No.GenesInModule = rep(0, No.Mods);
     for (mod in 1:No.Mods) No.GenesInModule[mod] = sum(MergedColors==RawModuleColors[mod]);
   
     # print(paste("No.GenesInModule: ", paste(No.GenesInModule, collapse=", ")));
     
     SortedRawModuleColors = RawModuleColors[order(-No.GenesInModule)]

     # print(paste("SortedRawModuleColors:", paste(SortedRawModuleColors, collapse=", ")));
     
     # Change the color names to the standard sequence, but leave grey grey (that's why rank in general does
     # not equal color)
 
     MergedNewColors = MergedColors;
 
     if (verbose>3) print(paste(spaces, "   Changing original colors:"));
     rank = 0;
     for (color in 1:length(SortedRawModuleColors)) if (SortedRawModuleColors[color]!="grey")
     {
       rank = rank + 1;
       if (verbose>3) print(paste(spaces, "      ", SortedRawModuleColors[color], 
                                  "to ", StandardColors[rank]));
       MergedNewColors[MergedColors==SortedRawModuleColors[color]] = StandardColors[rank];
     }
     # print("Table of new MergedColors:");
     # print(table(as.factor(MergedNewColors)));
  } else {
     MergedNewColors = MergedColors; 
  }

  if (CalculateNewPCs)
  {
     if (verbose>0) print.flush(paste(spaces, "  Calculating new PCs..."));
     NewPCs = ModulePCs(datE, Subsets, ModuleColors = MergedNewColors,
                        OnlySet = OnlySet, 
                        verbose = verbose-1, print.level = print.level+1);
  } else {
     NewPCs = NULL;
  }

  if (PCsInSingleFrame) 
  {
    NewPCs = NewPCs[[1]]$data;
    PCs = PCs[[1]]$data;
  }

  list(Colors = MergedNewColors, ClustTree = METree, CutHeight = CutHeight, OldPCs = PCs, NewPCs = NewPCs);
}





###############################################################################################################################
# I) GENERAL STATISTICAL FUNCTIONS

mean1=function(x) mean(x,na.rm=T)
var1=function(x) var(x,na.rm=T)



if(exists("scatterplot1") ) rm(scatterplot1);
scatterplot1=function(x,y, title1="",col1="black",xlab1=NA,ylab1=NA, cex1=1, cex.axis1=1.5,cex.lab1=1.5, cex.main1=1.5 ,ylim1=-1,correlationmethod="p" ){
if ( is.na(xlab1) ) xlab1= as.character(match.call(expand.dots = FALSE)$x)
if ( is.na(ylab1) ) ylab1= as.character(match.call(expand.dots = FALSE)$y)
x= as.numeric(as.character(x))
y= as.numeric(as.character(y))
cor1=signif(cor(x,y,use="p",method=correlationmethod),2)
corp=signif(cor.test(x,y,use="p",method=correlationmethod)$p.value,2)
if (corp<10^(-20) ) corp="<10^{-20}"
if ( length(ylim1)==2) {plot(x,y, main=paste(title1,"cor=", cor1,"p=",corp),col=col1,xlab=xlab1,ylab=ylab1, cex=cex1, 
cex.axis=cex.axis1,cex.lab=cex.lab1, cex.main=cex.main1,ylim=ylim1)} else {
plot(x,y, main=paste(title1,"cor=", cor1,"p=",corp),col=as.character(col1),xlab=xlab1,ylab=ylab1, cex=cex1, cex.axis=cex.axis1,cex.lab=cex.lab1, cex.main=cex.main1)
}
}



if  (exists("boxplot1") ) rm(boxplot1);
boxplot1=function(x,g,title1=" ",xlab1=NA,ylab1=NA, cex.axis1=1.5,cex.lab1=1.5, cex.main1=1.5,cex1=1.5,col1="white") {
if ( is.na(xlab1) ) xlab1= as.character(match.call(expand.dots = FALSE)$g)
print(xlab1)
if ( is.na(ylab1) ) ylab1= as.character( match.call(expand.dots = FALSE)$x)
print(ylab1)
p1=signif(kruskal.test(x, factor(g) )$p.value,2)
if (p1< 5.0*10^(-22) ) p1="<5.0x10^{-22}"
boxplot(x~factor(g), notch=T,varwidth=T, main=paste(title1,", p=",as.character(p1) ),col=col1,xlab=xlab1,ylab=ylab1, cex=cex1, cex.axis=cex.axis1,cex.lab=cex.lab1, cex.main=cex.main1)
}





# this function compute an asymptotic p-value for a given correlation (r) and sample size (n) 
if (exists("FisherTransformP") ) rm(FisherTransformP); 
FisherTransformP=function(r, n) {
#Z=sqrt(n-3) * 0.5*log((1+r)/(1-r)) 
#2*(1-pnorm(abs(Z) ))
# the following is implemented in the cor.test function
T=sqrt(n-2) * r/sqrt(1-r^2)
2*(1-pt(abs(T),n-2))
}


goodmargins1= c(0, 5, 2.3, 1)


# This function can be used to create an average linkage hierarchical clustering tree
# or the microarray samples. The rows of datExpr correspond to the samples and the columns to the genes
# You can optionally input a quantitative microarray sample trait.
if(exists("PlotClusterTreeSamples"))  rm(PlotClusterTreeSamples);
PlotClusterTreeSamples=function(datExpr, y=NULL) {
hierSamples=hclust( dist( datExpr  ), method="average" )
if (is.null(y) ) {plot(hierSamples, main="Cluster Tree of the Samples (average linkage, Euclidean dist.)", sub="")}
if (!is.null(y) ) {
if (!is.numeric(y) ) {warning("The microarray sample trait y is not numeric. Therefore, the function PlotClusterTreeSamples will transform it to a numeric variable");
y=as.numeric(as.character(y))
} # end of if (!is.numeric(y) )
if (  dim(as.matrix(datExpr))[[1]] != length(y) ) stop("Input Error: The number of microarray sample arrays does not match the length of the sample trait. Hint: consider transposing datExpr. Also make sure that y is a vector. dim(as.matrix(datExpr))[[1]] != length(y)")
par(mar=c(2,2,2,2), mfrow=c(2,1))
plot(hierSamples, main="Cluster Tree of the Samples (Euclidean dist.)", sub="")
if ( is.integer(y) ) {
if (min(y, na.rm=T)<0 ) y=y- min(y, na.rm=T)+1
hclustplotn(hierSamples, data.frame(y),main="Array Samples Colored by the Sample Trait")
} # end of if ( is.integer(y) )
if ( !is.integer(y) ) {
yBinary =as.integer(ifelse( y>=mean(y, na.rm=T), max(y,na.rm=T) , min(y,na.rm=T)))
if (min(yBinary, na.rm=T)<0 ) yBinary=yBinary- min(yBinary, na.rm=T)+1
hclustplotn(hierSamples, data.frame(yBinary),main="Samples Colored by the Dichotomized (Binary) Trait")
} # end of if ( !is.integer(y) )
} # end of if
}# end of function PlotClusterTreeSamples









