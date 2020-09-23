## ----globalOptions, include=FALSE----------------------------------------
# comment / uncomment to show / hide R code
knitr::opts_chunk$set(echo = TRUE)

## ----loadLib, results='hide', message=FALSE------------------------------
## Do Not Run
# how to install required packages?
# install.packages(c("RColorBrewer","mixOmics"))
# source("http://bioconductor.org/biocLite.R")
# biocLite("edgeR")
## End DNR
library(edgeR)
library(limma)
library(RColorBrewer)
library(mixOmics)
library(VennDiagram)
library(HTSFilter)

## ----importData----------------------------------------------------------
rawCountTable <- read.table("countData.txt", header=TRUE, sep="\t", row.names=1)
sampleInfo <- read.table("design.csv", header=TRUE, sep=",", row.names=1)

## ----checkData-----------------------------------------------------------
head(rawCountTable)
nrow(rawCountTable)

## ----checkDesign---------------------------------------------------------
head(sampleInfo)
sampleInfo$genotype <- as.factor(sampleInfo$genotype)

## ----DGEListCreation-----------------------------------------------------
dgeFull <- DGEList(rawCountTable, remove.zeros = TRUE)
dgeFull

## ----DGEListSamples------------------------------------------------------
dgeFull$samples$condition <- relevel(sampleInfo$condition, ref = "WT")
dgeFull$samples$genotype <- sampleInfo$genotype
dgeFull

## ----pseudoCounts--------------------------------------------------------
pseudoCounts <- log2(dgeFull$counts + 1)
head(pseudoCounts)

## ----histoPseudoCounts---------------------------------------------------
hist(pseudoCounts[ ,"Cond.WT.Rep.1"], main = "", xlab = "counts")

## ----boxplotPseudoCounts-------------------------------------------------
par(mar = c(8,4,1,2))
boxplot(pseudoCounts, col = "gray", las = 3, cex.names = 1)

## ----maPlotPseudoCounts, fig.width=10------------------------------------
limma::plotMA(pseudoCounts[ ,1:2], xlab = "M", ylab = "A", main = "")
abline(h = 0, col = "red")

## ----MDSPseudoCounts-----------------------------------------------------
colConditions <- brewer.pal(3, "Set2")
colConditions <- colConditions[match(sampleInfo$condition,
                                     levels(sampleInfo$condition))]
pchGenotypes <- c(8, 15, 16)[match(sampleInfo$genotype,
                                   levels(sampleInfo$genotype))]
plotMDS(pseudoCounts, pch = pchGenotypes, col = colConditions)
legend("topright", lwd = 2, col = brewer.pal(3, "Set2")[1:2], 
       legend = levels(sampleInfo$condition))
legend("bottomright", pch = c(8, 15, 16), 
       legend = levels(sampleInfo$genotype))

## ----cimPseudoCounts, fig.width=10, fig.height=10------------------------
sampleDists <- as.matrix(dist(t(pseudoCounts)))
sampleDists
cimColor <- colorRampPalette(rev(brewer.pal(9, "Reds")))(16)
cim(sampleDists, color = cimColor, symkey = FALSE, row.cex = 0.7, col.cex = 0.7)

## ----estimateNormFactors-------------------------------------------------
dgeFull <- calcNormFactors(dgeFull, method="TMM")
dgeFull

## ----countsAfterNorm-----------------------------------------------------
head(dgeFull$counts)

## ----sampleInfoAfterNorm-------------------------------------------------
dgeFull$samples

## ----estimateNormCounts, fig.width=10, fig.height=10---------------------
normCounts <- cpm(dgeFull)
pseudoNormCounts <- cpm(dgeFull, log = TRUE, prior.count = 1)
par(mar = c(8,4,1,2))
boxplot(pseudoNormCounts, col = "gray", las = 3, cex.names = 1)

## ----normMDS-------------------------------------------------------------
plotMDS(pseudoNormCounts, pch = pchGenotypes, col = colConditions)
legend("topright", lwd = 2, col = brewer.pal(3, "Set2")[1:2], 
       legend = levels(sampleInfo$condition))
legend("bottomright", pch = c(8, 15, 16), 
       legend = levels(sampleInfo$genotype))

## ----DGEListCreationG----------------------------------------------------
dgeFull.group <- DGEList(rawCountTable, remove.zeros = TRUE, 
                         group = dgeFull$samples$condition)
dgeFull.group$samples$norm.factors <- dgeFull$samples$norm.factors
dgeFull.group

## ----estimateDispersion, cache=TRUE--------------------------------------
dgeFull.group <- estimateCommonDisp(dgeFull.group)
dgeFull.group <- estimateTagwiseDisp(dgeFull.group)
dgeFull.group

## ----plotBCV-------------------------------------------------------------
plotBCV(dgeFull.group)

## ----fisherExact---------------------------------------------------------
dgeExactTest <- exactTest(dgeFull.group)
dgeExactTest

## ----topTag--------------------------------------------------------------
resExactTest <- topTags(dgeExactTest, n = nrow(dgeExactTest$table))
head(resExactTest$table)

## ----histPValExact, fig.width=10, fig.height=6---------------------------
par(mfrow = c(1,2))
hist(resExactTest$table$PValue, xlab = "p-value", main = "raw p-values")
hist(resExactTest$table$FDR, xlab = "p-value", main = "adjusted p-values")

## ----DEGExact------------------------------------------------------------
selectedET <- resExactTest$table$FDR < 0.05 & abs(resExactTest$table$logFC) > 1
selectedET <- resExactTest$table[selectedET, ]
nrow(selectedET)
head(selectedET)

## ----UDregExact----------------------------------------------------------
selectedET$updown <- factor(ifelse(selectedET$logFC > 0, "up", "down"))
head(selectedET)

## ----exportExact---------------------------------------------------------
write.table(selectedET, file = "tomatoDEG.csv", sep = ",")

## ----designMatrix--------------------------------------------------------
design.matrix <- model.matrix(~ dgeFull$samples$condition + 
                                dgeFull$samples$genotype)
design.matrix

## ----estimateDispersionGLM, cache=TRUE-----------------------------------
dgeFull <- estimateDisp(dgeFull, design.matrix)
dgeFull

## ----plotBCVGLM----------------------------------------------------------
plotBCV(dgeFull)

## ----GLMfit--------------------------------------------------------------
fit <- glmFit(dgeFull, design.matrix)
fit

## ----LRTtest-------------------------------------------------------------
dgeLRTtest <- glmLRT(fit, coef = 2)
dgeLRTtest

## ----LRTtest2------------------------------------------------------------
contrasts <- rep(0, ncol(design.matrix))
contrasts[3] <- 1
contrasts[4] <- -1
dgeLRTtest2 <- glmLRT(fit, contrast = contrasts)
dgeLRTtest2

## ----topTagsGLM----------------------------------------------------------
resLRT <- topTags(dgeLRTtest, n = nrow(dgeFull$counts))
head(resLRT$table)

## ----extractDEGGLM-------------------------------------------------------
selectedLRT <- resLRT$table$FDR < 0.05 & abs(resLRT$table$logFC) > 1
selectedLRT <- resLRT$table[selectedLRT, ]
nrow(selectedLRT)
head(selectedLRT)

## ----VennDiagram---------------------------------------------------------
vd <- venn.diagram(x = list("Exact test" = rownames(selectedET),
                            "GLM" = rownames(selectedLRT)),
                   fill = brewer.pal(3, "Set2")[1:2], filename = NULL)
grid.draw(vd)

## ----filter, cache=TRUE--------------------------------------------------
dgeFilt <- HTSFilter(dgeFull)$filteredData
dgeFilt

## ----GLMfilt-------------------------------------------------------------
fit <- glmFit(dgeFilt, design.matrix)
dgeLRTfilt <- glmLRT(fit, coef = 2)
resLRTfilt <- topTags(dgeLRTfilt, n = nrow(dgeFilt$counts))
selectedFilt <- resLRTfilt$table$FDR < 0.05 & abs(resLRTfilt$table$logFC) > 1
selectedFilt <- resLRTfilt$table[selectedFilt, ]
nrow(selectedFilt)
head(selectedFilt)

## ----VennDiagramFilt-----------------------------------------------------
vd <- venn.diagram(x = list("No filtering" = rownames(selectedLRT),
                            "Filtering" = rownames(selectedFilt)),
                   fill = brewer.pal(3, "Set2")[1:2], filename = NULL)
grid.draw(vd)

## ----MADEG---------------------------------------------------------------
dgeFilt$samples$group <- dgeFilt$samples$condition
plotSmear(dgeFilt, de.tags = rownames(selectedFilt))

## ----volcanoPlot---------------------------------------------------------
volcanoData <- cbind(resLRTfilt$table$logFC, -log10(resLRTfilt$table$FDR))
colnames(volcanoData) <- c("logFC", "negLogPval")
DEGs <- resLRTfilt$table$FDR < 0.05 & abs(resLRTfilt$table$logFC) > 1
point.col <- ifelse(DEGs, "red", "black")
plot(volcanoData, pch = 16, col = point.col, cex = 0.5)

## ----selectDEHeatmap, fig.width=10, fig.height=10------------------------
selY <- cpm(dgeFilt, log = TRUE, prior.count = 1)
selY <- selY[match(rownames(selectedFilt), rownames(dgeFilt$counts)), ]
finalHM <- cim(t(selY), color = cimColor, symkey = FALSE, row.cex = 0.7,
               col.cex = 0.7)

## ----plotDendoGenes, fig.width=10, fig.height=10-------------------------
plot(finalHM$ddc, leaflab="none")
abline(h=10, lwd=2, col="pink")

## ----cutDendoGenes-------------------------------------------------------
geneClust <- cutree(as.hclust(finalHM$ddc), h=10)
head(geneClust)

## ----nbClustGene---------------------------------------------------------
length(unique(geneClust))

## ----geneClust1----------------------------------------------------------
names(which(geneClust == 1))

## ----sessionInformation, echo=TRUE---------------------------------------
sessionInfo()

## ----deleteVD, echo=FALSE------------------------------------------------
system("rm VennDiagram*.log")

