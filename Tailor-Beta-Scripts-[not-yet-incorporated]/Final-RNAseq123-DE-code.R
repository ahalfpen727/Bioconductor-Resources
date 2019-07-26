

library(RNAseq123)
library(goseq)
groups
DGEobj=DGEList(G.rep.cnt.matrix,lib.size=colSums(G.rep.cnt.matrix),group=groups)
DGEobj
design <- model.matrix(~0+groups, data=DGEobj$samples)
colnames(design) <- levels(DGEobj$samples$group)
design
norm.factor <- calcNormFactors(DGEobj)
norm.factor
y <- estimateDisp(norm.factor,design)
y
est.common.disp <- estimateCommonDisp(y)
est.common.disp
est.tag.disp <- estimateTagwiseDisp(est.common.disp)
est.tag.disp

# perform quasi-likelihood F-tests:
fit <- glmQLFit(est.tag.disp,design)
fit
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)

################################
DGEobj=DGEList(g.rep.fpkm.ma,lib.size=colSums(g.rep.fpkm.ma),group=groups)
DGEobj
y <- estimateDisp(DGEobj, design, robust=TRUE)
y$common.dispersion
plotBCV(y)
# fit genewise glm's
fit <- glmFit(y, design)
# likelihood ratio tests for ctrl vs luts differences
lrt <- glmLRT(fit)
topTags(lrt)
o <- order(lrt$table$PValue)
# counts per million
cpm(y)[o[1:10],]
summary(decideTests(lrt))

plotMD(lrt)
abline(h=c(-1, 1), col="blue")
go <- goana(lrt)
topGO(go, ont="BP", sort="Up", n=30, truncate=30)

################################
DGEobj=DGEList(G.rep.cpm.ma,lib.size=colSums(G.rep.cpm.ma),group=groups)
DGEobj
y <- estimateDisp(DGEobj, design, robust=TRUE)
y$common.dispersion
plotBCV(y)
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)
# perform likelihood ratio tests:
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)

FDR <- p.adjust(qlf$table$PValue, method="BH")
sum(FDR < 0.05)

fit <- glmQLFit(y, design, robust=TRUE)
qlf <- glmQLFTest(fit)
topTags(qlf)
FDR <- p.adjust(qlf$table$PValue, method="BH")
sum(FDR < 0.05)
top <- rownames(topTags(qlf))
cpm(y)[top,]
summary(decideTests(qlf))

plotMD(qlf)
abline(h=c(-1,1), col="blue")
# Fold-change threshold for significance
fit <- glmQLFit(y, design)
tr <- glmTreat(fit, coef=2, lfc=1)
topTags(tr)

# Gene Ontology
fit <- glmQLFit(y, design)
fit
qlf <- glmQLFTest(fit, coef=2)
qlf
go <- goana(qlf, species="Hs")
go
topGO(go, sort="up")
keg <- kegga(qlf, species="Hs")
topKEGG(keg, sort="up")

# Classic Approach
DGEobj=DGEList(G.rep.cnt.matrix,group=grp)
DGEobj
y <- estimateDisp(DGEobj,design)

levels(y$samples$group)
et<-exactTest(y, pair=c(under,over))
topTags(et)

design <- model.matrix(~group, data=y$samples)
design

filter <- apply(g.cnt.matrix,1,function(x) mean(x)>5)
table(filter)

design <- model.matrix(~0 + groups, data=dgeObj$samples)
colnames(design) <- levels(conditions)
design


# edger Filter
keep <- rowSums(G.rep.cnt.matrix)>=10
nkeep <- sum(keep)
counts2 <- G.rep.cnt.matrix[keep,]
str(counts2)
# edgeR glm
#i <- 9

exactTst <- exactTest(est.tag.disp)
exactTst$table
o <- order(exactTst$table$PValue)
o.sig <- subset(exactTst$table, (PValue < 0.05))

z <- estimateGLMTrendedDisp(G.rep.cnt.matrix,design)
z

fite <- glmFit(G.rep.cnt.matrix,design,dispersion=z)
lrt <- glmLRT(fite)
o <- order(lrt$table$PValue)
o.sig <- subset(lrt$table, (PValue < 0.01))

sum(p.adjust(lrt$table$PValue,method="BH")<0.01)
topTags(exactTst)
# voom - ranked by lods
i <- 1
y <- voom(counts2,design,plot=FALSE)
fit <- lmFit(y,design)
fit <- eBayes(fit)
o <- order(fit$lods[,2], decreasing=TRUE)
sum(p.adjust(fit$p.value[,2],method="BH")<0.1)


# limma trend - ranked by lods
y <- cpm(counts2,log=TRUE,prior.count=1)
fit <- lmFit(y,design,weights=NULL)
fit <- eBayes(fit,trend=TRUE)
o <- order(fit$lods[,2], decreasing=TRUE)

# limma notrend - ranked by lods
fit <- eBayes(fit,trend=FALSE)
o <- order(fit$lods[,2], decreasing=TRUE)

# t-test
t.ord <- fit$coef[,2]/fit$stdev.unscaled[,2]/fit$sigma
p.ord <- pt(abs(t.ord),df=4,lower.tail=FALSE)*2
fdr.ord <- p.adjust(p.ord,method="BH")
o <- order(p.ord)
sum(fdr.ord<0.1)


# Note: using the devel versions of both packages!
library(DESeq) # version 1.9.11
library(edgeR) # version 2.99.8
library(VennDiagram)

# Read in data ------------------------------------------------------------

## Use pasilla data
datafile = system.file( "extdata/pasilla_gene_counts.tsv", package="pasilla" )
datafile

## Read in the data making the row names the first column
counttable <- read.table(datafile, header=T, row.names=1)
head(counttable)

## Make metadata data.frame
meta <- data.frame(
   row.names=colnames(counttable),
   condition=c("untreated", "untreated", "untreated", "untreated", "treated", "treated", "treated"),
   libType=c("single", "single", "paired", "paired", "single", "paired", "paired"))
meta$condition <- relevel(meta$condition, ref="untreated")
meta

## Independent filtering?
# keep_cpm <- rowSums(cpm(counttable)>2) >=2
# keep_quantile <- rowSums(counttable)>quantile(rowSums(counttable), probs=.5)
# addmargins(table(keep_cpm, keep_quantile))
# counttable <- counttable[keep_cpm, ]

# DESeq -------------------------------------------------------------------

## Make a new countDataSet
d <- newCountDataSet(counttable, meta)

## Estimate library size and dispersion
d <- estimateSizeFactors(d)
d <- estimateDispersions(d)
plotDispEsts(d, main="DESeq: Per-gene dispersion estimates")

## Principal components biplot on variance stabilized data, color-coded by condition-librarytype
print(plotPCA(varianceStabilizingTransformation(d), intgroup=c("condition", "libType")))

## Fit full and reduced models, get p-values
dfit1 <- fitNbinomGLMs(d, count~libType+condition)
dfit0 <- fitNbinomGLMs(d, count~libType)
dpval <- nbinomGLMTest(dfit1, dfit0)
dpadj <- p.adjust(dpval, method="BH")

## Make results table with pvalues and adjusted p-values
dtable <- transform(dfit1, pval=dpval, padj=dpadj)
dtable <- dtable[order(dtable$padj), ]
head(dtable)

# edgeR -------------------------------------------------------------------

## Make design matrix
condition <- relevel(factor(meta$condition), ref="untreated")
libType <- factor(meta$libType)
edesign <- model.matrix(~libType+condition)

## Make new DGEList, normalize by library size, and estimate dispersion allowing possible trend with average count size
e <- DGEList(counts=counttable)
e <- calcNormFactors(e)
e <- estimateGLMCommonDisp(e, edesign)
e <- estimateGLMTrendedDisp(e, edesign)
e <- estimateGLMTagwiseDisp(e, edesign)

## MDS Plot
plotMDS(e, main="edgeR MDS Plot")

## Biological coefficient of variation plot
plotBCV(e, cex=0.4, main="edgeR: Biological coefficient of variation (BCV) vs abundance")

## Fit the model, testing the coefficient for the treated vs untreated comparison
efit <- glmFit(e, edesign)
efit <- glmLRT(efit, coef="conditiontreated")

## Make a table of results
etable <- topTags(efit, n=nrow(e))$table
etable <- etable[order(etable$FDR), ]
head(etable)

## ~MA Plot
with(etable, plot(logCPM, logFC, pch=20, main="edgeR: Fold change vs abundance"))
with(subset(etable, FDR<0.05), points(logCPM, logFC, pch=20, col="red"))
abline(h=c(-1,1), col="blue")

# Comparison --------------------------------------------------------------

head(etable)
head(dtable)

addmargins(table(sig.edgeR=etable$FDR<0.05, sig.DESeq=dtable$padj<0.05))

merged <- merge(etable, dtable, by='row.names')
with(                     merged, plot(logFC, conditiontreated, xlab="logFC edgeR", ylab="logFC DESeq", pch=20, col="black", main="Fold change for DESeq vs edgeR"))
with(subset(merged, FDR<0.05),  points(logFC, conditiontreated, xlab="logFC edgeR", ylab="logFC DESeq", pch=20, col="red"))
with(subset(merged, padj<0.05), points(logFC, conditiontreated, xlab="logFC edgeR", ylab="logFC DESeq", pch=20, col="green"))
legend("topleft", xjust=1, yjust=1, legend=c("FDR<0.05 edgeR only", "FDR<0.05 DESeq & edgeR", "FDR>0.05"), pch=20, col=c("red", "green", "black"), bty="n")

