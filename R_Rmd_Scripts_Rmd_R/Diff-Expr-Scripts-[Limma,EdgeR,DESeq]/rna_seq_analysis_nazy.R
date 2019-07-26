source("http://bioconductor.org/biocLite.R")

rpackage.dir = "/project/umb_triley/Rpackages/"

dir = "/Users/celizabeth/Documents/UMB/Riley_Lab/HiSeq data/nazy/"
cts = read.table(paste(dir, sep = "", "counts_nazy.txt"))
# cts.s10 = read.table(paste(dir, sep = "", "counts_nazy_s10.txt"))
# cts = cbind(cts, cts.s10)

sel.rn = rowSums(cts) != 0 
sel.rn10 = rowSums(cts) > 10 

cts0 = cts[sel.rn,]
cts10 = cts[sel.rn10,]

# biocLite("edgeR", lib.loc=rpackage.dir)
library(edgeR, lib.loc = rpackage.dir)
# dgel = DGEList(counts=cts, group=factor(c(rep("HF_SAMP6", 22), rep(c("LF_SAMP6", "HF_AKR/J"), each = 24), rep("LF_AKR/J", times = 16))))
dgel10 = DGEList(counts=cts10, group=factor(c(rep("HF_SAMP6", 22), rep(c("LF_SAMP6", "HF_AKR/J"), each = 24), rep("LF_AKR/J", times = 16))))

# MA plots
hf_samp6_avg_count = apply(dgel$counts[,dgel$samples$group == "HF_SAMP6"], 1, function(x) mean(x))
lf_samp6_avg_count = apply(dgel$counts[,dgel$samples$group == "LF_SAMP6"], 1, function(x) mean(x))
hf_akrj_avg_count = apply(dgel$counts[,dgel$samples$group == "HF_AKR/J"], 1, function(x) mean(x))
lf_akrj_avg_count = apply(dgel$counts[,dgel$samples$group == "LF_AKR/J"], 1, function(x) mean(x))
 
# before normalization
maPlot(hf_samp6_avg_count,lf_samp6_avg_count,normalize=F,
       lowess=TRUE, ylim=c(-8,8),pch=19, cex=0.1)
abline(h=0, lty=2)
title("High-Fat SAMP6 vs. Low Fat SAMP6 before Normalization")

maPlot(hf_akrj_avg_count,lf_akrj_avg_count,normalize=F,
       lowess=TRUE, ylim=c(-8,8),pch=19, cex=0.1)
abline(h=0, lty=2)
title("High-Fat AKR/J vs. Low Fat AKR/J before Normalization")
maPlot(hf_samp6_avg_count,hf_akrj_avg_count,normalize=F,
       lowess=TRUE, ylim=c(-8,8),pch=19, cex=0.1)
abline(h=0, lty=2)
title("High-Fat SAMP6 vs. High-Fat AKR/J before Normalization")
maPlot(hf_samp6_avg_count,lf_akrj_avg_count,normalize=F,
       lowess=TRUE, ylim=c(-8,8),pch=19, cex=0.1)
abline(h=0, lty=2)
title("High-Fat SAMP6 vs. Low-Fat AKR/J before Normalization")
maPlot(lf_samp6_avg_count,lf_akrj_avg_count,normalize=F,
       lowess=TRUE, ylim=c(-8,8),pch=19, cex=0.1)
abline(h=0, lty=2)
title("Low-Fat SAMP6 vs. Low-Fat AKR/J before Normalization")
maPlot(lf_samp6_avg_count,hf_akrj_avg_count,normalize=F,
       lowess=TRUE, ylim=c(-8,8),pch=19, cex=0.1)
abline(h=0, lty=2)
title("Low-Fat SAMP6 vs. High-Fat AKR/J before Normalization")

# after normalization
maPlot(hf_samp6_avg_count,lf_samp6_avg_count,normalize=T,
       lowess=TRUE, ylim=c(-8,8),pch=19, cex=0.1)
abline(h=0, lty=2)
title("High-Fat SAMP6 vs. Low Fat SAMP6 after Normalization")

maPlot(hf_akrj_avg_count,lf_akrj_avg_count,normalize=T,
       lowess=TRUE, ylim=c(-8,8),pch=19, cex=0.1)
abline(h=0, lty=2)
title("High-Fat AKR/J vs. Low Fat AKR/J after Normalization")
maPlot(hf_samp6_avg_count,hf_akrj_avg_count,normalize=T,
       lowess=TRUE, ylim=c(-8,8),pch=19, cex=0.1)
abline(h=0, lty=2)
title("High-Fat SAMP6 vs. High-Fat AKR/J after Normalization")
maPlot(hf_samp6_avg_count,lf_akrj_avg_count,normalize=T,
       lowess=TRUE, ylim=c(-8,8),pch=19, cex=0.1)
abline(h=0, lty=2)
title("High-Fat SAMP6 vs. Low-Fat AKR/J after Normalization")
maPlot(lf_samp6_avg_count,lf_akrj_avg_count,normalize=T,
       lowess=TRUE, ylim=c(-8,8),pch=19, cex=0.1)
abline(h=0, lty=2)
title("Low-Fat SAMP6 vs. Low-Fat AKR/J after Normalization")
maPlot(lf_samp6_avg_count,hf_akrj_avg_count,normalize=T,
       lowess=TRUE, ylim=c(-8,8),pch=19, cex=0.1)
abline(h=0, lty=2)
title("Low-Fat SAMP6 vs. High-Fat AKR/J after Normalization")


# normalize
dgel10 = calcNormFactors(dgel10,method="TMM")

# estimate dispersion
dgel10 = estimateCommonDisp(dgel10)
dgel10 = estimateTagwiseDisp(dgel10)


# # Biological coefficient of variation vs. read counts
plotBCV(dgel10, cex=0.2, main = "BCV vs. Abundance")

# MDS plot
plotMDS(dgel)
title("Multidimenensional Scaling Plot after Normalization")


# Exact tests for differential gene expression.
# Groups are numbered as follows:
# 1 = HF_SAMP6
# 2 = LF_SAMP6
# 3 = HF_AKRJ
# 4 = LF_AKRJ
et1vs2.thresh10= exactTest(dgel10, pair = c("LF_SAMP6", "HF_SAMP6"))
et1vs3.thresh10= exactTest(dgel10, pair = c("HF_SAMP6", "HF_AKR/J"))
et2vs3.thresh10= exactTest(dgel10, pair = c("LF_SAMP6", "HF_AKR/J"))
et3vs4.thresh10= exactTest(dgel10, pair = c("LF_AKR/J", "HF_AKR/J"))
et1vs4.thresh10= exactTest(dgel10, pair = c("LF_AKR/J","HF_SAMP6"))
et2vs4.thresh10= exactTest(dgel10, pair = c("LF_SAMP6", "LF_AKR/J"))


# Make a table of each set of differential expression data
table12 = data.frame(topTags(et1vs2, n = nrow(et1vs2$table), adjust.method = "BY"))
table12 = cbind(rownames(table12), table12)
rownames(table12) = NULL
names(table12) = c("ID", "logFC", "logCPM", "PValue", "FDR")

table34 = data.frame(topTags(et3vs4.thresh10, n = nrow(et3vs4.thresh10$table), adjust.method = "BY"))
table34 = cbind(rownames(table34), table34)
rownames(table34) = NULL
names(table34) = c("ID", "logFC", "logCPM", "PValue", "FDR")

table34 = data.frame(topTags(et3vs4, n = nrow(et3vs4$table), adjust.method = "BY"))
table34 = cbind(rownames(table34), table34)
rownames(table34) = NULL
names(table34) = c("ID", "logFC", "logCPM", "PValue", "FDR")

table13= data.frame(topTags(et1vs3.thresh10, n = nrow(et1vs3.thresh10$table), adjust.method = "BY"))
table13 = cbind(rownames(table13), table13)
rownames(table13) = NULL
names(table13) = c("ID", "logFC", "logCPM", "PValue", "FDR")

table23 = data.frame(topTags(et2vs3.thresh10, n = nrow(et2vs3.thresh10$table), adjust.method = "BY"))
table23 = cbind(rownames(table23), table23)
rownames(table23) = NULL
names(table23) = c("ID", "logFC", "logCPM", "PValue", "FDR")

table14 = data.frame(topTags(et1vs4.thresh10, n = nrow(et1vs4.thresh10$table), adjust.method = "BY"))
table14 = cbind(rownames(table14), table14)
rownames(table14) = NULL
names(table14) = c("ID", "logFC", "logCPM", "PValue", "FDR")

table24 = data.frame(topTags(et2vs4.thresh10, n = nrow(et2vs4.thresh10$table), adjust.method = "BY"))
table24 = cbind(rownames(table24), table24)
rownames(table24) = NULL
names(table24) = c("ID", "logFC", "logCPM", "PValue", "FDR")


# plot smears
de.1vs2.001 = decideTestsDGE(et1vs2.thresh10, adjust.method="BY", p.value=0.001)
summary(de.1vs2.001)
plotSmear(et1vs2.thresh10, de.tags=rownames(et1vs2.thresh10$table)[as.logical(de.1vs2.001)])
# # identify(cpm(et1vs2.thresh20$table$logFC, log=T), et1vs2.thresh20$table$logFC,rownames(et1vs2.thresh20$table))
title("High-Fat SAMP6 vs Low-Fat SAMP6, p < 0.001")

de.3vs4.001 = decideTestsDGE(et3vs4.thresh10, adjust.method="BY", p.value=0.001)
plotSmear(et3vs4.thresh10, de.tags=rownames(et3vs4.thresh10$table)[as.logical(de.3vs4.001)])
title("High-Fat AKR/J vs Low-Fat AKR/J, p < 0.001")

de.1vs3.001 = decideTestsDGE(et1vs3.thresh10, adjust.method="BY", p.value=0.001)
plotSmear(et1vs3.thresh10, de.tags=rownames(et1vs3.thresh10$table)[as.logical(de.1vs3.001)])
title("High-Fat SAMP6 vs High-Fat AKR/J, p < 0.001")

de.2vs3.001 = decideTestsDGE(et2vs3.thresh10, adjust.method="BY", p.value=0.001)
plotSmear(et2vs3.thresh10, de.tags=rownames(et2vs3.thresh10$table)[as.logical(de.2vs3.001)])
title("Low-Fat SAMP6 vs High-Fat AKR/J, p < 0.001")

de.1vs4.001 = decideTestsDGE(et1vs4.thresh10, adjust.method="BY", p.value=0.001)
plotSmear(et1vs4.thresh10, de.tags=rownames(et1vs4.thresh10$table)[as.logical(de.1vs4.001)])
title("High-Fat SAMP6 vs. Low-Fat AKR/J, p < 0.001")

de.2vs4.001 = decideTestsDGE(et2vs4.thresh10, adjust.method="BY", p.value=0.001)
plotSmear(et2vs4.thresh10, de.tags=rownames(et2vs4.thresh10$table)[as.logical(de.1vs4.001)])
title("Low-Fat SAMP6 vs. Low-Fat AKR/J, p < 0.001")


# Add Fold Change column (exponentiated logFC) to expression tables
table12$FC = 2 ^ table12$logFC
table13$FC = 2 ^ table13$logFC
table23$FC = 2 ^ table23$logFC
table34$FC = 2 ^ table34$logFC
table14$FC = 2 ^ table14$logFC
table24$FC = 2 ^ table24$logFC


# create smaller tables of genes with pvalue < 0.001
table12.001 = table12[table12$PValue < .001,]
table13.001 = table13[table13$PValue < .001,]
table23.001 = table23[table23$PValue < .001,]
table34.001 = table34[table34$PValue < .001,]
table14.001 = table14[table14$PValue < .001,]
table24.001 = table24[table24$PValue < .001,]



biocLite("org.Mm.eg.db")
library("org.Mm.eg.db")
ls("package:org.Mm.eg.db")

getGeneSymbolAcc = function (table) 
{
  table$Symbol = NA
  table$Accession = NA
  for (i in 1:length(table$ID))
  {
    a = as.character(table$ID[i])
    b = tryCatch(get(a, org.Mm.egSYMBOL), error = function(e) e)
    if(inherits(b, "error")) next
    c = tryCatch(get(a, org.Mm.egREFSEQ), error = function(e) e)
    if(inherits(c, "error")) next
    table$Symbol[i] = b
    table$Accession[i] = c
  }
  return(table)
}


biocLite("GO.db")
library("GO.db")
ls("package:GO.db")

getGOTerms = function (table) {
  go = apply(table, 1, function (x) tryCatch(get(x, org.Mm.egGO), error = function(e) e))
  names(go) = table$ID
  
  table$GOTerms = NA
  for (i in 1:nrow(table))
  {
    eid = as.character(table$ID[i])
    print(paste("next eid: ", eid))
    if (!is.na(go[eid])) {
      for (j in 1:length(go[eid][[1]]))
      {
        term = get(go[eid][[1]][[j]]$GOID, GOTERM)
        print(Term(term))  
        if (is.na(table$GOTerms[i]))
          table$GOTerms[i] = paste0(Term(term), ", ")
        else
          table$GOTerms[i] = paste0(table$GOTerms[i], Term(term), sep =", ")
      }
    }
  }
  return (table)
}


table12.001 = getGeneSymbolAcc(table12.001)
table13.001 = getGeneSymbolAcc(table13.001)
table23.001 = getGeneSymbolAcc(table23.001)
table34.001 = getGeneSymbolAcc(table34.001)
table14.001 = getGeneSymbolAcc(table14.001)
table24.001 = getGeneSymbolAcc(table24.001)

table12.001 = getGOTerms(table12.001)
table13.001 = getGOTerms(table13.001)
table23.001 = getGOTerms(table23.001)
table34.001 = getGOTerms(table34.001)
table14.001 = getGOTerms(table14.001)
table24.001 = getGOTerms(table24.001)

# Excel file of results
library(xlsx)
wb = createWorkbook()
sheet1 = createSheet(wb, sheetName = "HF_SAMP6 vs. LF_SAMP6")
addDataFrame(table12, sheet1)
sheet2 = createSheet(wb, sheetName = "HF_SAMP6 vs. HF_AKR/J")
addDataFrame(table13, sheet2)
sheet3 = createSheet(wb, sheetName = "LF_SAMP6 vs. HF_AKR/J")
addDataFrame(table23, sheet3)
sheet4 = createSheet(wb, sheetName = "LF_AKR/J vs. HF_AKR/J")
addDataFrame(table34, sheet4)
sheet5 = createSheet(wb, sheetName = "HF_SAMP6 vs. LF_AKR/J")
addDataFrame(table14, sheet5)
sheet6 = createSheet(wb, sheetName = "LF_SAMP6 vs. LF_AKR/J")
addDataFrame(table24, sheet4)
saveWorkbook(wb, paste(dir, sep = "", "thresh10/edger_table_.001_thresh10.xlsx"))


wb = createWorkbook()
sheet1 = createSheet(wb, sheetName = "HF_SAMP6 vs. LF_SAMP6")
addDataFrame(table12, sheet1)
sheet2 = createSheet(wb, sheetName = "HF_SAMP6 vs. HF_AKR/J")
addDataFrame(table13, sheet2)
sheet3 = createSheet(wb, sheetName = "LF_SAMP6 vs. HF_AKR/J")
addDataFrame(table23, sheet3)
sheet4 = createSheet(wb, sheetName = "LF_AKR/J vs. HF_AKR/J")
addDataFrame(table34, sheet4)
sheet5 = createSheet(wb, sheetName = "HF_SAMP6 vs. LF_AKR/J")
addDataFrame(table14, sheet5)
sheet6 = createSheet(wb, sheetName = "LF_SAMP6 vs. LF_AKR/J")
addDataFrame(table24, sheet4)
saveWorkbook(wb, paste(dir, sep = "", "thresh10/edger_table_.0001_thresh10.xlsx"))

# txt table of results
# write.table(table12.001, paste(dir, sep = "", "thresh10/edger_HF_SAMP6_vs_LF_SAMP6_.001_thresh10.txt"))
# write.table(table13.001, paste(dir, sep = "", "thresh10/edger_HF_SAMP6_vs_HF_AKRJ_.001_thresh10.txt"))
# write.table(table23.001, paste(dir, sep = "", "thresh10/edger_LF_SAMP6_vs_HF_AKRJ_.001_thresh10.txt"))
# write.table(table34.001, paste(dir, sep = "", "thresh10/edger_HF_AKRJ_vs_LF_AKRJ_.001_thresh10.txt"))
# write.table(table14.001, paste(dir, sep = "", "thresh10/edger_HF_SAMP6_vs_LF_AKRJ_.001_thresh10.txt"))
# write.table(table24.001, paste(dir, sep = "", "thresh10/edger_LF_SAMP6_vs_LF_AKRJ_.001_thresh10.txt"))
# 
# write.table(table12.001[table12.001$PValue < .0001,], paste(dir, sep = "", "thresh10/edger_HF_SAMP6_vs_LF_SAMP6_.0001_thresh10.txt"))
# write.table(table13.001[table13.001$PValue < .0001,], paste(dir, sep = "", "thresh10/edger_HF_SAMP6_vs_HF_AKRJ_.0001_thresh10.txt"))
# write.table(table23.001[table23.001$PValue < .0001,], paste(dir, sep = "", "thresh10/edger_LF_SAMP6_vs_HF_AKRJ_.0001_thresh10.txt"))
# write.table(table34.001[table34.001$PValue < .0001,], paste(dir, sep = "", "thresh10/edger_HF_AKRJ_vs_LF_AKRJ_.0001_thresh10.txt"))
# write.table(table14.001[table14.001$PValue < .0001,], paste(dir, sep = "", "thresh10/edger_HF_SAMP6_vs_LF_AKRJ_.0001_thresh10.txt"))
# write.table(table24.001[table24.001$PValue < .0001,], paste(dir, sep = "", "thresh10/edger_LF_SAMP6_vs_LF_AKRJ_.0001_thresh10.txt"))
