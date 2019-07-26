# perform GO analysis, write HTML output to file for each set of expression data and GO category
dir = "/Users/celizabeth/Documents/UMB/Riley_Lab/HiSeq data/nazy/thresh10/"

# read in tables
library(xlsx)

# increase jvm memory, otherwise it crashes half the time
# options(java.parameters = "-Xmx3072m")

table12 = read.xlsx2(paste(dir, "edger_table_.0001_thresh10.xlsx", sep = ""),
                     sheetIndex = 1, colIndex = 2:10)
table13 = read.xlsx2(paste(dir, "edger_table_.0001_thresh10.xlsx", sep = ""),
                     sheetIndex = 2, colIndex = 2:10)
table23 = read.xlsx2(paste(dir, "edger_table_.0001_thresh10.xlsx", sep = ""),
                     sheetIndex = 3, colIndex = 2:10)
table34 = read.xlsx2(paste(dir, "edger_table_.0001_thresh10.xlsx", sep = ""),
                     sheetIndex = 4, colIndex = 2:10)
table14 = read.xlsx2(paste(dir, "edger_table_.0001_thresh10.xlsx", sep = ""),
                     sheetIndex = 5, colIndex = 2:10)
table24 = read.xlsx2(paste(dir, "edger_table_.0001_thresh10.xlsx", sep = ""),
                     sheetIndex = 6, colIndex = 2:10)

# biocLite("org.Mm.eg.db")
library("org.Mm.eg.db")

# biocLite("GOstats")
library(GOstats)

entrez_object = org.Mm.egGO
# Get the entrez gene identifiers that are mapped to a GO ID
mapped_genes = mappedkeys(entrez_object)

# Cellular component
params12.cc = new('GOHyperGParams',
                   geneIds=as.character(table12$ID),
                   universeGeneIds=mapped_genes,
                   ontology=c('CC'),
                   pvalueCutoff=0.001,
                   conditional=F,
                   testDirection='over',
                   annotation="org.Mm.eg.db"
)
hgOver12.cc = hyperGTest(params12.cc)
result12.cc = summary(hgOver12.cc)
result12.cc$Adj.Pvalue = p.adjust(result12.cc$Pvalue, method = 'BH')
result12.cc = cbind(result12.cc[1:2], result12.cc["Adj.Pvalue"], result12.cc[3:7])
head(result12.cc)

params13.cc = new('GOHyperGParams',
                   geneIds=as.character(table13$ID),
                   universeGeneIds=mapped_genes,
                   ontology=c('CC'),
                   pvalueCutoff=0.001,
                   conditional=F,
                   testDirection='over',
                   annotation="org.Mm.eg.db"
)
hgOver13.cc = hyperGTest(params13.cc)
result13.cc = summary(hgOver13.cc)
result13.cc$Adj.Pvalue = p.adjust(result13.cc$Pvalue, method = 'BH')
result13.cc = cbind(result13.cc[1:2], result13.cc["Adj.Pvalue"], result13.cc[3:7])
head(result13.cc)

params23.cc = new('GOHyperGParams',
                   geneIds=as.character(table23$ID),
                   universeGeneIds=mapped_genes,
                   ontology=c('CC'),
                   pvalueCutoff=0.001,
                   conditional=F,
                   testDirection='over',
                   annotation="org.Mm.eg.db"
)
hgOver23.cc = hyperGTest(params23.cc)
result23.cc = summary(hgOver23.cc)
result23.cc$Adj.Pvalue = p.adjust(result23.cc$Pvalue, method = 'BH')
result23.cc = cbind(result23.cc[1:2], result23.cc["Adj.Pvalue"], result23.cc[3:7])
head(result23.cc)

params34.cc = new('GOHyperGParams',
                   geneIds=as.character(table34$ID),
                   universeGeneIds=mapped_genes,
                   ontology=c('CC'),
                   pvalueCutoff=0.001,
                   conditional=F,
                   testDirection='over',
                   annotation="org.Mm.eg.db"
)
hgOver34.cc = hyperGTest(params34.cc)
result34.cc = summary(hgOver34.cc)
result34.cc$Adj.Pvalue = p.adjust(result34.cc$Pvalue, method = 'BH')
result34.cc = cbind(result34.cc[1:2], result34.cc["Adj.Pvalue"], result34.cc[3:7])
head(result34.cc)

params14.cc = new('GOHyperGParams',
                   geneIds=as.character(table14$ID),
                   universeGeneIds=mapped_genes,
                   ontology=c('CC'),
                   pvalueCutoff=0.001,
                   conditional=F,
                   testDirection='over',
                   annotation="org.Mm.eg.db"
)
hgOver14.cc = hyperGTest(params14.cc)
result14.cc = summary(hgOver14.cc)
result14.cc$Adj.Pvalue = p.adjust(result14.cc$Pvalue, method = 'BH')
result14.cc = cbind(result14.cc[1:2], result14.cc["Adj.Pvalue"], result14.cc[3:7])
head(result14.cc)

params24.cc = new('GOHyperGParams',
                   geneIds=as.character(table24$ID),
                   universeGeneIds=mapped_genes,
                   ontology=c('CC'),
                   pvalueCutoff=0.001,
                   conditional=F,
                   testDirection='over',
                   annotation="org.Mm.eg.db"
)
hgOver24.cc = hyperGTest(params24.cc)
result24.cc = summary(hgOver24.cc)
result24.cc$Adj.Pvalue = p.adjust(result24.cc$Pvalue, method = 'BH')
result24.cc = cbind(result24.cc[1:2], result24.cc["Adj.Pvalue"], result24.cc[3:7])
head(result24.cc)

# Biological process
params12.bp = new('GOHyperGParams',
                geneIds=as.character(table12$ID),
                universeGeneIds=mapped_genes,
                ontology=c('BP'),
                pvalueCutoff=0.001,
                conditional=F,
                testDirection='over',
                annotation="org.Mm.eg.db"
)
hgOver12.bp = hyperGTest(params12.bp)
result12.bp = summary(hgOver12.bp)
result12.bp$Adj.Pvalue = p.adjust(result12.bp$Pvalue, method = 'BH')
result12.bp = cbind(result12.bp[1:2], result12.bp["Adj.Pvalue"], result12.bp[3:7])
head(result12.bp)

params13.bp = new('GOHyperGParams',
                geneIds=as.character(table13$ID),
                universeGeneIds=mapped_genes,
                ontology=c('BP'),
                pvalueCutoff=0.001,
                conditional=F,
                testDirection='over',
                annotation="org.Mm.eg.db"
)
hgOver13.bp = hyperGTest(params13.bp)
result13.bp = summary(hgOver13.bp)
result13.bp$Adj.Pvalue = p.adjust(result13.bp$Pvalue, method = 'BH')
result13.bp = cbind(result13.bp[1:2], result13.bp["Adj.Pvalue"], result13.bp[3:7])
head(result13.bp)

params23.bp = new('GOHyperGParams',
                geneIds=as.character(table23$ID),
                universeGeneIds=mapped_genes,
                ontology=c('BP'),
                pvalueCutoff=0.001,
                conditional=F,
                testDirection='over',
                annotation="org.Mm.eg.db"
)
hgOver23.bp = hyperGTest(params23.bp)
result23.bp = summary(hgOver23.bp)
result23.bp$Adj.Pvalue = p.adjust(result23.bp$Pvalue, method = 'BH')
result23.bp = cbind(result23.bp[1:2], result23.bp["Adj.Pvalue"], result23.bp[3:7])
head(result23.bp)

params34.bp = new('GOHyperGParams',
                geneIds=as.character(table34$ID),
                universeGeneIds=mapped_genes,
                ontology=c('BP'),
                pvalueCutoff=0.001,
                conditional=F,
                testDirection='over',
                annotation="org.Mm.eg.db"
)
hgOver34.bp = hyperGTest(params34.bp)
result34.bp = summary(hgOver34.bp)
result34.bp$Adj.Pvalue = p.adjust(result34.bp$Pvalue, method = 'BH')
result34.bp = cbind(result34.bp[1:2], result34.bp["Adj.Pvalue"], result34.bp[3:7])
head(result34.bp)

params14.bp = new('GOHyperGParams',
                geneIds=as.character(table14$ID),
                universeGeneIds=mapped_genes,
                ontology=c('BP'),
                pvalueCutoff=0.001,
                conditional=F,
                testDirection='over',
                annotation="org.Mm.eg.db"
)
hgOver14.bp = hyperGTest(params14.bp)
result14.bp = summary(hgOver14.bp)
result14.bp$Adj.Pvalue = p.adjust(result14.bp$Pvalue, method = 'BH')
result14.bp = cbind(result14.bp[1:2], result14.bp["Adj.Pvalue"], result14.bp[3:7])
head(result14.bp)

params24.bp = new('GOHyperGParams',
                geneIds=as.character(table24$ID),
                universeGeneIds=mapped_genes,
                ontology=c('BP'),
                pvalueCutoff=0.001,
                conditional=F,
                testDirection='over',
                annotation="org.Mm.eg.db"
)
hgOver24.bp = hyperGTest(params24.bp)
result24.bp = summary(hgOver24.bp)
result24.bp$Adj.Pvalue = p.adjust(result24.bp$Pvalue, method = 'BH')
result24.bp = cbind(result24.bp[1:2], result24.bp["Adj.Pvalue"], result24.bp[3:7])
head(result24.bp)


# Molecular function
params12.mf = new('GOHyperGParams',
                   geneIds=as.character(table12$ID),
                   universeGeneIds=mapped_genes,
                   ontology=c('MF'),
                   pvalueCutoff=0.001,
                   conditional=F,
                   testDirection='over',
                   annotation="org.Mm.eg.db"
)
hgOver12.mf = hyperGTest(params12.mf)
result12.mf = summary(hgOver12.mf)
result12.mf$Adj.Pvalue = p.adjust(result12.mf$Pvalue, method = 'BH')
result12.mf = cbind(result12.mf[1:2], result12.mf["Adj.Pvalue"], result12.mf[3:7])
head(result12.mf)

params13.mf = new('GOHyperGParams',
                   geneIds=as.character(table13$ID),
                   universeGeneIds=mapped_genes,
                   ontology=c('MF'),
                   pvalueCutoff=0.001,
                   conditional=F,
                   testDirection='over',
                   annotation="org.Mm.eg.db"
)
hgOver13.mf = hyperGTest(params13.mf)
result13.mf = summary(hgOver13.mf)
result13.mf$Adj.Pvalue = p.adjust(result13.mf$Pvalue, method = 'BH')
result13.mf = cbind(result13.mf[1:2], result13.mf["Adj.Pvalue"], result13.mf[3:7])
head(result13.mf)

params23.mf = new('GOHyperGParams',
                   geneIds=as.character(table23$ID),
                   universeGeneIds=mapped_genes,
                   ontology=c('MF'),
                   pvalueCutoff=0.001,
                   conditional=F,
                   testDirection='over',
                   annotation="org.Mm.eg.db"
)
hgOver23.mf = hyperGTest(params23.mf)
result23.mf = summary(hgOver23.mf)
result23.mf$Adj.Pvalue = p.adjust(result23.mf$Pvalue, method = 'BH')
result23.mf = cbind(result23.mf[1:2], result23.mf["Adj.Pvalue"], result23.mf[3:7])
head(result23.mf)

params34.mf = new('GOHyperGParams',
                   geneIds=as.character(table34$ID),
                   universeGeneIds=mapped_genes,
                   ontology=c('MF'),
                   pvalueCutoff=0.001,
                   conditional=F,
                   testDirection='over',
                   annotation="org.Mm.eg.db"
)
hgOver34.mf = hyperGTest(params34.mf)
result34.mf = summary(hgOver34.mf)
result34.mf$Adj.Pvalue = p.adjust(result34.mf$Pvalue, method = 'BH')
result34.mf = cbind(result34.mf[1:2], result34.mf["Adj.Pvalue"], result34.mf[3:7])
head(result34.mf)

params14.mf = new('GOHyperGParams',
                   geneIds=as.character(table14$ID),
                   universeGeneIds=mapped_genes,
                   ontology=c('MF'),
                   pvalueCutoff=0.001,
                   conditional=F,
                   testDirection='over',
                   annotation="org.Mm.eg.db"
)
hgOver14.mf = hyperGTest(params14.mf)
result14.mf = summary(hgOver14.mf)
result14.mf$Adj.Pvalue = p.adjust(result14.mf$Pvalue, method = 'BH')
result14.mf = cbind(result14.mf[1:2], result14.mf["Adj.Pvalue"], result14.mf[3:7])
head(result14.mf)

params24.mf = new('GOHyperGParams',
                   geneIds=as.character(table24$ID),
                   universeGeneIds=mapped_genes,
                   ontology=c('MF'),
                   pvalueCutoff=0.001,
                   conditional=F,
                   testDirection='over',
                   annotation="org.Mm.eg.db"
)
hgOver24.mf = hyperGTest(params24.mf)
result24.mf = summary(hgOver24.mf)
result24.mf$Adj.Pvalue = p.adjust(result24.mf$Pvalue, method = 'BH')
result24.mf = cbind(result24.mf[1:2], result24.mf["Adj.Pvalue"], result24.mf[3:7])
head(result24.mf)


# HTML doc of results
htmlReport(hgOver12.cc, file = paste(dir, sep = "", "GO enrichment analysis/CC_HF_SAMP6_vs_LF_SAMP6.html"))
htmlReport(hgOver13.cc, file = paste(dir, sep = "", "GO enrichment analysis/CC_HF_SAMP6_vs_HF_AKRJ.html"))
htmlReport(hgOver23.cc, file = paste(dir, sep = "", "GO enrichment analysis/CC_LF_SAMP6_vs_HF_AKRJ.html"))
htmlReport(hgOver34.cc, file = paste(dir, sep = "", "GO enrichment analysis/CC_HF_AKRJ_vs_LF_AKRJ.html"))
htmlReport(hgOver14.cc, file = paste(dir, sep = "", "GO enrichment analysis/CC_HF_SAMP6_vs_LF_AKRJ.html"))
htmlReport(hgOver24.cc, file = paste(dir, sep = "", "GO enrichment analysis/CC_LF_SAMP6_vs_LF_AKRJ.html"))

htmlReport(hgOver12.bp, file = paste(dir, sep = "", "GO enrichment analysis/BP_HF_SAMP6_vs_LF_SAMP6.html"))
htmlReport(hgOver13.bp, file = paste(dir, sep = "", "GO enrichment analysis/BP_HF_SAMP6_vs_HF_AKRJ.html"))
htmlReport(hgOver23.bp, file = paste(dir, sep = "", "GO enrichment analysis/BP_LF_SAMP6_vs_HF_AKRJ.html"))
htmlReport(hgOver34.bp, file = paste(dir, sep = "", "GO enrichment analysis/BP_HF_AKRJ_vs_LF_AKRJ.html"))
htmlReport(hgOver14.bp, file = paste(dir, sep = "", "GO enrichment analysis/BP_HF_SAMP6_vs_LF_AKRJ.html"))
htmlReport(hgOver24.bp, file = paste(dir, sep = "", "GO enrichment analysis/BP_LF_SAMP6_vs_LF_AKRJ.html"))

htmlReport(hgOver12.mf, file = paste(dir, sep = "", "GO enrichment analysis/CC_HF_SAMP6_vs_LF_SAMP6.html"))
htmlReport(hgOver13.mf, file = paste(dir, sep = "", "GO enrichment analysis/CC_HF_SAMP6_vs_HF_AKRJ.html"))
htmlReport(hgOver23.mf, file = paste(dir, sep = "", "GO enrichment analysis/CC_LF_SAMP6_vs_HF_AKRJ.html"))
htmlReport(hgOver34.mf, file = paste(dir, sep = "", "GO enrichment analysis/CC_HF_AKRJ_vs_LF_AKRJ.html"))
htmlReport(hgOver14.mf, file = paste(dir, sep = "", "GO enrichment analysis/CC_HF_SAMP6_vs_LF_AKRJ.html"))
htmlReport(hgOver24.mf, file = paste(dir, sep = "", "GO enrichment analysis/CC_LF_SAMP6_vs_LF_AKRJ.html"))


# txt table of results
write.table(result12.mf, paste(dir, sep = "", "GO enrichment analysis/MF_HF_SAMP6_LF_SAMP6.txt"))
write.table(result13.mf, paste(dir, sep = "", "GO enrichment analysis/MF_HF_SAMP6_HF_AKRJ.txt"))
write.table(result23.mf, paste(dir, sep = "", "GO enrichment analysis/MF_LF_SAMP6_HF_AKRJ.txt"))
write.table(result34.mf, paste(dir, sep = "", "GO enrichment analysis/MF_HF_AKRJ_LF_SAMP6.txt"))
write.table(result14.mf, paste(dir, sep = "", "GO enrichment analysis/MF_HF_SAMP6_LF_AKRJ.txt"))
write.table(result24.mf, paste(dir, sep = "", "GO enrichment analysis/MF_LF_SAMP6_LF_AKRJ.txt"))

write.table(result12.bp, paste(dir, sep = "", "GO enrichment analysis/BP_HF_SAMP6_LF_SAMP6.txt"))
write.table(result13.bp, paste(dir, sep = "", "GO enrichment analysis/BP_HF_SAMP6_HF_AKRJ.txt"))
write.table(result23.bp, paste(dir, sep = "", "GO enrichment analysis/BP_LF_SAMP6_HF_AKRJ.txt"))
write.table(result34.bp, paste(dir, sep = "", "GO enrichment analysis/BP_HF_AKRJ_LF_SAMP6.txt"))
write.table(result14.bp, paste(dir, sep = "", "GO enrichment analysis/BP_HF_SAMP6_LF_AKRJ.txt"))
write.table(result24.bp, paste(dir, sep = "", "GO enrichment analysis/BP_LF_SAMP6_LF_AKRJ.txt"))

write.table(result12.cc, paste(dir, sep = "", "GO enrichment analysis/CC_HF_SAMP6_LF_SAMP6.txt"))
write.table(result13.cc, paste(dir, sep = "", "GO enrichment analysis/CC_HF_SAMP6_HF_AKRJ.txt"))
write.table(result23.cc, paste(dir, sep = "", "GO enrichment analysis/CC_LF_SAMP6_HF_AKRJ.txt"))
write.table(result34.cc, paste(dir, sep = "", "GO enrichment analysis/CC_HF_AKRJ_LF_SAMP6.txt"))
write.table(result14.cc, paste(dir, sep = "", "GO enrichment analysis/CC_HF_SAMP6_LF_AKRJ.txt"))
write.table(result24.cc, paste(dir, sep = "", "GO enrichment analysis/CC_LF_SAMP6_LF_AKRJ.txt"))