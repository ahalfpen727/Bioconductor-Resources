# perform GO analysis, write HTML output to file for each set of expression data and GO category
dir = "/Users/celizabeth/Documents/UMB/Riley_Lab/HiSeq data/jose/"

# read in tables
library(xlsx)

# increase jvm memory, otherwise it crashes half the time
options(java.parameters = "-Xmx3072m")

table12 = read.xlsx2(paste(dir, "thresh100/edger_table_.0001_thresh100.xlsx", sep = ""),
                     sheetIndex = 1, colIndex = 2:10)
table13 = read.xlsx2(paste(dir, "thresh100/edger_table_.0001_thresh100.xlsx", sep = ""),
                     sheetIndex = 2, colIndex = 2:10)
table23 = read.xlsx2(paste(dir, "thresh100/edger_table_.0001_thresh100.xlsx", sep = ""),
                     sheetIndex = 3, colIndex = 2:10)
table32 = read.xlsx2(paste(dir, "thresh100/edger_table_.0001_thresh100.xlsx", sep = ""),
                     sheetIndex = 4, colIndex = 2:10)

biocLite("org.Hs.eg.db")
library("org.Hs.eg.db")

biocLite("GOstats")
library(GOstats)

entrez_object = org.Hs.egGO
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
                   annotation="org.Hs.eg.db"
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
                   annotation="org.Hs.eg.db"
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
                   annotation="org.Hs.eg.db"
)
hgOver23.cc = hyperGTest(params23.cc)
result23.cc = summary(hgOver23.cc)
result23.cc$Adj.Pvalue = p.adjust(result23.cc$Pvalue, method = 'BH')
result23.cc = cbind(result23.cc[1:2], result23.cc["Adj.Pvalue"], result23.cc[3:7])
head(result23.cc)

# Biological process
params12.bp = new('GOHyperGParams',
                   geneIds=as.character(table12$ID),
                   universeGeneIds=mapped_genes,
                   ontology=c('BP'),
                   pvalueCutoff=0.001,
                   conditional=F,
                   testDirection='over',
                   annotation="org.Hs.eg.db"
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
                   annotation="org.Hs.eg.db"
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
                   annotation="org.Hs.eg.db"
)
hgOver23.bp = hyperGTest(params23.bp)
result23.bp = summary(hgOver23.bp)
result23.bp$Adj.Pvalue = p.adjust(result23.bp$Pvalue, method = 'BH')
result23.bp = cbind(result23.bp[1:2], result23.bp["Adj.Pvalue"], result23.bp[3:7])
head(result23.bp)

# Molecular function
params12.mf = new('GOHyperGParams',
                   geneIds=as.character(table12$ID),
                   universeGeneIds=mapped_genes,
                   ontology=c('MF'),
                   pvalueCutoff=0.001,
                   conditional=F,
                   testDirection='over',
                   annotation="org.Hs.eg.db"
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
                   annotation="org.Hs.eg.db"
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
                   annotation="org.Hs.eg.db"
)
hgOver23.mf = hyperGTest(params23.mf)
result23.mf = summary(hgOver23.mf)
result23.mf$Adj.Pvalue = p.adjust(result23.mf$Pvalue, method = 'BH')
result23.mf = cbind(result23.mf[1:2], result23.mf["Adj.Pvalue"], result23.mf[3:7])
head(result23.mf)

# HTML doc results
htmlReport(hgOver12.bp, paste(dir, sep = "", "thresh100/bp_untreated_cxcl12_.001_thresh100.html"))
htmlReport(hgOver13.bp, paste(dir, sep = "", "thresh100/bp_untreated_tgfb_.001_thresh100.html"))
htmlReport(hgOver23.bp, paste(dir, sep = "", "thresh100/bp_cxcl12_tgfb_.001_thresh100.html"))

htmlReport(hgOver12.cc, paste(dir, sep = "", "thresh100/cc_untreated_cxcl12_.001_thresh100.html"))
htmlReport(hgOver13.cc, paste(dir, sep = "", "thresh100/cc_untreated_tgfb_.001_thresh100.html"))
htmlReport(hgOver23.cc, paste(dir, sep = "", "thresh100/cc_cxcl12_tgfb_.001_thresh100.html"))

htmlReport(hgOver12.mf, paste(dir, sep = "", "thresh100/mf_untreated_cxcl12_.001_thresh100.html"))
htmlReport(hgOver13.mf, paste(dir, sep = "", "thresh100/mf_untreated_tgfb_.001_thresh100.html"))
htmlReport(hgOver23.mf, paste(dir, sep = "", "thresh100/mf_cxcl12_tgfb_.001_thresh100.html"))

# txt table of results
write.table(result12.bp, paste(dir, sep = "", "thresh100/go_bp_untreated_cxcl12_.001_thesh20.txt"))
write.table(result13.bp, paste(dir, sep = "", "thresh100/go_bp_untreated_tgfb_.001_thresh100.txt"))
write.table(result23.bp, paste(dir, sep = "", "thresh100/go_bp_cxcl12_tgfb_.001_thresh100.txt"))

write.table(result12.cc, paste(dir, sep = "", "thresh100/go_cc_untreated_cxcl12_.001_thesh20.txt"))
write.table(result13.cc, paste(dir, sep = "", "thresh100/go_cc_untreated_tgfb_.001_thresh100.txt"))
write.table(result23.cc, paste(dir, sep = "", "thresh100/go_cc_cxcl12_tgfb_.001_thresh100.txt"))

write.table(result12.mf, paste(dir, sep = "", "thresh100/go_mf_untreated_cxcl12_.001_thesh20.txt"))
write.table(result13.mf, paste(dir, sep = "", "thresh100/go_mf_untreated_tgfb_.001_thresh100.txt"))
write.table(result23.mf, paste(dir, sep = "", "thresh100/go_mf_cxcl12_tgfb_.001_thresh100.txt"))