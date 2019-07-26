source("http://bioconductor.org/biocLite.R")

dir = "~/Documents/UMB/Riley_Lab/HiSeq data/jose/"

# read in tables
library(xlsx)

# increase jvm memory, otherwise it crashes half the time
options(java.parameters = "-Xmx3072m")

table12 = read.xlsx2(paste(dir, "thresh100/edger_table_.0001_thresh100.xlsx", sep = ""),
                     sheetIndex = 1)
table13 = read.xlsx2(paste(dir, "thresh100/edger_table_.0001_thresh100.xlsx", sep = ""),
                     sheetIndex = 2, colIndex = 2:10)
table23 = read.xlsx2(paste(dir, "thresh100/edger_table_.0001_thresh100.xlsx", sep = ""),
                     sheetIndex = 3, colIndex = 2:10)
table32 = read.xlsx2(paste(dir, "thresh100/edger_table_.0001_thresh100.xlsx", sep = ""),
                     sheetIndex = 4, colIndex = 2:10)

biocLite("BioNet")
biocLite("DLBCL")
library(BioNet)
library(DLBCL)
data(interactome)

# reformat Symbol column to match node labels in the graph
reformatSym = function(table)
{
  table$Symbol = as.character(table$Symbol)
  table$ID = as.character(table$ID)
  for (node1 in interactome@nodes)
  {
    sym = interactome@nodeData@data[node1][[1]]$geneSymbol
    if (!is.null(table["Symbol"][table["Symbol"] == sym]))
    {
      table["Symbol"][table["Symbol"] == sym] = node1
    }
  }
  return(table)
}

# reformat the 'Symbol' column so that genes in the data can be matched with nodes in the network
table12 = reformatSym(table12)

# perform Bonferonni correction on p values, which are too low for BUM algorithm
table12$PValue = as.numeric(levels(table12$PValue))[table12$PValue]
table12$PValue = p.adjust(table12$PValue, method = "bonferroni", n = nrow(table12))
table12 = table12[table12$PValue < 1e-200,]

# get the subnetwork
subnet = subNetwork(table12$Symbol, interactome)
subnet = rmSelfLoops(subnet)
subnet

# get the FDR from diff. expression data
fdr12 = as.numeric(levels(table12$FDR))[table12$FDR]
fdr12 = fdr12[fdr12 > 0]
names(fdr12) = levels(table12$ID)

# get p values from differential expression data
# pval12 = as.numeric(levels(table12$PValue))[table12$PValue]
pval12 = table12$PValue
pval12 = pval12 + (1 * 10 ^ -300) # add minimum acceptable p value to all values so that model can be created
names(pval12) = table12$Symbol

# fit Beta-uniform mixture model
bum12 = fitBumModel(pval12, plot = F)

scores12 = scoreNodes(subnet, bum12, fdr = 0.001)
module12 = runFastHeinz(subnet, scores12)
logFC12 = as.numeric(levels(table12$logFC))[table12$logFC]
names(logFC12) = table12$Symbol
plotModule(module12, scores = scores12, diff.expr = logFC12)
