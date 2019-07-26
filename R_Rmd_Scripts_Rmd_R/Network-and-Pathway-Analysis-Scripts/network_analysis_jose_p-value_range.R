source("http://bioconductor.org/biocLite.R")

# dir = "~/Documents/UMB/Riley_Lab/HiSeq_data/jose/"
## options(echo=TRUE)
## args = commandArgs(trailingOnly = T)
args = commandArgs(TRUE)

for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
}

p.value.start = 0
p.value.end = (1 * 10 ^ -100)

# read in tables
library(xlsx)

# increase jvm memory, otherwise it crashes half the time
# options(java.parameters = "-Xmx3072m")

## table = read.xlsx2(paste(inDir, filename, sep = ""), startRow = 2)
table = read.table(paste(inDir, filename, sep = ""), header=T)

# biocLite("BioNet")
# biocLite("DLBCL")
library(BioNet)
library(DLBCL)
data(interactome)

# reformat gene_id column to match node labels in the graph
reformatSym = function(table)
{
    table$gene_id = as.character(table$gene_id)
    table$ID = as.character(table$ID)
    for (node1 in interactome@nodes) {
        sym = interactome@nodeData@data[node1][[1]]$genegene_id
        if (!is.null(table["gene_id"][table["gene_id"] == sym])) {
            table["gene_id"][table["gene_id"] == sym] = node1
        }
    }
    return(table)
}

# reformat the 'gene_id' column so that genes in the data can be matched with nodes in the network
table = reformatSym(table)

# perform Bonferroni correction on p values, which are too low for BUM algorithm
table$q_value = as.numeric(levels(table$q_value))[table$q_value]
table$q_value = p.adjust(table$q_value, method = "bonferroni", n = nrow(table))

# remove all but the lowest p-value genes, so that the generated graph is readable
table = table[table$q_value > p.value.start & table$q_value < p.value.end,]
table
# get the subnetwork
subnet = subNetwork(table$gene_id, interactome)
subnet = rmSelfLoops(subnet)

# get the FDR from diff. expression data
# fdr = as.numeric(levels(table$FDR))[table$FDR]
# fdr = fdr[fdr > 0]
# names(fdr) = levels(table$ID)

# get p values from differential expression data
pval = table$q_value
pval = pval + (1 * 10 ^ -300) # add minimum acceptable p value to all values so that model can be created
names(pval) = table$gene_id

# fit Beta-uniform mixture model
bum = fitBumModel(pval, plot = F)

# score nodes based on model
scores = scoreNodes(subnet, bum, fdr = 0.001)

# create module
module = runFastHeinz(subnet, scores)
log2_fold_change) = as.numeric(levels(table$"log2(fold_change)"))[table$"log2(fold_change)"]
names(log2_fold_change)) = table$gene_id

# display subnetwork
## thresh = paste(strsplit(filename, "/")[[1]][1], "/", sep = "")
outfile = paste(outDir,"/", inFile, "_p", p.value.start, "-", p.value.end, ".png", sep = "")

png(file = outfile)
plotModule(module, scores = scores, diff.expr = log2(fold_change))
dev.off()
library("RWeka")
