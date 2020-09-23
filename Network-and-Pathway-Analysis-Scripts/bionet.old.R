source("http://bioconductor.org/biocLite.R")
options(error=traceback) # causes a traceback to appear if there is an error, and the traceback has the line number, prefixed by #

# dir = "~/Documents/UMB/Riley_Lab/HiSeq_data/jose/"
## options(echo=TRUE)
## args = commandArgs(trailingOnly = T)
args = commandArgs(TRUE)

for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
}

## p.value.start = 0
## p.value.end = 1 * (10 ^ -1)

# read in tables
## library(xlsx)

# increase jvm memory, otherwise it crashes half the time
# options(java.parameters = "-Xmx3072m")

## table = read.xlsx2(paste(inDir, inFile, sep = ""), startRow = 2)
tablePathName = paste(inDir, "/", inFile, sep = "")
print(paste("tablePathName=", tablePathName, sep=""))

cuffdiff.table = read.table(tablePathName, header=TRUE)

print("cuffdiff.table:");
cuffdiff.table[0:10, ]

# biocLite("BioNet")
# biocLite("DLBCL")
library(BioNet)
library(DLBCL)  # diffuse large B-cell lymphomas (DLBCL) data; also contains HPRD PPI

# Literature-curated human protein-protein interactions obtained from HPRD [11].
# Altogether the entire network used here comprises 9386 nodes and 36504 edges
print("Here 0")

data(interactome)

# reformat gene_id column to match node labels in the graph
reformatSym = function(table)
{
    table$gene_id = as.character(table$gene_id)
    ## table$ID = as.character(table$ID)
    for (node1 in interactome@nodes) {
        sym = interactome@nodeData@data[node1][[1]]$geneSymbol
        if (!is.null(table["gene_id"][table["gene_id"] == sym])) {
            table["gene_id"][table["gene_id"] == sym] = node1
            ## print(paste("sym=", sym, "; node1=",node1, sep=""))
        }
    }
    return(table)
}

# reformat the 'gene_id' column so that genes in the data can be matched with nodes in the network
cuffdiff.table = reformatSym(cuffdiff.table)

# perform Bonferroni correction on p values, which are too low for BUM algorithm
print("Here 1")

## cuffdiff.table$q_value = as.numeric(levels(cuffdiff.table$q_value))[cuffdiff.table$q_value]

print("Here 1.5")

cuffdiff.table[0:10, ]

cuffdiff.table$q_value = p.adjust(cuffdiff.table$q_value, method = "bonferroni", n = nrow(cuffdiff.table))

print("Here 2")

# remove all but the lowest p-value genes, so that the generated graph is readable
cuffdiff.table[0:10, ]

## cuffdiff.table = cuffdiff.table[(cuffdiff.table$q_value > p.value.start) & (cuffdiff.table$q_value < p.value.end), ]

## cuffdiff.table[0:10, ]

# get the subnetwork
print("Here 3")
subnet = subNetwork(cuffdiff.table$gene_id, interactome)
print("Here 4")
subnet = rmSelfLoops(subnet)
print("Here 5")

# get the FDR from diff. expression data
# fdr = as.numeric(levels(cuffdiff.table$FDR))[cuffdiff.table$FDR]
# fdr = fdr[fdr > 0]
# names(fdr) = levels(cuffdiff.table$ID)

# get p values from differential expression data
pval = cuffdiff.table$q_value
pval = pval + (1 * 10 ^ -300) # add minimum acceptable p value to all values so that model can be created
names(pval) = cuffdiff.table$gene_id

# fit Beta-uniform mixture model
bum = fitBumModel(pval, plot = FALSE)

# score nodes based on model
scores = scoreNodes(subnet, bum, fdr = 0.01)

# create module
module = runFastHeinz(subnet, scores)
## log2_fold_change = as.numeric(levels(cuffdiff.table$"log2(fold_change)"))[cuffdiff.table$"log2(fold_change)"]
## log2_fold_change = cuffdiff.table$"log2(fold_change)"
log2_fold_change = cuffdiff.table$log2.fold_change.
names(log2_fold_change) = cuffdiff.table$gene_id

# display subnetwork
## thresh = paste(strsplit(inFile, "/")[[1]][1], "/", sep = "")
## outFile = paste(outDir,"/", inFile, "_p", p.value.start, "-", p.value.end, ".png", sep = "")
outFile = paste(outDir,"/", inFile, ".png", sep = "")

png(file = outFile)
plotModule(module, scores = scores, diff.expr = log2_fold_change)
dev.off()
library("RWeka")
