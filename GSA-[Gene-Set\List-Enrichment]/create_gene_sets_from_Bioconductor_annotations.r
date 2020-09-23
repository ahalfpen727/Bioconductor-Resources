########################################################################
## create_gene_sets_from_Bioconductor_annotations.R
## Weil Lai
## wlai@alum.mit.edu
## August 13, 2007
##
## This R script lets users convert human, mouse, or rat annotations
## from Bioconductor annotation packages to a "G" list, which is
## required for analyzing gene sets in sigPathway. Although the script
## looks a bit convoluted, it takes less than one minute to run the script
## on a Pentium M 1.86 Mhz computer.
##
## The values for egidDir, biocAnnot, and perhaps other variables will
## need to be changed so that the script will run properly. Please
## make sure the Bioconductor annotation package of interest has already
## been installed in R before running this script.
##
## I have tested this script for hgu133a, hgu95av2, hs25kresogen, rat2302,
## and mgu74av2. It should work for other human, mouse, and rat
## Bioconductor annotations. More details regarding the "G" list format
## can be found in the sigPathway vignette.
########################################################################

## "Genesets_EntrezGeneIDs.RData" is an R workspace containing Entrez Gene
## IDs corresponding to gene sets from GO, KEGG, and other pathway
## databases listed in the sigPathway vignette PDF. It can be found at
## http://compbio.hms.harvard.edu/files/parklab/files/genesets_entrezgeneids_rdata.zip
##
## egidDir stands for the directory name where Genesets_EntrezGeneIDs.RData
## is located
 egidDir <- "W:/sigPathway development/pathway lists/egids"
 load(file.path(egidDir, "Genesets_EntrezGeneIDs.RData"))

## Give the Bioconductor annotation package name here
 biocAnnot <- "rat2302"

## get Entrez Gene IDs for each "probe set" (Affymetrix terminology)
## on the array
 library(biocAnnot, character.only = TRUE)
 xx <- as.list(get(paste(biocAnnot, "ENTREZID", sep = "")))
 xx <- xx[!is.na(xx)]
 xx <- unlist(xx)

xxUnique <- unique(xx)

yy <- vector("list", length(xxUnique))

for(i in 1:length(yy))
 yy[[i]] <- names(xx)[xx == xxUnique[i]]

## Match probe sets (by Entrez Gene IDs) to master gene list
 zz <- vector("list", length(G.EGIDs))

for(i in 1:length(zz)) {
 m <- match(G.EGIDs[[i]]$probes, xxUnique)
 zz[[i]] <- unlist(yy[m])

if(i %% 1000 == 0)
 cat("i = ", i, "\n", sep = "")
 }

## Disregard gene sets that did not get represented on the array
 idx <- which(sapply(zz, length) > 0)
 G.biocAnnot <- G.EGIDs[idx]

for(i in 1:length(idx))
 G.biocAnnot[[i]]$probes <- zz[[idx[i]]]

## Save mapped gene set to the directory specified in egidDir
 save(G.biocAnnot, file = file.path(egidDir, "G_biocAnnot.RData"))