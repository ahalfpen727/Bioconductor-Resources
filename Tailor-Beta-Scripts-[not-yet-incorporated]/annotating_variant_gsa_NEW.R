#  GSA for gene set enrinchment with a section from vcftools
# for variant annotation and ID conversion
library("GSEABase")
library(GSA)
library(refGenome)
vignette("readGTF")
vignette("ExpressionSetIntroduction")
vignette("refGenome")
###################################################
readGTF 
GSAsummaryTable(gsaRes, save=TRUE, file="gsaResTab.xls")
 nw <- networkPlot(gsaRes,class="non")
 nw$geneSets
 # An example usage:
myGsc <- loadGSC("myModel.xml")
#####################################################
# To load iTO977:
metMap <- system.file("extdata", "probe2metabolites_iTO977.txt.gz", 
                      package="piano")
# To load iIN800:
metMap <- system.file("extdata", "probe2metabolites_iIN800.txt.gz", 
                      package="piano")
# Convert into piano format:
metMap <- read.delim(metMap)
myGsc <- loadGSC(metMap)
gsaRes <- runGSA(myPval,myFC,gsc=myGsc,
                 geneSetStat="reporter",
                 signifMethod="nullDist", nPerm=1000,
                 gsSizeLim=c(5,100))
###################################################
library(piano)
library(VariantAnnotation)
library(cgdv17)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
sourceUrl(BSgenome)
biocLite(PolyPhen.Hsapiens.dbSNP131)

oxygen <- c("aerobic","anaerobic")
limitation <- c("Clim","Nlim")
mySetup <- cbind(oxygen[c(1,1,1,2,2,2,1,1,1,2,2,2)],
               limitation[c(1,1,1,1,1,1,2,2,2,2,2,2)])

# The rownames correspond to the CEL-file names (CAE1.CEL etc):
rownames(mySetup) <- c("CAE1","CAE2","CAE3","CAN1","CAN2","CAN3",
                     "NAE1","NAE2","NAE3","NAN1","NAN2","NAN3")
colnames(mySetup) <- c("oxygen","limitation")

# The final setup object can look like this:
mySetup
# Get path to example data and setup files:
dataPath <- system.file("extdata", package="piano")

# Load pre-normalized data:
myArrayData <- loadMAdata(datadir=dataPath, 
                          dataNorm="norm_data.txt.gz", 
                          platform="yeast2")
myArrayData
myArrayData$setup
# Check the annotation (top 10 rows):
myArrayData$annotation[1:10,]
runQC(myArrayData)
# To only run the PCA:
runQC(myArrayData, rnaDeg=FALSE, nuseRle=FALSE, hist=FALSE,
      boxplot=FALSE, pca=TRUE)
# Additionally, for the PCA you can specify other colors:
runQC(myArrayData, rnaDeg=FALSE, nuseRle=FALSE, hist=FALSE,
      boxplot=FALSE, pca=TRUE, colors=c("cyan","orange"))
extractFactors(myArrayData)


pfc <- diffExp(myArrayData, contrasts=c("aerobic_Clim - anaerobic_Clim", 
                                        "aerobic_Nlim - anaerobic_Nlim"))
# Sort genes in "aerobic_Clim - anaerobic_Clim" according to adjusted p-value:
ii <- sort(pfc$resTable[[1]]$adj.P.Val, index.return=TRUE)$ix
pfc$resTable[[1]][ii[1:5],]
# Sort genes in "aerobic_Nlim - anaerobic_Nlim" according to adjusted p-value:
ii <- sort(pfc$resTable[[2]]$adj.P.Val, index.return=TRUE)$ix
pfc$resTable[[2]][ii[1:5],]

# Get p-values from the aerobic_Clim vs anaerobic_Clim comparison:
myPval <- pfc$pValues["aerobic_Clim - anaerobic_Clim"]
# Display the first values and gene IDs:
head(myPval)
# Custom gene to gene set mapping:
genes2genesets <- cbind(paste("gene",c("A","A","A","B","B","C","C","C","D"),sep=""),
                        paste("set",c(1,2,3,1,3,2,3,4,4),sep=""))
genes2genesets
# Load into correct format:
myGsc <- loadGSC(genes2genesets)
# View summary:
myGsc
# View all gene sets:
myGsc$gsc

myStats <- c(-1.5,-0.5,1,2)
names(myStats) <- paste("gene",c("A","B","C","D"),sep="")
myStats
###################################################
gsaRes <- runGSA(myStats, gsc=myGsc)

gsaRes
names(gsaRes)
###################################################
#  GSA WITH GO AND BIOMART
#  ANNOTATION WITH VCF
###################################################
library("biomaRt")
# Select ensembl database and S. cerevisiae dataset:
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="scerevisiae_gene_ensembl",host="www.ensembl.org")
# Map Yeast 2.0 microarray probeset IDs to GO:
mapGO <- getBM(attributes=c('affy_yeast_2', 'name_1006'), 
               filters = 'affy_yeast_2', 
               values = rownames(myPval), 
               mart = ensembl)
# Remove blanks ("")
mapGO <- mapGO[mapGO[,2]!="",]
# Check the 10 first rows to see what we got:
mapGO[1:10,]
myGsc <- loadGSC(mapGO)

#GSAsummaryTable(gsaRes, save=TRUE, file="gsaResTab.xls")
#nw <- networkPlot(gsaRes,class="non")
## nw <- networkPlot(gsaRes,class="non")

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("mypackage")
library(VariantAnnotation)
library(cgdv17)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(PolyPhen.Hsapiens.dbSNP131)

library(VariantAnnotation)
library(cgdv17)
file <- system.file("vcf", "NA06985_17.vcf.gz", package = "cgdv17")

#Examine header data in a vcf file 
hdr <- scanVcfHeader(file)
info(hdr) 
geno(hdr) 
## DataFrame with 12 rows and 3 columns
# Variants in the VCF have been aligned NCBI genome build GRCh37:
meta(hdr)$META
#############################################################
#   Convert gene symbols to gene ids
#
#Use the org.Hs.eg.db package to convert gene symbols to gene ids.
##########################################################3333

## get entrez ids from gene symbols
library(org.Hs.eg.db)
genesym <- c("TRPV1", "TRPV2", "TRPV3")
geneid <- select(org.Hs.eg.db, keys="ENSEMBL", keytype="SYMBOL",
                 columns="ENTREZID")
geneid 

#Create gene ranges
#We use the hg19 known gene track from UCSC
#to identify the TRPV gene ranges. 
#ranges will be used to extract variants from the VCF file.
#Load the annotation package.

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb

# VCF file was aligned to a genome from NCBI 
#known gene track was from UCSC. 
#NCBI AND UCSC have different naming conventions for chromosomes.
# In order to use these two pieces of data in
#a matching or overlap operation the chromosome names
#(also called sesqlevels) need to match. We will modify the txdb 
# to match the VCF file.

txdb <- renameSeqlevels(txdb, gsub("chr", "", seqlevels(txdb)))
txdb <- keepSeqlevels(txdb, "17")
#Create a list of transcripts by gene:
txbygene = transcriptsBy(txdb, "gene")
#Create the gene ranges for the TRPV genes
gnrng <- unlist(range(txbygene[geneid$ENTREZID]), use.names=FALSE)
names(gnrng) <- geneid$SYMBOL
# Extract variant subsets
#A ScanVcfParam object is used to retrieve data subsets. 
#This object can specify genomic coordinates (ranges)
#or individual VCF elements. Extractions of ranges (vs fields) 
#requires a tabix index. See ?indexTabix for details.

param <- ScanVcfParam(which = gnrng, info = "DP",
                      geno = c("GT", "cPd"))
param
## Extract the TRPV ranges from the VCF file 
vcf <- readVcf(file, "hg19", param)
## Inspect the VCF object with the 'fixed', 'info' and 'geno' accessors
vcf
head(fixed(vcf))
#  Variant location in the gene model
# The locateVariants function
# identifies where a variant falls with respect to gene structure
# e.g., exon, utr, splice site, etc.
# We use the gene model from the TxDb.Hsapiens.UCSC.hg19.knownGene  
# package loaded eariler.

## Use the 'region' argument to define the region
## of interest. See ?locateVariants for details.
cds <- locateVariants(vcf, txdb, CodingVariants())
five <- locateVariants(vcf, txdb, FiveUTRVariants())
splice <- locateVariants(vcf, txdb, SpliceSiteVariants())
intron <- locateVariants(vcf, txdb, IntronVariants())

all <- locateVariants(vcf, txdb, AllVariants())

# Each row in cds represents a variant-transcript match 
# so multiple rows per variant are possible. 
#If interested in gene-centric questions the data can be 
# summarized by gene regardless of transcript.
## variants match more than one gene?
table(sapply(split(mcols(all)$GENEID, mcols(all)$QUERYID), 
      function(x) length(unique(x)) > 1))

## Summarize the number of variants by gene:
idx <- sapply(split(mcols(all)$QUERYID, mcols(all)$GENEID), unique)
sapply(idx, length)
## Summarize variant location by gene:
sapply(names(idx), 
    function(nm) {
        d <- all[mcols(all)$GENEID %in% nm, c("QUERYID", "LOCATION")]
        table(mcols(d)$LOCATION[duplicated(d) == FALSE])
    })


# Amino acid coding changes in non-synonymous variants
# predictCoding function
# copute Amino acid coding for non-synonymous variants
# BSgenome.Hsapiens.UCSC.hg19 package is used as the source 
#of the reference alleles. 
# Variant alleles are provided by the user
library(BSgenome.Hsapiens.UCSC.hg19)
seqlevelsStyle(vcf) <- "UCSC"
seqlevelsStyle(txdb) <- "UCSC"
aa <- predictCoding(vcf, txdb, Hsapiens)

## Did any variants match more than one gene?
table(sapply(split(mcols(aa)$GENEID, mcols(aa)$QUERYID), 
        function(x) length(unique(x)) > 1))
## Summarize the number of variants by gene:
idx <- sapply(split(mcols(aa)$QUERYID, mcols(aa)$GENEID, drop=TRUE), unique)
sapply(idx, length)
## Summarize variant consequence by gene:
sapply(names(idx), 
       function(nm) {
           d <- aa[mcols(aa)$GENEID %in% nm, c("QUERYID","CONSEQUENCE")]
           table(mcols(d)$CONSEQUENCE[duplicated(d) == FALSE])
       })
# Annotating with the ensemblVEP package 
# Load ensemblVEP:
# The ‘file’ argument to ensemblVEP must be a vcf file on disk. 
# write VCF object with TRPV variants + submit to ensemblVEP.

library(ensemblVEP)
dest <- tempfile()
writeVcf(vcf, dest)
# Call ensemblVEP with the file containing only the TRPV variants 
# and the custom VEPParam object:
gr <- ensemblVEP(file = dest)
head(gr, 3)

# Exploring Package Content

#Packages have extensive help pages, 
#and include vignettes highlighting common use cases. 
#The help pages and vignettes are available from within R. 
#After loading a package, use syntax like
#to obtain an overview of help on the VariantAnnotation package, 
#and the predictCoding function. View the package vignette with
help(package="VariantAnnotation")
?predictCoding
browseVignettes(package="VariantAnnotation")

#To view vignettes providing a more comprehensive introduction 
help.start()

###################################################
### code chunk number 23: piano-vignette.Rnw:439-442
###################################################
tmp<-c(100,95,5)
names(tmp)<-c("Genes (tot)","Genes (up)","Genes (down)")
tmp
tmp<-c(0.0001,0.0001,0.9999,0.0001,0.005)
names(tmp) <- c("p (non.dir)","p (dist.dir.up)","p (dist.dir.dn)","p (mix.dir.up)","p (mix.dir.dn)")
tmp
# Get the p-values for the contrast aerobic_Clim - anaerobic_Clim
myPval <- pfc$pValues["aerobic_Clim - anaerobic_Clim"]
head(myPval)
# Get the fold-changes for the contrast aerobic_Clim - anaerobic_Clim
myFC <- pfc$foldChanges["aerobic_Clim - anaerobic_Clim"]
head(myFC)

gsaRes

nw$geneSets
gsaResTab <- GSAsummaryTable(gsaRes)
# Which columns contain p-values:
grep("p \\(",colnames(gsaResTab),value=T)
grep("p \\(",colnames(gsaResTab))
ii <- which(gsaResTab[,10]<0.0001)
gsaResTab$Name[ii]
# Get minimum p-value for each gene set:
minPval <- apply(gsaResTab[,c(4,7,10,14,18)],1,min,na.rm=TRUE)
# Select significant gene sets:
ii <- which(minPval<0.0001)
gsaResTabSign <- gsaResTab[ii,c(1,4,7,10,14,18)]
# Look at the first 10 gene sets:
gsaResTabSign[1:10,]


nw<- networkPlot(gsaRes, class="distinct", direction="both",
                  significance=0.005, label="numbers")

###  (eval = FALSE)
###################################################
## GSAsummaryTable(gsaRes, save=TRUE, file="gsaResTab.xls")
## nw <- networkPlot(gsaRes,class="non")
## nw$geneSets
## # An example usage:
## myGsc <- loadGSC("myModel.xml")
##
###################################################
# To load iTO977:
metMap <- system.file("extdata", "probe2metabolites_iTO977.txt.gz", 
                      package="piano")
# To load iIN800:
metMap <- system.file("extdata", "probe2metabolites_iIN800.txt.gz", 
                      package="piano")
# Convert into piano format:
metMap <- read.delim(metMap)
myGsc <- loadGSC(metMap)
gsaRes <- runGSA(myPval,myFC,gsc=myGsc,
                 geneSetStat="reporter",
                 signifMethod="nullDist", nPerm=1000,
                 gsSizeLim=c(5,100))
###################################################


###################################################
### code chunk number 39: piano-vignette.Rnw:666-668 (eval = FALSE)
###################################################
## myTval <- pfc$resTable[["aerobic_Clim - anaerobic_Clim"]]$t
## names(myTval) <- pfc$resTable[["aerobic_Clim - anaerobic_Clim"]]$ProbesetID
###################################################
### code chunk number 40: piano-vignette.Rnw:671-686 (eval = FALSE)
###################################################
## myGsc <- loadGSC(mapGO)
## gsaRes1 <- runGSA(myTval,geneSetStat="mean",gsc=myGsc,
##                   nPerm=1000,gsSizeLim=c(10,800))
## gsaRes2 <- runGSA(myTval,geneSetStat="median",gsc=myGsc,
##                   nPerm=1000,gsSizeLim=c(10,800))
## gsaRes3 <- runGSA(myTval,geneSetStat="sum",gsc=myGsc,
##                   nPerm=1000,gsSizeLim=c(10,800))
## gsaRes4 <- runGSA(myTval,geneSetStat="maxmean",gsc=myGsc,
##                   nPerm=1000,gsSizeLim=c(10,800))
## gsaRes5 <- runGSA(myPval,myFC,geneSetStat="fisher",gsc=myGsc,
##                   nPerm=1000,gsSizeLim=c(10,800))
## gsaRes6 <- runGSA(myPval,myFC,geneSetStat="stouffer",gsc=myGsc,
##                   nPerm=1000,gsSizeLim=c(10,800))
## gsaRes7 <- runGSA(myPval,myFC,geneSetStat="tailStrength",gsc=myGsc,
##                   nPerm=1000,gsSizeLim=c(10,800))
###################################################
## resList <- list(gsaRes1,gsaRes2,gsaRes3,gsaRes4,gsaRes5,gsaRes6,gsaRes7)
## names(resList) <- c("mean","median","sum","maxmean","fisher",
##                     "stouffer","tailStrength")
## ch <- consensusHeatmap(resList,cutoff=30,method="mean")
## ch$pMat
ensembl_us_west = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org")
head(listDatasets(ensembl_us_west))

library(biomaRt)
listEnsembl()
listEnsembl("GRCh=37")
listEnsembl(version=78)
 listMarts(host="uswest.ensembl.org")
ensembl = useEnsembl(biomart="ensembl")
head(listDatasets(ensembl))
grch37 = useEnsembl(biomart="ensembl",GRCh=37)
listDatasets(grch37)[31:35,]
 
ensembl78 = useEnsembl(biomart="ensembl",version=78)
listDatasets(ensembl78)[31:35,]
ensembl = useEnsembl(biomart="ensembl",
                     dataset="hsapiens_gene_ensembl")
ensembl = useEnsembl(biomart="ensembl", 
                     dataset="hsapiens_gene_ensembl", GRCh=37)
ensembl = useEnsembl(biomart="ensembl", 
                     dataset="hsapiens_gene_ensembl", version=78)

head(listFilters(ensembl))
head(listAttributes(ensembl))
chr1_genes <- getBM(attributes=c('ensembl_gene_id',
                    'ensembl_transcript_id','hgnc_symbol',
                    'chromosome_name','start_position',
                    'end_position'),filters ='chromosome_name', 
                    values ="1", mart = ensembl)
head(chr1_gene)
hgnc_swissprot <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol','hgnc_id','uniprot_swissprot'),filters = 'ensembl_gene_id', values = 'ENSG00000139618', mart = ensembl)
hgnc_swissprot
 
 
 
 
 
 
 
 