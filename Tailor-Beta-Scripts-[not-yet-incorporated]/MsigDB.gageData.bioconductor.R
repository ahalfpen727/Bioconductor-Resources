# using one of the annotation packges
library(AnnotationDbi)
library(org.Hs.eg.db) # or, e.g. Homo.sapiens
columns(org.Hs.eg.db)
k<-keytypes(org.Hs.eg.db)
k
sym<-keys(org.Hs.eg.db, keytype="SYMBOL")
egid<-keys(org.Hs.eg.db, keytype="ENTREZID")
refseq<-keys(org.Hs.eg.db, keytype="REFSEQ")
go.all.ids<-keys(org.Hs.eg.db, keytype="GOALL")
go.ids<-keys(org.Hs.eg.db, keytype="GO")
ontid<-keys(org.Hs.eg.db, keytype="ONTOLOGYALL")
ont<-keys(org.Hs.eg.db, keytype="ONTOLOGY")
pathid<-keys(org.Hs.eg.db, keytype="PATH")
ucscid<-keys(org.Hs.eg.db, keytype="UCSCKG")
# returns a named character vector, see ?mapIds for multiVals options
sym.res <- mapIds(org.Hs.eg.db, keys=k, column="SYMBOL",
                  keytype="ENTREZID")
# generates warning for 1:many mappings
res <- select(org.Hs.eg.db, keys=k, columns=c("ENTREZID",
                  "ENSEMBL","SYMBOL"), keytype="ENTREZID")

# map from one annotation to another using biomart
library(biomaRt)
listMarts()
listMarts(archive = TRUE)
entrez=c("673","837")
goids = getBM(attributes = c('entrezgene', 'go_id'), 
              filters = 'entrezgene', 
              values = entrez, 
              mart = ensembl)
head(goids)
m<-useMart("ensembl",, dataset = "hsapiens_gene_ensembl")
listDatasets(m)
m<-useMart("ensembl", dataset = "hsapiens_gene_ensembl")
map <- getBM(mart = m,
  attributes = c("ensembl_gene_id", "entrezgene"),
  filters = "ensembl_gene_id", 
  values = some.ensembl.genes)

filters = listFilters(m)
attributes = listAttributes(m)
goids = getBM(attributes = c('hgnc_symbol', 'go_id'), 
              filters = 'hgnc_symbol', 
              values = entrez, 
              mart = m)
head(goids)

goids = getBM(attributes = c('entrezgene', 'go_id'), 
              filters = 'entrezgene', 
              values = entrez, 
              mart = m)
head(goids)


goids = getBM(attributes = c('refseq_mrna', 'go_id'), 
              filters = 'refseq_mrna', 
              values = entrez, 
              mart = m)
head(goids)

snpmart = useMart(biomart = "ENSEMBL_MART_SNP", dataset="hsapiens_snp")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")


library(GenomicRanges)
z <- GRanges("chr1",IRanges(1000001,1001000),strand="+")
start(z)
end(z)
width(z)
strand(z)
mcols(z) # the 'metadata columns', any information stored alongside each range
ranges(z) # gives the IRanges
seqnames(z) # the chromosomes for each ranges
seqlevels(z) # the possible chromosomes
seqlengths(z) # the lengths for each chromosome
library(BSgenome.Hsapiens.UCSC.hg19)
dnastringset <- getSeq(Hsapiens, granges) # returns a DNAStringSet
# also Views() for Bioconductor >= 3.1

library(Biostrings)
dnastringset <- readDNAStringSet("transcripts.fa")

substr(dnastringset, 1, 10) # to character string
subseq(dnastringset, 1, 10) # returns DNAStringSet
Views(dnastringset, 1, 10) # lightweight views into object
complement(dnastringset)
reverseComplement(dnastringset)
matchPattern("ACGTT", dnastring) # also countPattern, also works on Hsapiens/genome
vmatchPattern("ACGTT", dnastringset) # also vcountPattern
letterFrequecy(dnastringset, "CG") # how many C's or G's
# also letterFrequencyInSlidingView
alphabetFrequency(dnastringset, as.prob=TRUE)
# also oligonucleotideFrequency, dinucleotideFrequency, trinucleotideFrequency
# transcribe/translate for imitating biological processes

Sequencing data

Rsamtools scanBam returns lists of raw values from BAM files

library(Rsamtools)
which <- GRanges("chr1",IRanges(1000001,1001000))
what <- c("rname","strand","pos","qwidth","seq")
param <- ScanBamParam(which=which, what=what)
# for more BamFile functions/details see ?BamFile
# yieldSize for chunk-wise access
bamfile <- BamFile("/path/to/file.bam")
reads <- scanBam(bamfile, param=param)
res <- countBam(bamfile, param=param) 
# for more sophisticated counting modes
# see summarizeOverlaps() below

# quickly check chromosome names
seqinfo(BamFile("/path/to/file.bam"))

# DNAStringSet is defined in the Biostrings package
# see the Biostrings Quick Overview PDF
dnastringset <- scanFa(fastaFile, param=granges)
library(GenomicAlignments)
ga <- readGAlignments(bamfile) # single-end
ga <- readGAlignmentPairs(bamfile) # paired-end

# get a transcript database, which stores exon, trancript, and gene information
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# or build a txdb from GTF file (e.g. downloadable from Ensembl FTP site)
txdb <- makeTranscriptDbFromGFF("file.GTF", format="gtf")

# or build a txdb from Biomart (however, not as easy to reproduce later)
txdb <- makeTranscriptDbFromBiomart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# in Bioconductor >= 3.1, also makeTxDbFromGRanges

# saving and loading
saveDb(txdb, file="txdb.sqlite")
loadDb("txdb.sqlite")

# extracting information from txdb
g <- genes(txdb) # GRanges, just start to end, no exon/intron information
tx <- transcripts(txdb) # GRanges, similar to genes()
e <- exons(txdb) # GRanges for each exon
ebg <- exonsBy(txdb, by="gene") # exons grouped in a GRangesList by gene
ebt <- exonsBy(txdb, by="tx") # similar but by transcript

# then get the transcript sequence
txSeq <- extractTranscriptSeqs(Hsapiens, ebt)

Summarizing information across ranges and experiments

The SummarizedExperiment is a storage class for high-dimensional information tied to the same GRanges or GRangesList across experiments (e.g., read counts in exons for each gene).

library(GenomicAlignments)
fls <- list.files(pattern="*.bam$")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
ebg <- exonsBy(txdb, by="gene")
# see yieldSize argument for restricting memory
bf <- BamFileList(fls)
library(BiocParallel)
register(MulticoreParam(4))
# lots of options in the man page
# singleEnd, ignore.strand, inter.features, fragments, etc.
se <- summarizeOverlaps(ebg, bf)

# operations on SummarizedExperiment
assay(se) # the counts from summarizeOverlaps
colData(se)
rowRanges(se)
library(Rsubread)
res <- featureCounts(files, annot.ext="annotation.gtf",
  isGTFAnnotationFile=TRUE,
  GTF.featureType="exon",
  GTF.attrType="gene_id")
res$counts


library(EGSEAdata)
egsea.data("mouse")

## ----collectionlist------------------------------------------------------
info = egsea.data("mouse", returnInfo = TRUE)
names(info)
info$msigdb$info$collections

## ----loadegsea, message=FALSE, warning=FALSE-----------------------------
library(EGSEA)

## ----indexing------------------------------------------------------------
gs.annots = buildIdx(entrezIDs=v$genes$ENTREZID, species="mouse", 
           msigdb.gsets=c("c2", "c5"), go.part = TRUE)
names(gs.annots)

## ----exploresets---------------------------------------------------------
class(gs.annots$c2)
summary(gs.annots$c2)
show(gs.annots$c2)
s = getSetByName(gs.annots$c2, "SMID_BREAST_CANCER_LUMINAL_A_DN")
class(s)
names(s)
names(s$SMID_BREAST_CANCER_LUMINAL_A_DN)

## ----indexclass----------------------------------------------------------
slotNames(gs.annots$c2)

## ----symbolmap-----------------------------------------------------------
colnames(v$genes)
symbolsMap = v$genes[, c(1, 2)]
colnames(symbolsMap) = c("FeatureID", "Symbols")
symbolsMap[, "Symbols"] = as.character(symbolsMap[, "Symbols"])

## ----base----------------------------------------------------------------
egsea.base()

## ----selectbasemethods---------------------------------------------------
baseMethods = egsea.base()[-2]
baseMethods

## ----combine-------------------------------------------------------------
egsea.combine()

## ----sort----------------------------------------------------------------
egsea.sort()

## ----egseatest-----------------------------------------------------------
gsa = egsea(voom.results=v, contrasts=contr.matrix,  
         gs.annots=gs.annots, symbolsMap=symbolsMap,
         baseGSEAs=baseMethods, sort.by="med.rank",
         num.threads = 8, report = FALSE)

## ----showegsea-----------------------------------------------------------
show(gsa)

## ----summariseegsea------------------------------------------------------
summary(gsa)

## ----topsets-------------------------------------------------------------
topSets(gsa, gs.label="c2", contrast = "comparison", names.only=TRUE)

## ----topsetslim----------------------------------------------------------
t = topSets(gsa, contrast = "comparison",
             names.only=FALSE, number = Inf, verbose = FALSE)
t[grep("LIM_", rownames(t)), c("p.adj", "Rank", "med.rank", "vote.rank")]

## ----topsets2------------------------------------------------------------
topSets(gsa, gs.label="kegg", contrast="BasalvsLP", sort.by="med.rank")
topSets(gsa, gs.label="kegg", contrast="comparison", sort.by="med.rank")

## ----heatmaps------------------------------------------------------------
plotHeatmap(gsa, gene.set="LIM_MAMMARY_STEM_CELL_UP", gs.label="c2",
         contrast = "comparison", file.name = "hm_cmp_LIM_MAMMARY_STEM_CELL_UP", format="png")
plotHeatmap(gsa, gene.set="LIM_MAMMARY_STEM_CELL_DN", gs.label="c2",
         contrast = "comparison", file.name = "hm_cmp_LIM_MAMMARY_STEM_CELL_DN", format="png")

## ----pathwayplot1, eval=FALSE--------------------------------------------
#  plotPathway(gsa, gene.set = "Vascular smooth muscle contraction",
#               contrast = "BasalvsLP", gs.label = "kegg",
#               file.name = "Vascular_smooth_muscle_contraction")

## ----pathwayplot2, eval=FALSE--------------------------------------------
#  plotPathway(gsa, gene.set = "Vascular smooth muscle contraction",
#               contrast = "comparison", gs.label = "kegg",
#               file.name = "Vascular_smooth_muscle_contraction_cmp")

## ----mdsplot-------------------------------------------------------------
plotMethods(gsa, gs.label = "c2", contrast = "BasalvsLP", 
         file.name = "mds_c2_BasalvsLP", format="png")
plotMethods(gsa, gs.label = "c5BP", contrast = "BasalvsLP", 
         file.name = "mds_c5_BasalvsLP", format="png")

## ----keggsummaryplot1----------------------------------------------------
plotSummary(gsa, gs.label = 3, contrast = 3, 
         file.name = "summary_kegg_LPvsML", format="png")

## ----c2summaryplot2------------------------------------------------------
plotSummary(gsa, gs.label = 1, contrast = 3, 
         file.name = "summary_c2_LPvsML", 
         x.axis = "med.rank", format="png")

## ----c2summaryplot3------------------------------------------------------
plotSummary(gsa, gs.label = 1, contrast = 3, 
         file.name = "summary_sig_c2_LPvsML", 
         x.axis = "med.rank", x.cutoff=300, format="png")

## ----summaryplotkegg1and2------------------------------------------------
plotSummary(gsa, gs.label = "kegg", contrast = c(1,2), 
         file.name = "summary_kegg_1vs2", format="png")

## ----gographs------------------------------------------------------------
plotGOGraph(gsa, gs.label="c5BP", contrast = 1, file.name="BasalvsLP-c5BP-top-", format="png")
plotGOGraph(gsa, gs.label="c5CC", contrast = 1, file.name="BasalvsLP-c5CC-top-", format="png")

## ----summarybarplot------------------------------------------------------
plotBars(gsa, gs.label = "c2", contrast = "comparison", file.name="comparison-c2-bars", format="png")

## ----summaryheatmap------------------------------------------------------
plotSummaryHeatmap(gsa, gs.label="c2", hm.vals = "avg.logfc.dir",
         file.name="summary_heatmaps_c2", format="png")
plotSummaryHeatmap(gsa, gs.label="kegg", hm.vals = "avg.logfc.dir",
         file.name="summary_heatmaps_kegg", format="png")

## ----toptable------------------------------------------------------------
t = limmaTopTable(gsa, contrast=1)
head(t)

## ----htmlreport, warning=FALSE-------------------------------------------
generateReport(gsa, number = 20, report.dir="./mam-rnaseq-egsea-report") 

## ----mareadidats---------------------------------------------------------
library(limma)
url = "http://bioinf.wehi.edu.au/EGSEA/arraydata.zip"
utils::download.file(url, destfile="arraydata.zip", mode="wb")
utils::unzip("arraydata.zip", exdir = ".")
targets = read.delim("targets.txt", header=TRUE, sep=" ")
data = read.idat(as.character(targets$File), 
                   bgxfile="GPL6887_MouseWG-6_V2_0_R0_11278593_A.bgx",
                   annotation=c("Entrez_Gene_ID","Symbol", "Chromosome"))
data$other$Detection = detectionPValues(data)
data$targets = targets
colnames(data) = targets$Sample

## ----manormalize---------------------------------------------------------
data = neqc(data)

## ----mafilter------------------------------------------------------------
table(targets$Celltype)
keep.exprs = rowSums(data$other$Detection<0.05)>=5
table(keep.exprs)
data = data[keep.exprs,]
dim(data)

## ----malinearmodel-------------------------------------------------------
head(data$genes)
sum(is.na(data$genes$Entrez_Gene_ID))
data1 = data[!is.na(data$genes$Entrez_Gene_ID), ]
dim(data1)
expr = data1$E
group = as.factor(data1$targets$Celltype)
probe.annot = data1$genes[, 2:4]
head(probe.annot)
head(data1$targets)
experiment = as.character(data1$targets$Experiment)
design = model.matrix(~0 + group + experiment)
colnames(design) = gsub("group", "", colnames(design))
design
contr.matrix = makeContrasts(
         BasalvsLP = Basal-LP,
         BasalvsML = Basal-ML,
         LPvsML = LP-ML,
         levels = colnames(design))
contr.matrix

## ----maindex-------------------------------------------------------------
library(EGSEA)
library(EGSEAdata)
gs.annots = buildIdx(entrezIDs=unique(probe.annot[, 2]), 
             species="mouse", 
             msigdb.gsets=c("c2", "c5"), go.part = TRUE)
names(gs.annots)

## ----maegsea-------------------------------------------------------------
baseMethods = egsea.base()[-2]
baseMethods

gsam = egsea.ma(expr=expr, group=group, 
 		  	probe.annot = probe.annot,
 		  	design = design,
         contrasts=contr.matrix,  
         gs.annots=gs.annots,
         baseGSEAs=baseMethods, sort.by="med.rank",
         num.threads = 8, report = FALSE)

## ----mareport, warning=FALSE---------------------------------------------
generateReport(gsam, number = 20, report.dir="./mam-ma-egsea-report")

## ----matopsets-----------------------------------------------------------
topSets(gsam, gs.label="c2", contrast = "comparison", names.only=TRUE)

## ----softwareinfo--------------------------------------------------------
sessionInfo()

