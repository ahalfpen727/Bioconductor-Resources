# gtf_cleanup_script
#######################################
## load libraries
##########################################
library(ReactomePA);library(DOSE);library(SGSeq)
library(cummeRbund);library(plyr);library(clusterProfiler)
library(edgeR);library(limma);library(ggbio)
library(gage);library(pathview);library(EGSEA)
library(EGSEAdata);library(msigdbr);library(GSA)
library(STRINGdb);library(Rsubread);library(biomaRt)
library(EnrichmentBrowser);library(rtracklayer)
library(Rsamtools);library(msigdbr)
library(BioNet); library(DLBCL)
data(interactome)
browseVignettes("msigdbr")
#######################################
## initialize
##########################################

cuff <- readCufflinks(dir='/media/andrew/Seagate4TbExtHDD/umb_triley/urine1/cuffdiff_results_hg38_default/LUTS-over-CTRL/',genome="hg38",rebuild=T)
cuff
rna.seq.sum<-runInfo(cuff)
RNAseq.run.info<-rna.seq.sum[1,2]
RNAseq.run.info
replicates.info<-replicates(cuff)
replicates.info
groups<-factor(replicates.info$sample_name)
samples<-replicates.info$rep_name
groups
samples
over="LUTS"
under="CTRL"

#######################################
## initialize
##########################################
ref.gtf<-"/home/drew/umb_triley/Reference-Genomes/Human/UCSC_hg38/genes.gtf"
hg38.gff<-readGFF(ref.gtf)
hg38.gr<-readGFFAsGRanges(ref.gtf)
library(biomaRt) #functions to create a genetable from a gff3
hg38.gtf<-clusterProfiler::Gff2GeneTable("/home/drew/umb_triley/Reference-Genomes/Human/UCSC_hg38/genes.gtf")
load("geneTable.rda")
head(hg38.gff)
edb<-geneTable$GeneID
head(geneTable$GeneName)

library(GenomicRanges)
z <- GRanges("chr1",IRanges(1000001,1001000),strand="+")
start(z);end(z);width(z);strand(z);mcols(z) # the 'metadata columns', any information stored alongside each range
ranges(z);seqnames(z); seqlevels(z); seqlengths(z)
library(BSgenome.Hsapiens.UCSC.hg38)
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

#Rsamtools scanBam returns lists of raw values from BAM files

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# extracting information from txdb
g <- genes(txdb) # GRanges, just start to end, no exon/intron information
tx <- transcripts(txdb) # GRanges, similar to genes()
e <- exons(txdb) # GRanges for each exon
ebg <- exonsBy(txdb, by="gene") # exons grouped in a GRangesList by gene
ebt <- exonsBy(txdb, by="tx") # similar but by transcript
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
ebg <- exonsBy(txdb, by="gene")
# see yieldSize argument for restricting memory
se <- summarizeOverlaps(ebg, bf)
saveDb(txdb, file="txdb.sqlite")
loadDb("txdb.sqlite")
# in Bioconductor >= 3.1, also makeTxDbFromGRanges
setwd("~/media/drew/Seagate4TbExtHDD/umb_triley/urine1/cuffcompare_results_hg38_gtf_guided")
annot.files <- list.files(pattern = "./Sample_[3:20]_out/transcripts.gtf")
setwd(dir = "/home/drew/umb_triley/urine1/cuffcompare_results_hg38_gtf_guided")
dir()
for (i in seq_along(annot.files)) {
    assign(paste("cuffcompare", annot.files[i], sep = "_"), readGFF(annot.files[i]))}

#setwd("../cuffcompare_results_hg38_gtf_guided")
annot.files <- list.files(pattern = ".combined.gtf",all.files = F)

for (i in seq_along(annot.files)) {
   assign(annot.files[i],readGFF(annot.files[i]))}

for (i in annot.files){
    assign(paste("cuffgranges",i, sep = "_"), readGFFAsGRanges(i))}



cuffcmpdir<-c(".", pattern = "^[Sample]", full.names = TRUE, ignore.case = F)
cuffgtfs<-dir(cuffcmpdir, pattern = "transcripts.gtf", full.names = TRUE, ignore.case = TRUE)
cuff.notmap<-grep(cuffgtfs,pattern = "map",invert = T)
annot.files<-cuffgtfs[cuff.notmap]
cuffgtf.1<-readGFF(annot.files[1])


cuffdirs<-dir(".", pattern = "^[Sample]", full.names = TRUE, ignore.case = F)
cuffgtfs<-dir(cuffdirs, pattern = "transcripts.gtf", full.names = TRUE, ignore.case = TRUE)
cuff.notmap<-grep(cuffgtfs,pattern = "map",invert = T)
annot.files<-cuffgtfs[cuff.notmap]
cuffgtf.1<-readGFF(annot.files[1])
for (i in seq_along(annot.files)) {
    x<-gsub(x = annot.files[i], pattern="_out/",replacement = "_")
        assign(paste(x, annot.files[i], sep = ""), readGFF(annot.files[i]))}

for (i in seq_along(annot.files)) {
    assign(paste("cuffcompare", annot.files[i], sep = "_"), readGFF(annot.files[i]))}
for (i in seq_along(annot.files)){
    assign(paste("cuffcompare.granges", annot.files[i], sep = "_"), readGFFAsGRanges(annot.files[i]))}

for (i in seq_along(annot.files)){
  print(table(annot.files[[i["class_code"]]]))}

    assign(paste("cuffcompare.granges", annot.files[i], sep = "_"), readGFFAsGRanges(annot.files[i]))}
novel.cuffcompare_cuffcmp.combined.gtf<-cuffcompare_cuffcmp.combined.gtf["class_code"!="="]
head(novel.cuffcompare_cuffcmp.combined.gtf)

# then get the transcript sequence
txSeq <- extractTranscriptSeqs(Hsapiens, ebt)


tophatdirs<-dir(".", pattern = "^[Sample]", full.names = TRUE, ignore.case = F)
topbams<-dir(tophatdirs, pattern = "accepted_hits.bam", full.names = TRUE, ignore.case = TRUE)
top.notbai<-grep(topbams,pattern = ".bai",invert = T)
bam.files<-cuffgtfs[top.notbai]
annot.files<-bam.files
bf <- BamFileList(bam.files)

library(GenomicAlignments)
fls <- list.files(pattern="*.bam$")
bflst <- BamFileList(bf)
library(BiocParallel)
register(MulticoreParam(4))
# lots of options in the man page
# singleEnd, ignore.strand, inter.features, fragments, etc.

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
table()
table(annot.files[1][,"class_code"],annot.files[1][,"exon_number"])
)
      library(Rsamtools)
which <- GRanges("chr1",IRanges(1000001,1001000))
what <- c("rname","strand","pos","qwidth","seq")
param <- ScanBamParam(which=which, what=what)
# for more BamFile functions/details see ?BamFile
# yieldSize for chunk-wise access
bamfile <- BamFile("/path/to/file.bam")
bamfile <- BamFile("./")
reads <- scanBam(bamfile, param=param)
res <- countBam(bamfile, param=param)
# for more sophisticated counting modes
# see summarizeOverlaps() below
dirs <- dir(path=".Sample_*_out")
fls <- list.files(pattern="*.bam$")

cuffgtf.1<-readGFF(annot.files[1])
for (i in seq_along(annot.files)) {
    x<-gsub(x = annot.files[i], pattern="_out/",replacement = "_")
        assign(paste(x, annot.files[i], sep = ""), readGFF(annot.files[i]))}

for (i in seq_along(annot.files)) {
    assign(paste("cuffcompare", annot.files[i], sep = "_"), readGFF(annot.files[i]))}
for (i in seq_along(annot.files)){
    assign(paste("cuffcompare.granges", annot.files[i], sep = "_"), readGFFAsGRanges(annot.files[i]))}

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
cuffgtf.1<-readGFF(hg38.gtf)

# or build a txdb from GTF file (e.g. downloadable from Ensembl FTP site)
txdb <- makeTranscriptDbFromGFF("file.GTF", format="gtf")

# or build a txdb from Biomart (however, not as easy to reproduce later)
txdb <- makeTranscriptDbFromBiomart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

##########################################
## Genomic features
##########################################
## ---------------------------------------------------------------------
## The upstream sequences located in
##   http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/
## are based on RefSeq genes (RefSeq Genes track in the Genome Browser).
## Upstream sequences based on UCSC genes (UCSC Genes track in the
## Genome Browser) can easily be extracted from the full genome
## sequences with:
setwd(dir = "/media/drew/2f1bee59-e3b3-4af3-9055-b88f7a7084e9/umb_triley/urine1/cuffcompare_results_hg38_gtf_guided")
annot.files <- list.files(pattern = ".combined.gtf")
for (i in seq_along(annot.files)) {
    assign(paste("cuffcompare", annot.files[i], sep = "_"), readGFF(annot.files[i]))}
for (i in seq_along(annot.files)) {
    assign(paste(paste("cuffcompare.granges", annot.files[i], sep = "_"), readGFFAsGRanges(annot.files[i]))}

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
hg38.genome <- BSgenome.Hsapiens.UCSC.hg38
seqlengths(hg38.genome)
hg38.knownGene_txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
hg38.knownGene_up1000seqs <- extractUpstreamSeqs(hg38.genome, hg38.knownGene_txdb)
hg38.knownGene_txdb
hg38.knownGene_up1000seqs
hg38.knownGene_up1000seqs

hg38.refGene_txdb <- makeTxDbFromUCSC("hg38", "refGene")
refGene_up1000seqs <- extractUpstreamSeqs(genome, refGene_txdb)
refGene_txdb
refGene_up1000seqs
library(RMariaDB)
######################################################
#makeTxDbFromUCSC(), makeTxDbFromBiomart(), makeTxDbFromGFF()
######################################################
library(BSgenome.Hsapiens.NCBI.GRCh38)
genome <- BSgenome.Hsapiens.NCBI.GRCh38
seqlengths(genome)
library(TxDb.Hsapiens.NCBI.GRCh38.knownGene)
##############################################################################################
# code chunk for debugging purposes
# Must have all cuffdiff output (from *-over-* directory)
# as well as ref genome fasta and annotation gtf file
# in the current working directory
#############################################################################################
geneListString<-c("CXCL12","TGFB","MYC","RAS","CXCR4","IL8","IL6")
# cuff<-readCufflinks(dir=".", gtfFile="cuffcmp.combined.gtf",genome="genome.fa",rebuild=F)
cuff<-readCufflinks(dir=".", gtfFile="GRCh38_genes.gtf",rebuild=F)
cuff
inDir="LUTS";outDir="CTRL"
over="LUTS";under="CTRL"
pdf(file.path("CummeRbund_GRCh38_gtf_only.pdf"))
gtfFile="GRCh38_genes.gtf"

hg38gtfdb <- makeTxDbFromGFF(file="GRCh38_genes.gtf",format="gtf",
                             organism="Homo sapiens")

gtfdir<-as.matrix(dir())
gtfdir
gtfdir<-as.data.frame(dir(), stringsAsFactors = T)
gtf.files<-grep(pattern="combined.gtf", gtfdir)\
gtffiles<-list.files(pattern="combined.gtf")
gtf.dir<-list()
sample.gtfs<-gtfdir[gtf.files]

sapply(filez,FUN=function(eachPath){
  file.rename(from=eachPath,to=sub(pattern=”xxx”,replacement=”newTextString”,eachPath))

})
for (i in sample.gtfs){
  i<-readGFF(i)
  names(gtf.dif)[i] <- sample.gtfs
  #Matching the names of the elements to the filenames  file.rename(from = i, to=eval(i))
}
x <- as.list(rnorm(10000))
names(x) <- paste("a", 1:length(x), sep = "")
list2env(x , envir = .GlobalEnv)


###########################################################################
setwd(dir="umb_triley/urine1/cuffdiff_results_hg38_gtf_guided")
setwd(dir="LUTS-over-CTRL")

GTFFile="../../cuffmerge_results_hg38_gtf_guided/merged.gtf"
mergedgtf <- readGFF(GTFFile,)
merged.gtf<-as.data.frame(mergedgtf,na.rm=F)
head(merged.gtf)
dim(merged.gtf)
merged.granges<-makeGRangesFromDataFrame(merged.gtf, keep.extra.columns=TRUE)
merged.granges
merged.gtf.db <- makeTxDbFromGFF(file = "../../cuffmerge_results_hg38_gtf_guided/merged.gtf",
                                   format="gtf", organism="Homo sapiens")
merged.gtf.db
# saving and loading

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


############################################################################
# geneFPKMexpression-matrix ####
############################################################################

gene_exp.diff<-diffData(genes(cuff))
head(gene_exp.diff)
dim(gene_exp.diff)

g.rep.matrix<-repFpkmMatrix(genes(cuff))
head(g.rep.matrix);dim(g.rep.matrix)

g.cnt.matrix<-repCountMatrix(genes(cuff))
head(g.cnt.matrix);dim(g.cnt.matrix)

iso.rep.matrix<-as.data.frame(repFpkmMatrix(isoforms(cuff)))
row.names(iso.rep.matrix)->tracking_id
head(iso.rep.matrix)

gene.rep.ma<-cbind(tracking_id, iso.rep.matrix)
head(gene.rep.ma)

gene.xloc.matrix<-featureNames(isoforms(cuff))
head(gene.xloc.matrix)

gene.list<-getGenes(cuff,geneId = gene.xloc.matrix)
gene.list

gene_annotation_data<-featureNames(gene.list@isoforms)
head(gene_annotation_data)

gene.features<-annotation(isoforms(cuff))
head(gene.features)

iso_exp.diff<-diffData(isoforms(cuff))
head(iso_exp.diff)
dim(iso_exp.diff)

############################################################################
# Sig-geneFPKMexpression-matrix ####
############################################################################

mySigGenes<-getSig(cuff,x=over,y=under,alpha=.05,level='genes')
head(mySigGenes)
length(mySigGenes)
sigGenes<-getGenes(cuff, mySigGenes)
length(sigGenes)

sig_genes_annot<-sigGenes@annotation
head(sig_genes_annot)

sig_genes_exp.diff<-diffData(sigGenes)
head(sig_genes_exp.diff);dim(sig_genes_exp.diff)
gene.annot<-cbind(gene_symbol=sig_genes_annot$gene_short_name, xloc.ids=sig_genes_annot$gene_id)
head(gene.annot)
sig.genes_exp.diff<-as.data.frame(cbind(logFC=sig_genes_exp.diff$log2_fold_change,q_value=sig_genes_exp.diff$q_value),
											 row.names=sig_genes_exp.diff$gene_id)

sig.gene.annot<-cbind(gene.annot, sig.genes_exp.diff)
dim(sig.gene.annot);head(sig.gene.annot)
sig.gene.annot<-sig.gene.annot[,-2]
sig.u.gene.annot<-unique(sig.gene.annot$gene_symbol)
length(sig.u.gene.annot)
head(sig.gene.annot[sig.u.gene.annot,]);dim(sig.gene.annot[sig.u.gene.annot,])
sig.gene.df<-sig.gene.annot[sig.u.gene.annot,]
head(sig.gene.df);dim(sig.gene.df)
sig.genes_exp.df<-as.data.frame(cbind(logFC=sig.gene.df$logFC,q_value=sig.gene.df$q_value),
											 row.names=sig.gene.df$gene_symbol)

head(sig.genes_exp.df);dim(sig.genes_exp.df)
head(sig.genes_exp.diff);dim(sig.genes_exp.diff)
gene.features<-gene.features[,c(1,2,4)]
head(gene.features)

g.exp.df<-as.data.frame(merge(x=gene.features, y=gene.rep.ma, by.x="isoform_id", by.y="tracking_id"))#,row.names = "gene_short_name")
head(g.exp.df)

G.exp.df<-as.data.frame(g.exp.df[,-c(1:2)],row.names = g.exp.df[,3])
g.rep.matrix<-G.exp.df[,-1]
head(g.rep.matrix)
dim(g.rep.matrix)

head(g.rep.matrix[unique(row.names(g.rep.matrix)),])
g.rep.ma<-g.rep.matrix[unique(row.names(g.rep.matrix)),]
head(g.rep.ma);dim(g.rep.ma)
g.rep.ma<-as.matrix(g.rep.ma)

over.group<-grep(pattern = over, x = colnames(g.rep.ma),ignore.case = T)
under.group<-grep(pattern = under, x = colnames(g.rep.ma),ignore.case = T)

g.under.matrix<-g.rep.ma[,under.group]
head(g.under.matrix)
dim(g.under.matrix)

g.over.matrix<-g.rep.ma[,over.group]
head(g.over.matrix)
dim(g.over.matrix)

sig_over_gene_data<-subset(sig_genes_exp.diff, (significant =='yes') & (log2_fold_change > 0))
sig_over_gene_data<-subset(sig.genes_exp.df,logFC > 0)
head(sig_over_gene_data);dim(sig_over_gene_data)

sig_under_gene_data<-subset(sig_genes_exp.diff, (significant =='yes') & (log2_fold_change < 0))
sig_under_gene_data<-subset(sig.genes_exp.df,logFC < 0)
head(sig_under_gene_data);dim(sig_under_gene_data)

head(g.rep.ma);dim(g.rep.ma)
dim(sig_genes_exp.diff)

sig.genes<-row.names(g.rep.ma) %in% row.names(sig.genes_exp.df)
head(g.rep.ma[sig.genes,]);dim(g.rep.ma[sig.genes,])
s.g.rep.matrix<-g.rep.ma[sig.genes,]
head(g.rep.ma);dim(g.rep.ma)

g.under.genes<-row.names(g.rep.ma) %in% sig_under_gene_data$gene_id
g.under.genes<-row.names(g.rep.ma) %in% row.names(sig_under_gene_data)
g.u.rep.ma<-g.rep.ma[g.under.genes,]
head(g.u.rep.ma);dim(g.u.rep.ma)

g.over.genes<-row.names(g.rep.ma) %in% sig_over_gene_data$gene_id
g.over.genes<-row.names(g.rep.ma) %in% row.names(sig_over_gene_data)
g.o.rep.ma<-g.rep.ma[g.over.genes,]
head(g.o.rep.ma);dim(g.o.rep.ma)

s.g.over.matrix<-g.o.rep.ma[,over.group]
head(s.g.over.matrix);dim(s.g.over.matrix)

s.g.under.matrix<-g.u.rep.ma[,under.group]
head(s.g.under.matrix);dim(s.g.under.matrix)


##########################################
## Install Packages
##########################################
#BSgenome.Hsapiens.UCSC.hg38","EnrichmentBrowser","GSA","EGSEA"))
#source("https://bioconductor.org/biocLite.R")
#biocLite(""FDb.InfiniumMethylation.hg19"")
#biocLite("GenomeInfoDbData")
#devtools::install_github("stephenturner/msigdf")
# Get the data and build the vignette (requires plyr, dplyr, tidyr, knitr, rmarkdown)
#devtools::install_github("stephenturner/msigdf", build_vignettes = TRUE, force=T)
#biocLite(c("ALL","AllSorts"))
#install_github("Oshlack/AllSorts")
#biocLite("SGSeq")
### R code from vignette source 'GOstatsForUnsupportedOrganisms.Rnw'

###################################################
### code chunk number 1: available Schemas
###################################################
library("AnnotationForge");library("org.Hs.eg.db")
library("GSEABase");library("GOstats")
available.dbschemas()
###################################################
### code chunk number 5: <make parameter
###################################################

org.hs.GO.db = toTable(org.Hs.egGO)
org.hs.go.frameData = data.frame(org.hs.GO.db$go_id, org.hs.GO.db$Evidence, org.hs.GO.db$gene_id)
head(org.hs.go.frameData)
org.hs.goFrame=GOFrame(org.hs.go.frameData,organism="Homo sapiens")
org.hs.goAllFrame=GOAllFrame(org.hs.goFrame)
#library("GSEABase")
org.hsa.gsc <- GeneSetCollection(org.hs.goAllFrame, setType = GOCollection())
head(org.hsa.gsc)

save(org.hsa.gsc,file = "../../Reference-Genomes/Org.hsa.eg.GO.db.RData")

universe = Lkeys(org.Hs.egGO)
genes=sig.gene.df$gene_symbol
#genes = universe[1:500]
params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",
                             geneSetCollection=org.hsa.gsc,
                             geneIds = genes,
                             universeGeneIds = universe,
                             ontology = "MF",
                             pvalueCutoff = 0.05,
                             conditional = FALSE,
                             testDirection = "over")


###################################################
### code chunk number 6: call HyperGTest
###################################################
Over <- hyperGTest(params)
head(summary(Over))
frame = toTable(org.Hs.egPATH)
keggframeData = data.frame(frame$path_id, frame$gene_id)
head(keggframeData)
keggFrame=KEGGFrame(keggframeData,organism="Homo sapiens")


###################################################
### code chunk number 8: KEGG Parameters
###################################################
gsc <- GeneSetCollection(keggFrame, setType = KEGGCollection())
universe = Lkeys(org.Hs.egGO)
genes = universe[1:500]
kparams <- GSEAKEGGHyperGParams(name="My Custom GSEA based annot Params",
                               geneSetCollection=gsc,
                               geneIds = genes,
                               universeGeneIds = universe,
                               pvalueCutoff = 0.05,
                               testDirection = "over")
kOver <- hyperGTest(params)
head(summary(kOver))


###################################################
### code chunk number 9: info
###################################################
toLatex(sessionInfo())


###################################################
library(cummeRbund)
#source('openGraphSaveGraph.R');
options(error=traceback)
# Gets arguments that were passed in via command line
args = commandArgs(TRUE)
for (i in 1:length(args)) {
  eval(parse(text=args[[i]]))
}

###################################################
### code chunk number 1: args
###################################################
source("https://bioconductor.org/biocLite.R")
# These arguments are passed in on the command line via launchR.eachDiffDir.sh
# args:
# diffDir=\"${diffDir}\" inDir=\"${INPUT}/${diffDir}\" outDir=\"${OUTPUT}/${diffDir}\"  under=\"${under}\" over=\"${over}\" Rplots=\"$\{OUTPUT}/${diffDir}/${R_PLOTS}\" FPKMmatrix=\"${OUTPUT}/${diffDir}/${FPKM_MATRIX}\" DiffTable=\"${OUTPUT}/${diffDir}/${DIFF_TABLE}\"
##########################################
## genomic and genic features
##########################################
# Genome Data
library(GenomeInfoDbData)
library(TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts)
library(GenomicRanges);library(rtracklayer)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
# Gene Models
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# Annotation Maps
library(org.Hs.eg.db);library(biomaRt)
library(GenomicAlignments);library(GenomicFeatures)
library(GenomicRanges);library(ReactomePA);library(reactome.db)
library(KEGG.db);library(GO.db);library(AnnotationDbi)
library(IRanges);library(annotate)
biocLite(c("TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts","SNPlocs.Hsapiens.dbSNP144.GRCh38", "BSgenome.Hsapiens.NCBI.GRCh38"))
library(transcriptR)
# DBsnp 137
library(RITANdata);library(RITAN)
library(RCy3)
library(R)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
# Methylation Data
library(FDb.InfiniumMethylation.hg19)
library(FDb.UCSC.snp137common.hg19)

length(colnames(merged.gtf))
mergedGTFa<-merged.gtf[,1:11]
mergedGTFa<-merged.gtf[,1:8]
head(mergedGTFa)
mergedGTFb<-merged.gtf[,10:18]
head(mergedGTFb)

library(org.Hs.eg.db)
class(org.Hs.eg.db)
head(keys(org.Hs.eg.db))
class(keys(org.Hs.eg.db))

OrgHsEG <- org.Hs.egSYMBOL2EG
OrgHsKEGG <- org.Hs.egPATH2EG
OrgHseHG <- org.Hs.egSYMBOL
OrgHsNR <- org.Hs.egREFSEQ
OrgHsGO <- org.Hs.egGO2EG

#store the first six keys
my_keys <- head(keys(org.Hs.eg.db))
keytypes(org.Hs.eg.db)
#same as above
columns(org.Hs.eg.db)
#selecting
merged.org.df <- select(org.Hs.eg.db, keys = mergedGTFb$gene_name, columns=c("SYMBOL","ENTREZID", "REFSEQ"), keytype = "REFSEQ")
colnames(merged.org.df)<-c("gene_name", "symbol", "entrez_id")

head(merged.org.df)
head(mergedGTFb)
head(mergedGTFa)

dim(merged.org.df)
dim(mergedGTFb)
dim(mergedGTFa)
fixed.gtf<-cbind(mergedGTFa, merged.org.df, mergedGTFb)
colnames(fixed.gtf)<-c("seqid","source","type","start", "end", "score", "strand","phase","gene_id",
                       "symbol","entrez_id", "transcript_id", "exon_number", "gene_name", "oId", "nearest_ref",
                       "class_code", "tss_id", "p_id", "contained_in")
head(fixed.gtf)
dim(fixed.gtf)

fixedGTF<-file.path("../../cuffmerge_results_hg38_gtf_guided/fixed.merged.gtf")
write.table(fixed.gtf,file = fixedgtf,)






"../../cuffmerge_results_hg38_gtf_guided/fixed.merged.gtf"
select(org.Hs.eg.db,
       keys = my_keys,
       columns=c("ENTREZID","SYMBOL","GENENAME"),
       keytype="ENTREZID")

a <- select(org.Hs.eg.db,
            keys = my_genes,
            columns=c("ENTREZID", "SYMBOL","OMIM"),
            keytype="SYMBOL")
a
#symbol and transcript location
b <- select(TxDb.Hsapiens.UCSC.hg19.knownGene,
            keys = a$ENTREZID,
            columns=c('GENEID', 'TXCHROM', 'TXSTART', 'TXEND', 'TXID'),
            keytype="GENEID")
b
c <- select(org.Hs.eg.db, keys = mygenes, columns=c("REFSEQ", "SYMBOL"), keytype = "REFSEQ")
names(b) <- c('ENTREZID', 'TXID', 'TXCHROM', 'TXSTART', 'TXEND')
c <- merge(a, b, 'ENTREZID')
c
columns(org.Hs.eg.db)
#find out the columns we can use to search the database
keytypes(org.Hs.eg.db)
#which chromosome does the gene TP53 reside?
select(org.Hs.eg.db, keys="TP53", cols=c("SYMBOL", "CHR"), keytype="SYMBOL")
#find out the keys for the keytype CHR
NRids<-keys(org.Hs.eg.db, "REFSEQ")
#this question was asked in the lecture
#How many genes are there on chromosome 22 are in this annotation database?
#store all symbols
symbol <- keys(org.Hs.eg.db, "SYMBOL")
#how many gene symbols
length(symbol)
length(NRids)
#distribution of gene symbols along the chromosomes
symbol_chr <- select(org.Hs.eg.db, keys=symbol, cols=c("CHR","SYMBOL"), keytype="SYMBOL")
#the above warning is for duplicated rows
#double check how many gene symbols are in symbol_chr
length(symbol_chr$SYMBOL)
#unique ones
length(unique(symbol_chr$SYMBOL))
#distribution of symbols on chromosomes
table(symbol_chr$CHR)
plot(table(symbol_chr$CHR), xlab="Chromosomes", ylab="Number of genes")

##########################################
## load biomart
##########################################
library("biomaRt"); listMarts()
ensembl=useMart("ensembl");ensembl
listDatasets(ensembl)
listDatasets(ensembl)$version
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensembl
filters = listFilters(ensembl)
filters
attributes = listAttributes(ensembl)
attributes
genes
GeneSyms<-getBM(attributes=c("refseq_mrna", "hgnc_symbol"),
                filters = 'refseq_mrna', values = genes, mart = ensembl)
GeneSyms
dim(GeneSyms)
length(genes)


GeneNameAnnots<-getBM(attributes=c("refseq_mrna", "entrezgene","hgnc_symbol"),
                      filters = 'refseq_mrna', values = merged.gtf$gene_name, mart = ensembl)
GeneNameAnnots
dim(GeneNameAnnots)

PathAnnots<-getBM(attributes=c("refseq_mrna","go_id", "reactome","kegg_enzyme"),
                  filters = 'refseq_mrna', values = merged.gtf$gene_name, mart = ensembl)
PathAnnots
dim(PathAnnots)
length(PathAnnots)


library("biomaRt")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
#look for filters with RefSeq
grep("refseq", filters$name, ignore.case=T, value=T)
attributes = listAttributes(ensembl)
#RefSeq for beta actin
my_refseq <- 'NM_001101'
getBM(attributes='ensembl_gene_id', filters = 'refseq_mrna', values = my_refseq , mart = ensembl)
getBM(attributes=c('ensembl_gene_id','description'), filters = 'refseq_mrna', values = my_refseq , mart = ensembl)
library("biomaRt")
ensembl<- useMart("ensembl")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
#use the ensembl mart and the human dataset
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

#create a filter for all assembled human chromosomes
my_chr <- c(1:22, 'M', 'X', 'Y')

#listAttributes shows all attributes
attributes <- listAttributes(ensembl)

#find entrez attribute name
grep(pattern="entrez", x=attributes$description, ignore.case=T)
#[1] 45
attributes[45,]
#name   description
#45 entrezgene EntrezGene ID

#find refseq attribute name
grep(pattern="refseq", x=attributes$description, ignore.case=T)
#[1] 65 66 67 68 69 70

attributes[65:70,]
grep(pattern="ucsc", x=attributes$description, ignore.case=T)
attributes[73,]
#   name description
#73 ucsc     UCSC ID

#find Ensembl gene name
head(attributes[grep(pattern="ensembl gene",
                     x=attributes$description,
                     ignore.case=T
),
]
)
my_refseq_mrna <- getBM(attributes = 'refseq_mrna',
                        filters = 'chromosome_name',
                        values = my_chr,
                        mart = ensembl)

my_entrez_gene <- getBM(attributes = 'entrezgene',
                        filters = 'chromosome_name',
                        values = my_chr,
                        mart = ensembl
)

my_ucsc <- getBM(attributes = 'ucsc',
                 filters = 'chromosome_name',
                 values = my_chr,
                 mart = ensembl
)

my_ensembl_gene_id <- getBM(attributes = 'ensembl_gene_id',
                            filters = 'chromosome_name',
                            values = my_chr,
                            mart = ensembl
)

#how many identifiers from each database
length(my_ensembl_gene_id[,1])
length(my_entrez_gene[,1])
length(my_refseq_mrna[,1])
length(my_ucsc[,1])
getBM(attributes=c('refseq_mrna', 'ucsc'),
      filters = 'refseq_mrna',
      values = 'NM_001195597',
      mart = ensembl)
my_annotation <- getBM(attributes = c('ucsc', 'ensembl_gene_id', 'refseq_mrna', 'entrezgene'),
                       filters = 'chromosome_name',
                       values = my_chr,
                       mart = ensembl)

head(my_annotation)
#substitute the blank lines with NAs
my_annotation <- sapply(my_annotation,gsub,pattern="^$",replacement=NA)
#label the rows as sequential numbers
row.names(my_annotation) <- 1:length(my_annotation[,1])
head(my_annotation)


###########################################################################

novelmerged<-merged.gtf[which(merged.gtf["class_code"] != "="),]
head(novelmerged)
dim(novelmerged)

oldGTFFile="../../cuffdiff_results_GRCh38_default/GRCh38_genes.gtf"
oldmergedgtf <- readGFF(oldGTFFile)
grch38_genes.gtf<-as.data.frame(oldmergedgtf)
head(grch38_genes.gtf)
dim(grch38_genes.gtf)

grch38.granges<-makeGRangesFromDataFrame(grch38_genes.gtf, keep.extra.columns=TRUE)
grch38.granges
grch38.gtf.db <- makeTxDbFromGFF(file = "../../cuffdiff_results_GRCh38_default/GRCh38_genes.gtf",
                                 format="gtf", organism="Homo sapiens")
grch38.gtf.db
###########################################################################

##########################################
## Gene-IDs to Gene Symbols
##########################################
library(annotate)
#find the GO attribute name
grep("go",attributes$name, ignore.case=T, value=T)
grep("ensembl", filters$name, ignore.case=T, value=T)
test <- 'ENSG00000075624'
getBM(attributes="go_id", filters="ensembl_gene_id", values = test, mart = ensembl)

Term()
library("GO.db")
goooo <- getBM(attributes="go_id", filters="ensembl_gene_id", values = test, mart = ensembl)
Term(goooo)
na.omit(as.data.frame(Term(goooo)))
for(go in goooo[1:5]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}
#two Ensembl ids
test <- c('ENSG00000206172','ENSG00000075624')
getBM(attributes=c('ensembl_gene_id','go_id'), filters="ensembl_gene_id", values = test, mart = ensembl)


session <- browserSession("UCSC")
##########################################
## Gene-IDs to Gene Symbols
##########################################
library("biomaRt")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)

test <- 'ENSG00000118473'
getBM(attributes=c('ensembl_gene_id', "hgnc_symbol"), filters = "ensembl_gene_id", values=test, mart=ensembl)
ensembl_gene_id hgnc_symbol
1 ENSG00000118473       SGIP1

test <- c('ENSG00000118473', 'ENSG00000162426')
getBM(attributes=c('ensembl_gene_id', "hgnc_symbol"), filters = "ensembl_gene_id", values=test, mart=ensembl)
ensembl_gene_id hgnc_symbol
1 ENSG00000118473       SGIP1
2 ENSG00000162426     SLC45A1

getBM(attributes=c('ensembl_gene_id', "hgnc_symbol", "description"), filters = "ensembl_gene_id", values=test, mart=ensembl)

##########################################
## working with snps
##########################################

library(biomaRt)
listMarts()
snp <- useMart("snp",dataset="hsapiens_snp")
start <- Sys.time()
out=getBM(attributes=c("refsnp_id","allele","chrom_start"),
          filters=c("chr_name","start","end"),
          values=list(8,148350, 158612),
          mart=snp)
end <- Sys.time()
print(end - start)
Time difference of 1.544547 secs

nrow(out)
[1] 52
library(GO.db)
GO.db # metadata
class(GO.db)
keytypes(GO.db)
my_go_keys <- head(keys(GO.db))
my_go_keys
select(GO.db,
       keys = my_go_keys,
       columns=c("GOID", "TERM", "ONTOLOGY"),
       keytype="GOID")


library(AnnotationHub)
ah = AnnotationHub()
query(ah, "HepG2")
query(ah, c("HepG2", "H4K5"))


library(KEGG.db)
library(KEGGREST)
brca2K = keggGet("hsa:675")
names(brca2K[[1]])
brpat = keggGet("path:hsa05212")
names(brpat[[1]])
brpat[[1]]$GENE[seq(1,132,2)] # entrez gene ids
library(png)
library(grid)
brpng = keggGet("hsa05212", "image")
grid.raster(brpng)



k5 = keys(GO.db)[1:5]
cgo = columns(GO.db)
select(GO.db, keys=k5, columns=cgo[1:3])

con = GO_dbconn()
dbListTables(con)
dbGetQuery(con,
           "select _id, go_id, term from go_term limit 5")
