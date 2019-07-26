# gtf_cleanup_script
R.home();path.expand("~")
library(cummeRbund);library(msigdf)
pdf(file.path("CummeRbund_GRCh38_gtf_only.pdf"))
#source('openGraphSaveGraph.R');
options(error=traceback)
# Gets arguments that were passed in via command line
args = commandArgs(TRUE)
for (i in 1:length(args)) {
  eval(parse(text=args[[i]]))}
##########################################
## Genomic features
## The upstream sequences located in
##   http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/
## are based on RefSeq genes (RefSeq Genes track in the Genome Browser).
## Upstream sequences based on UCSC genes (UCSC Genes track in the
## Genome Browser) can easily be extracted from the full genome
## sequences with:
browseUCSCtrack(genome="hg38", tablename="knownGene",
                url="http://genome.ucsc.edu/cgi-bin/")

###################################################

###################################################
### code chunk number 1: args
###################################################
# These arguments are passed in on the command line via launchR.eachDiffDir.sh
# args:
# diffDir=\"${diffDir}\" inDir=\"${INPUT}/${diffDir}\" outDir=\"${OUTPUT}/${diffDir}\"  under=\"${under}\" over=\"${over}\" Rplots=\"$\{OUTPUT}/${diffDir}/${R_PLOTS}\" FPKMmatrix=\"${OUTPUT}/${diffDir}/${FPKM_MATRIX}\" DiffTable=\"${OUTPUT}/${diffDir}/${DIFF_TABLE}\"
####################################################################################
## load libraries
####################################################################################
library(devtools);library(sleuth);library(msigdf)
library(tidyverse);library(plyr);library(ggbio)
library(edgeR);library(limma);library(DESeq2)
library(gage);library(gageData);library(pathview)
library(EGSEA);library(GSA);library(SGSeq)
library(DOSE);library(EnrichmentBrowser)
library(Rsubread);library(ReactomePA)
library(rtracklayer);library(Rsamtools)

##########################################
## genomic and genic features
##########################################
# Genome Data
library(biomaRt);library(gageData);library(msigdf)
library(GenomicRanges);library(GenomeInfoDbData)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
# Gene Models
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# Annotation Maps
library(org.Hs.eg.db);library(GenomicAlignments)
library(GenomicFeatures);library(GenomicRanges)
library(IRanges);library(annotate)
library(ReactomePA);library(reactome.db)
library(KEGG.db);library(GO.db);library(AnnotationDbi)
library(transcriptR)

##########################################
hg38.fa <- BSgenome.Hsapiens.UCSC.hg38
seqlengths(hg38.fa)
mySession = browserSession("UCSC")
genome(mySession) <- "hg38"
hg38.txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
reg.chrs <- c(as.character(seq(1:23)), "X", "Y")
seqlevels(reg.chrs)

hg38.up1000seqs <- extractUpstreamSeqs(hg38.fa, hg38.txdb)
hg38.up1000seqs
#DEFAULT_CIRC_SEQS
tx_by_gene <- transcriptsBy(hg38.txdb, by="gene")
tx_by_gene
seqlevels(tx_by_gene)
ucsc.seqls=seqlevels(tx_by_gene)

hg38.transcripts <- transcripts(hg38.txdb)
hg38.exons <- exons(hg38.txdb)
hg38.cds <- cds(hg38.txdb)
cds <- keepSeqlevels(hg38.cds,)

genome(hg38.cds)<- "hg38"
hg38txseqs <- extractTranscriptSeqs(Hsapiens, TxDb.Hsapiens.UCSC.hg38.knownGene,   use.names=TRUE)

####################################################################################
hg19.fa <- BSgenome.Hsapiens.UCSC.hg19
seqlengths(hg19.fa)
mySession = browserSession("UCSC")
genome(mySession) <- "hg19"
hg19.txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
hg19.up1000seqs <- extractUpstreamSeqs(hg19.fa, hg19.txdb)
hg19.up1000seqs
DEFAULT_CIRC_SEQS
hg19.up1000seqs <- extractUpstreamSeqs(hg19.fa, hg19.txdb)
hg19.up1000seqs
hg19txseqs <- extractTranscriptSeqs(Hsapiens, TxDb.Hsapiens.UCSC.hg19.knownGene, use.names=TRUE)

######################################################
#makeTxDbFromUCSC(), makeTxDbFromBiomart(), makeTxDbFromGFF()
######################################################
library(BiocInstaller)
# Or do it in R
download.file("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
  "hg19ToHg38.over.chain.gz")
library(R.utils)
gunzip("hg19ToHg38.over.chain.gz")

#We will use it to convert the HepG2 binding address to hg38.

# Or do it in R
download.file("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
  "hg19ToHg38.over.chain.gz")
library(R.utils)
gunzip("hg19ToHg38.over.chain.gz")

#We will use it to convert the HepG2 binding address to hg38.

library(ph5)
library(ERBS)
data(HepG2)
library(rtracklayer)
ch = import.chain("hg19ToHg38.over.chain")
nHepG2 = liftOver(HepG2, ch)data(HepG2)
library(rtracklayer)
ch = import.chain("hg19ToHg38.over.chain")
nHepG2 = liftOver(HepG2, ch)

#biocLite("BSgenome.Hsapiens.UCSC.hg19.masked")
library(BSgenome.Hsapiens.NCBI.GRCh38)
genome <- BSgenome.Hsapiens.NCBI.GRCh38
seqlengths(genome)
library(TxDb.Hsapiens.NCBI.GRCh38)
##############################################################################################
# code chunk for debugging purposes
# Must have all cuffdiff output (from *-over-* directory)
# as well as ref genome fasta and annotation gtf file
# in the current working directory

library(org.Hs.eg.db)
class(org.Hs.eg.db)
head(keys(org.Hs.eg.db))
class(keys(org.Hs.eg.db))

OrgHsEG <- org.Hs.egSYMBOL2EG
OrgHsKEGG <- org.Hs.egPATH2EG
OrgHseHG <- org.Hs.egSYMBOL
OrgHsNR <- org.Hs.egREFSEQ
OrgHsGO <- org.Hs.egGO2EG

#############################################################################################
hg19gtfdb <- makeTxDbFromGFF(file="/home/drew/umb_triley/urine1/cuffcompare_results_hg19_gtf_guided/cuffcmp.combined.gtf",format="auto", organism="Homo sapiens")

hg38gtfdb <- makeTxDbFromGFF(file="/home/drew/umb_triley/urine1/cuffcompare_results_hg38_gtf_guided/cuffcmp.combined.gtf",format="auto", organism="Homo sapiens")
hg38gtfdb
grch38gtfdb <- makeTxDbFromGFF(file="/home/drew/umb_triley/Reference-Genomes/Human/NCBI_GRCh38/genes.gtf",format="auto", organism="Homo sapiens")
grch38gtfdb
hg19.gtf.db <- makeTxDbFromGFF(file =  "/home/drew/umb_triley/Reference-Genomes/Human/UCSC_hg19/UCSC_hg19.gtf",format="gtf", organism="Homo sapiens")

hg19.r.gtf.db <- makeTxDbFromGFF(file =  "/home/drew/umb_triley/Reference-Genomes/Human/UCSC_hg19/hg19_ribosomal.gtf",format="gtf", organism="Homo sapiens")

hg19.genes.gtf.db <- makeTxDbFromGFF(file =  "/home/drew/umb_triley/Reference-Genomes/Human/UCSC_hg19/genes.gtf",format="gtf", organism="Homo sapiens")


geneListString<-c("CXCL12","TGFB","MYC","RAS","CXCR4","IL8","IL6")
grch38.gtfFile="/home/drew/umb_triley/urine1/cuffcompare_results_grch38_gtf_guided/cuffcmp.combined.gtf"
grch38.fa="/home/drew/umb_triley/Reference-Genomes/Human/NCBI_GRCh38/genome.fa"
grchg38gtftx <- makeTxDbFromGFF(file=grch38.gtfFile,format="auto", organism="Homo sapiens")
grchg38gtftx
cuff<-readCufflinks(dir="/home/drew/umb_triley/urine1/cuffdiff_results_grch38_gtf_guided/LUTS-over-CTRL/", gtfFile="cuffcmp.combined.gtf",genome="genome.fa",rebuild=F)
inDir="LUTS";outDir="CTRL"
over="LUTS";under="CTRL"

grch38cmp.gtf <- readGFF(grch38.gtfFile)
grch38cmp.gtf
grch38.merged.granges<-makeGRangesFromDataFrame(grch38cmp.gtf, keep.extra.columns=TRUE)
grch38.merged.granges
grch38.merged.gtf.db <- makeTxDbFromGFF(file = "/home/drew/umb_triley/urine1/cuffcompare_results_grch38_gtf_guided/cuffcmp.combined.gtf",format="auto", organism="Homo sapiens")
grch38.merged.gtf.db

novelmerged<-merged.gtf[which(merged.gtf["class_code"] != "="),]
head(novelmerged)
dim(novelmerged)

GTFFile="/home/drew/umb_triley/urine1/cuffcompare_results_grch38_gtf_guided/cuffcmp.combined.gtf"
mergedgtf <- readGFF(GTFFile)

grch38.granges<-makeGRangesFromDataFrame(mergedgtf, keep.extra.columns=TRUE)
grch38.granges
#store the first six keys
my_keys <- head(keys(org.Hs.eg.db))
keytypes(org.Hs.eg.db)
#same as above
columns(org.Hs.eg.db)
#selecting
merged.org.df <- select(org.Hs.eg.db, keys = mergedGTFb$gene_name, columns=c("SYMBOL","ENTREZID", "REFSEQ"), keytype = "REFSEQ")
colnames(merged.org.df)<-c("gene_name", "symbol", "entrez_id")

colnames(fixed.gtf)<-c("seqid","source","type","start", "end", "score", "strand","phase","gene_id",
                       "symbol","entrez_id", "transcript_id", "exon_number", "gene_name", "oId", "nearest_ref",
                       "class_code", "tss_id", "p_id", "contained_in")
head(fixed.gtf)
dim(fixed.gtf)

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
