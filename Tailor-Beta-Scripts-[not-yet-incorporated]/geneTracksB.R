library(cummeRbund);library(ReactomePA); library(pathview)
#library(org.Mm.eg.db)
library(limma);library(DOSE);library(GO.db);library(org.Hs.eg.db)
library(topGO);library(GSEABase);library(clusterProfiler)
library(biomaRt);library(KEGG.db);library(reactome.db)
###############################################################################
# Features, Counts, and all inclusive tables 
# (CDS-p_ids, TSS-TSS ids, isoforms- XLOC and TCONS ids, and Genes XLOC ids)
# Significantly differentially expressed features: method 1
############################################################################
## ---CummeRbund and GSEA ------------------------------------------------------------------------------------------------
#source('~/R_Py_scripts/dcEnrichmentGO.R', echo=TRUE)
cuff<-readCufflinks(gtfFile="cuffcmp.combined.gtf",genome="genome.fa",rebuild=F)
cuff
refgtf<-read.table("cuffcmp.combined.gtf", sep=c("\t"), header=F)
head(refgtf)
genome="genome.fa"
over="LUTS"
under="CTRL"
###############################################################################
# Features, Counts, and all inclusive tables 
# get gene symbols for CDS TSS Isoforms features
################################################################################
runInfo(cuff);
repdata<-replicates(cuff)
reps<-repdata$rep_name
mySigGenes<-getSig(cuff,x=over,y=under,alpha=.05,level='genes')
head(mySigGenes)
length(mySigGenes)
sigGenes<-getGenes(cuff, mySigGenes)
length(sigGenes)

sig_genes_exp.diff<-diffData(sigGenes)
sig_genes_annot<-sigGenes@annotation
head(sig_genes_exp.diff)
head(sig_genes_annot)
sig_genes_exp.diff["gene_id"]<-sig_genes_annot["gene_short_name"]
head(sig_genes_exp.diff)

genes_exp.diff<-diffData(isoforms(cuff))
dim(genes_exp.diff)
g.rep.matrix<-repFpkmMatrix(isoforms(cuff))
head(g.rep.matrix)
gene.xloc.matrix<-featureNames(isoforms(cuff))
gene.list<-getGenes(cuff,geneId = gene.xloc.matrix)
gene.list
gene_annotation_data<-featureNames(gene.list@isoforms)
dim(gene_annotation_data)
gene.rep.matrix<-cbind2(gene_annotation_data, g.rep.matrix)
gene_exp.diff<-cbind2(gene_annotation_data, genes_exp.diff)
genes_exp.diff<-as.data.frame(cbind(gene_annotation_data, 
								   logFC=gene_exp.diff$log2_fold_change,
								   pvalue=gene_exp.diff$p_value,
								   qvalue=gene_exp.diff$q_value))
library(BSgenome.Hsapiens.UCSC.hg19)
bsg <- BSgenome.Hsapiens.UCSC.hg19
## Get the genome version
unique(genome(bsg))
library(Gviz)
library(GenomicFeatures)
library(EnsDb.Hsapiens.v75)
library(Rsamtools)
edb <- EnsDb.Hsapiens.v75
## List all available columns in the database.
columns(edb)
seqnames(bsg)
## Get the FaFile with the genomic sequence matching the Ensembl version
listColumns(edb)
keytypes(edb)
columns(edb)
## Get the FaFile with the genomic sequence matching the Ensembl version
listColumns(edb)
keytypes(edb)

gids <- keys(edb, keytype = "GENEID")
gids <- keys(mySigGenes     , keytype = "GENENAME")
length(gids)
## using the AnnotationHub package.
Dna <- getGenomeFaFile(edb)
## Get start/end coordinates of all genes.
genes <- genes(Dna)
## Subset to all genes that are encoded on chromosomes for which
## we do have DNA sequence available.
genes <- edb[seqnames(genes) %in% seqnames(seqinfo(Dna))]
gnames1 <- keys(edb, keytype = "GENENAME", filter = SeqnameFilter("Y"))
gnames2 <- keys(edb, keytype = "GENENAME", filter = SeqnameFilter((sig_genes_annot$gene_id)))
head(gnames1)
## Use the /standard/ way to fetch data.
siggenes<-select(edb, keys = c(sig_genes_annot$gene_id), keytype = "GENENAME",
       columns = c("GENEID", "GENENAME", "TXID", "TXBIOTYPE"))
allgenes<-select(edb, keys = names(gene_annotation_data$gene_short_name),
            TxbiotypeFilter("protein_coding"),
       columns = c("GENEID", "GENENAME", "TXID", "TXBIOTYPE"))

## Note: ordering of the results might not match ordering of keys!
## Get the gene sequences, i.e. the sequence including the sequence of
## all of the gene's exons and introns.
geneSeqs <- getSeq(siggenes, edb)
geneSeqs <- getSeq(allgenes)

## Extract the sequences of all transcripts encoded on chromosome Y.
yTx <- extractTranscriptSeqs(Dna, edb ) #E filter = SeqnameFilter("Y"))
yTxSeqs <- extractTranscriptSeqs(Dna, yTx)
yTxSeqs

## Along these lines, we could use the method also to retrieve the coding sequence
## of all transcripts on the Y chromosome.
cdsY <- cdsBy(edb, filter = SeqnameFilter(allgenes$GENENAME))
excds<-extractTranscriptSeqs(cdsY, edb)
## Retrieving a Gviz compatible GRanges object with all genes
## encoded on chromosome Y.
gr <- getGeneRegionTrackForGviz(edb,filter = SeqnameFilter(siggenes$GENENAME),chromosome="chr1", start = 20400000, end = 21400000)
## Define a genome axis track
gat <- GenomeAxisTrack()

## We have to change the ucscChromosomeNames option to FALSE to enable Gviz usage
## with non-UCSC chromosome names.
options(ucscChromosomeNames = FALSE)

plotTracks(list(gat, GeneRegionTrack(gr)))
options(ucscChromosomeNames = TRUE)
gr <- getGeneRegionTrackForGviz(edb,filter = SeqnameFilter(siggenes$GENENAME))  # chromosome=10, start = 20400000, end = 21400000)
## Define a genome axis track
gat <- GenomeAxisTrack()
plotTracks(list(gat, GeneRegionTrack(gr)))
## We have to change the ucscChromosomeNames option to FALSE to enable Gviz usage
## with non-UCSC chromosome names.
options(ucscChromosomeNames = FALSE)

plotTracks(list(gat, GeneRegionTrack(gr)))
options(ucscChromosomeNames = TRUE)

seqlevelsStyle(edb) <- "UCSC"
## Retrieving the GRanges objects with seqnames corresponding to UCSC chromosome names.
gr <- getGeneRegionTrackForGviz(edb, chromosome = "chr10", start = 20400000, end = 21400000)
seqnames(gr)
## Define a genome axis track
gat <- GenomeAxisTrack()
plotTracks(list(gat, GeneRegionTrack(gr)))

gr <- getGeneRegionTrackForGviz(edb, chromosome = "chr10", start = 20400000, end = 21400000)
seqnames(gr)
## Define a genome axis track
gat <- GenomeAxisTrack()
plotTracks(list(gat, GeneRegionTrack(gr)))

gr <- getGeneRegionTrackForGviz(edb, chromosome =c("chr12"), start = 0, end = 21400000)
seqnames(gr)
## Define a genome axis track
gat <- GenomeAxisTrack()
plotTracks(list(gat, GeneRegionTrack(gr)))
annTrack <- AnnotationTrack(genome="hg19", 
                            id=paste("annTrack item", 1:4),
                            name="annotation track foo",
                            stacking="squish")

protCod <- getGeneRegionTrackForGviz(edb, chromosome = "chr1",
                     start = 0, end = 2140000000,
                     filter = GenebiotypeFilter("protein_coding"))
lincs <- getGeneRegionTrackForGviz(edb, chromosome = "chr1",
                   start = 0, end = 2140000000,
                   filter = GenebiotypeFilter("lincRNA"))

plotTracks(list(gat, GeneRegionTrack(protCod, name = "protein coding"),
        GeneRegionTrack(lincs, name = "lincRNAs")), transcriptAnnotation = "symbol")
seqlevelsStyle <- "Ensembl"

yCds <- cdsBy(edb, filter = SeqnameFilter("Y"))
yCds


## Define the filter
grf <- GRangesFilter(GRanges("11", ranges = IRanges(114000000, 114000050),
                 strand="+"), condition = "overlapping")

## Query genes:
gn <- genes(edb, filter = grf)
gn
## GRanges object with 1 range and 5 metadata columns:
##                   seqnames                 ranges strand |         gene_id
##                      <Rle>              <IRanges>  <Rle> |     <character>
##   ENSG00000109906       11 [113930315, 114121398]      + | ENSG00000109906
##                     gene_name    entrezid   gene_biotype seq_coord_system
##                   <character> <character>    <character>      <character>
##   ENSG00000109906      ZBTB16        7704 protein_coding       chromosome
##   -------
##   seqinfo: 1 sequence from GRCh37 genome
## Next we retrieve all transcripts for that gene so that we can plot them.
txs <- transcripts(edb, filter = GenenameFilter(gn$gene_name))
plot(3, 3, pch = NA, xlim = c(start(gn), end(gn)), ylim = c(0, length(txs)),
     yaxt = "n", ylab = "")
## Highlight the GRangesFilter region
rect(xleft = start(grf), xright = end(grf), ybottom = 0, ytop = length(txs),
     col = "red", border = "red")
for(i in 1:length(txs)) {
    current <- txs[i]
    rect(xleft = start(current), xright = end(current), ybottom = i-0.975,
     ytop = i-0.125, border = "grey")
    text(start(current), y = i-0.5, pos = 4, cex = 0.75, labels = current$tx_id)
}



library(AnnotationHub)
ah <- AnnotationHub()

## Query all available files for Ensembl release 77 for
## Mus musculus.
query(ah, c("Homo sapien", "release-77"))
# AH21979 | Homo_sapiens.GRCh38.cdna.all.fa       
#  AH21980 | Homo_sapiens.GRCh38.dna_rm.toplevel.fa
#  AH21981 | Homo_sapiens.GRCh38.dna_sm.toplevel.fa
#  AH21982 | Homo_sapiens.GRCh38.dna.toplevel.fa   
#  AH21983 | Homo_sapiens.GRCh38.ncrna.fa          
#  AH21984 | Homo_sapiens.GRCh38.pep.all.fa        
#  AH28812 | Homo_sapiens.GRCh38.77.gtf          
## Get the resource for the gtf file with the gene/transcript definitions.
Gtf <- ah["AH28822"]
## Create a EnsDb database file from this.
DbFile <- ensDbFromAH(Gtf)
## We can either generate a database package, or directly load the data
edb <- EnsDb(DbFile)


## Identify and get the FaFile object with the genomic DNA sequence matching
## the EnsDb annotation.
Dna <- getGenomeFaFile(edb)
library(Rsamtools)
## We next retrieve the sequence of all exons on chromosome Y.
exons <- exons(edb, filter = SeqnameFilter("Y"))
exonSeq <- getSeq(Dna, exons)

## Alternatively, look up and retrieve the toplevel DNA sequence manually.
Dna <- ah[["AH22042"]]
load(system.file("YGRanges.RData", package = "ensembldb"))
Y
DB <- ensDbFromGRanges(Y, path=tempdir(), version = 75,
               organism = "Homo_sapiens")
## Warning in ensDbFromGRanges(Y, path = tempdir(), version = 75, organism
## = "Homo_sapiens"): I'm missing column(s): 'entrezid'. The corresponding
## database column(s) will be empty!
edb <- EnsDb(DB)
library(ensembldb)

## the GTF file can be downloaded from
# ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/
gtffile <- "GRCh38_genes.gtf"
## generate the SQLite database file
DB <- ensDbFromGtf(gtf = gtffile)
## load the DB file directly
EDB <- EnsDb(DB)

## alternatively, build the annotation package
## and finally we can generate the package
makeEnsembldbPackage(ensdb = DB, version = "0.99.12",
             maintainer = "Johannes Rainer <johannes.rainer@eurac.edu>",
             author = "J Rainer")

