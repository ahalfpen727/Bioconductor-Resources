## ----warning=FALSE, message=FALSE----------------------------------------
library(EnsDb.Hsapiens.v75)

## Making a "short cut"
edb <- EnsDb.Hsapiens.v75
## print some informations for this package
edb

## for what organism was the database generated?
organism(edb)
    ## the resulting tables will be stored by default to the current working
    ## directory; if the correct Ensembl api (version 75) is defined in the
    ## PERL5LIB environment variable, the ensemblapi parameter can also be omitted.
    fetchTablesFromEnsembl(75,
                           ensemblapi="/home/bioinfo/ensembl/75/API/ensembl/modules",
                           species="human")

    ## containing the annotations
    DBFile <- makeEnsemblSQLiteFromTables()

    ## and finally we can generate the package
    makeEnsembldbPackage(ensdb=DBFile, version="0.0.1",
                         maintainer="Johannes Rainer <johannes.rainer@eurac.edu>",
                         author="J Rainer")

    ## Build an annotation database form a GFF file from Ensembl.
    ## ftp://ftp.ensembl.org/pub/release-83/gff3/rattus_norvegicus
    gff <- "Rattus_norvegicus.Rnor_6.0.83.gff3.gz"
    DB <- ensDbFromGff(gff=gff)
    edb <- EnsDb(DB)
    edb

    ## Build an annotation file from a GTF file.
    ## the GTF file can be downloaded from
    ## ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/
    gtffile <- "Homo_sapiens.GRCh38.10.gtf.gz"
    ## generate the SQLite database file
    DB <- ensDbFromGtf(gtf=paste0(ensemblhost, gtffile))

    ## load the DB file directly
    EDB <- EnsDb(DB)

    ## Alternatively, we could fetch a GTF file directly from AnnotationHub
    ## and build the database from that:
    library(AnnotationHub)
    ah <- AnnotationHub()
    ## Query for all GTF files from Ensembl for Ensembl version 81
    query(ah, c("Ensembl", "release-71", "GTF"))
    ## We could get the one from e.g. Bos taurus:
    DB <- ensDbFromAH(ah["AH47941"])
    edb <- EnsDb(DB)
    edb
#  AH133 | Homo_sapiens.GRCh37.69.cdna.all.fa       
#  AH134 | Homo_sapiens.GRCh37.69.dna.toplevel.fa   
#  AH135 | Homo_sapiens.GRCh37.69.dna_rm.toplevel.fa
#  AH136 | Homo_sapiens.GRCh37.69.dna_sm.toplevel.fa
#  AH137 | Homo_sapiens.GRCh37.69.ncrna.fa          
#  AH138 | Homo_sapiens.GRCh37.69.pep.all.fa        
## 
GRCh38.fa<-ah[["AH134"]]
GRCH38.GTF<-ah[["AH50844"]]
GRCh38.gtf<-query(ah, c("Ensembl", "9606", "GTF"))
GRCh38.fa<-query(ah, c("Ensembl", "9606", "fa"))
head(GRCh38.fa)
head(GRCH38.GTF)
head(GRCh38.fa)
head(GRCh38.gtf)


## Generate a sqlite database for genes encoded on chromosome Y
chrY <- system.file("chrY", package="ensembldb")
DBFile <- makeEnsemblSQLiteFromTables(path=chrY ,dbname=tempfile())
## load this database:
edb <- EnsDb(DBFile)

edb

## Generate a sqlite database from a GRanges object specifying
## genes encoded on chromosome Y
load(system.file("YGRanges.RData", package="ensembldb"))

Y

DB <- ensDbFromGRanges(Y, path=tempdir(), version=10,
                       organism="Homo_sapiens")
edb <- EnsDb(DB)


## ------------------------------------------------------------------------
Txsig <- transcripts(edb, filter = list(GenenameFilter(sig_genes_annot$gene_id)))
Txsig
Txall <- transcripts(edb, filter = list(GenenameFilter(gene_annotation_data$gene_short_name)))
Txall

## as this is a GRanges object we can access e.g. the start coordinates with
head(start(Txall))
head(start(Txsig))

## or extract the biotype with
head(Txall$tx_biotype)
head(Txsig$tx_biotype)
head(Txsig@ranges)
head(Txall@ranges)
head(Txsig@strand)
head(Txall@strand)
head(Txsig@elementMetadata)
head(Txall@seqnames)
## list all database tables along with their columns
listTables(edb)
## list columns from a specific table
listColumns(edb, "tx")
listColumns(edb, "gene")
listColumns(edb, "chromosome")
## ------------------------------------------------------------------------
Tx <- transcripts(edb,
		  columns = c(listColumns(edb , "gene"), "gene_name"),
		  filter = TxbiotypeFilter("nonsense_mediated_decay"),
		  return.type = "DataFrame")
nrow(Tx)
Tx
## ------------------------------------------------------------------------

yCds <-cdsBy(edb,columns=c(listColumns(edb , "gene"), "gene_name"), filter = GenenameFilter(genes_exp.diff$gene_short_name))
yCds
# ------------------------------------------------------------------------)
## Define the filter
gn <- GRangesFilter(GRanges("11", ranges = IRanges(100000000, 200000000),
			     strand="+"), condition = "overlapping")
grf <- GRangesFilter(GRanges("11", ranges = IRanges(0, 114000050),
			     strand="+"), condition = "overlapping")
grf
gn
## Next we retrieve all transcripts for that gene so that we can plot them.
txs <- transcripts(edb, filter = TxbiotypeFilter("nonsense_mediated_decay"),
				   columns = c(listColumns(edb , "gene"), "gene_name"))
txs
## ----tx-for-zbtb16, message=FALSE, fig.align='center', fig.width=7.5, fig.height=5----
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

## ------------------------------------------------------------------------
transcripts(edb, filter = grf)

## ------------------------------------------------------------------------
## Get all gene biotypes from the database. The GenebiotypeFilter
## allows to filter on these values.
listGenebiotypes(edb)
## Get all transcript biotypes from the database.
listTxbiotypes(edb)

## ------------------------------------------------------------------------
## We're going to fetch all genes which names start with BCL. To this end
## we define a GenenameFilter with partial matching, i.e. condition "like"
## and a % for any character/string.
BCLs <- genes(edb,
	      columns = c("gene_name", "entrezid", "gene_biotype"),
	      filter = list(GenenameFilter("BCL%", condition = "like")),
	      return.type = "DataFrame")
nrow(BCLs)
BCLs

## ------------------------------------------------------------------------
## determine the average length of snRNA, snoRNA and rRNA genes encoded on
## chromosomes X and Y.
mean(lengthOf(edb, of = "tx",
	      filter = list(GenebiotypeFilter(c("snRNA", "snoRNA", "rRNA")),
			    SeqnameFilter(c("X", "Y")))))

## determine the average length of protein coding genes encoded on the same
## chromosomes.
mean(lengthOf(edb, of = "tx",
	      filter = list(GenebiotypeFilter("protein_coding"),
			    SeqnameFilter(c("X", "Y")))))

## ------------------------------------------------------------------------
## Extract all exons 1 and (if present) 2 for all genes encoded on the
## Y chromosome
exons(edb, columns = c("tx_id", "exon_idx"),
      filter = list(SeqnameFilter("Y"),
		    ExonrankFilter(3, condition = "<")))

## ------------------------------------------------------------------------
TxByGns <- transcriptsBy(edb, by = "gene",
			 filter = list(SeqnameFilter(c("X", "Y"))))
TxByGns
ensDbFromAH(ah, outfile, path, organism, genomeVersion, version)

ensDbFromGRanges(x, outfile, path, organism, genomeVersion,
                 version)

ensDbFromGff(gff, outfile, path, organism, genomeVersion,
             version)

ensDbFromGtf(gtf, outfile, path, organism, genomeVersion,
             version)

fetchTablesFromEnsembl(version, ensemblapi, user="anonymous",
                       host="ensembldb.ensembl.org", pass="",
                       port=5306, species="human")

makeEnsemblSQLiteFromTables(path=".", dbname)

makeEnsembldbPackage(ensdb, version, maintainer, author,
                     destDir=".", license="Artistic-2.0")


## ----eval=FALSE----------------------------------------------------------
#  ## will just get exons for all genes on chromosomes 1 to 22, X and Y.
#  ## Note: want to get rid of the "LRG" genes!!!
  EnsGenes <- exonsBy(edb, by = "gene",
  		    filter = list(SeqnameFilter(c(1:22, "X", "Y")),
  				  GeneidFilter("ENSG%", "like")))

## ----eval=FALSE----------------------------------------------------------
#  ## Transforming the GRangesList into a data.frame in SAF format
  EnsGenes.SAF <- toSAF(EnsGenes)

## ----eval=FALSE----------------------------------------------------------
#  ## Create a GRanges of non-overlapping exon parts.
  DJE <- disjointExons(edb,
  		     filter = list(SeqnameFilter(c(1:22, "X", "Y")),
  				   GeneidFilter("ENSG%", "like")))

## ----eval=FALSE----------------------------------------------------------
  library(EnsDb.Hsapiens.v75)
  library(Rsamtools)
edb <- EnsDb.Hsapiens.v75
#  
#  ## Get the FaFile with the genomic sequence matching the Ensembl version
#  ## using the AnnotationHub package.
  Dna <- getGenomeFaFile(edb)
#  
#  ## Get start/end coordinates of all genes.
  genes <- genes(edb)
#  ## Subset to all genes that are encoded on chromosomes for which
#  ## we do have DNA sequence available.
  genes <- genes[seqnames(genes) %in% seqnames(seqinfo(Dna))]
#  
#  ## Get the gene sequences, i.e. the sequence including the sequence of
#  ## all of the gene's exons and introns.
  geneSeqs <- getSeq(Dna, genes)

## ----eval=FALSE----------------------------------------------------------
#  ## get all exons of all transcripts encoded on chromosome Y
  yTx <- exonsBy(edb, filter = SeqnameFilter("Y"))
#  
#  ## Retrieve the sequences for these transcripts from the FaFile.
  library(GenomicFeatures)
  yTxSeqs <- extractTranscriptSeqs(Dna, yTx)
  yTxSeqs
#  
#  ## Extract the sequences of all transcripts encoded on chromosome Y.
  yTx <- extractTranscriptSeqs(Dna, edb, filter = SeqnameFilter("Y"))
#  
#  ## Along these lines, we could use the method also to retrieve the coding sequence
#  ## of all transcripts on the Y chromosome.
  cdsY <- cdsBy(edb, filter = SeqnameFilter("Y"))
  extractTranscriptSeqs(Dna, cdsY)

## ----message=FALSE-------------------------------------------------------
## Change the seqlevels style form Ensembl (default) to UCSC:
seqlevelsStyle(edb) <- "UCSC"

## Now we can use UCSC style seqnames in SeqnameFilters or GRangesFilter:
genesY <- genes(edb, filter = list(SeqnameFilter("Y")))
  				   GeneidFilter("ENSG%", "like")))
## The seqlevels of the returned GRanges are also in UCSC style
seqlevels(genesY)

## ------------------------------------------------------------------------
seqlevelsStyle(edb) <- "UCSC"

## Getting the default option:
getOption("ensembldb.seqnameNotFound")

## Listing all seqlevels in the database.
seqlevels(edb)[1:22]

## Setting the option to NA, thus, for each seqname for which no mapping is available,
## NA is returned.
options(ensembldb.seqnameNotFound=NA)
seqlevels(edb)[1:23]

## Resetting the option.
options(ensembldb.seqnameNotFound = "ORIGINAL")

## ----warning=FALSE, message=FALSE----------------------------------------
library(BSgenome.Hsapiens.UCSC.hg19)
bsg <- BSgenome.Hsapiens.UCSC.hg19

## Get the genome version
unique(genome(bsg))
unique(genome(edb))
## Although differently named, both represent genome build GRCh37.

## Extract the full transcript sequences.
yTxSeqs <- extractTranscriptSeqs(bsg, exonsBy(edb, "tx", filter = SeqnameFilter("chrY")))

yTxSeqs

## Extract just the CDS
Test <- cdsBy(edb, "tx", filter = SeqnameFilter("chrY"))
yTxCds <- extractTranscriptSeqs(bsg, cdsBy(edb, "tx", filter = SeqnameFilter("chrY")))
yTxCds
byTxCds <- extractTranscriptSeqs(edb, cdsBy(edb, "tx", filter = SeqnameFilter("chrY")))
byTxCds

## ------------------------------------------------------------------------
seqlevelsStyle(edb) <- "Ensembl"

library(Gviz)
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
## encoded on chromosome Y.

gr <- getGeneRegionTrackForGviz(edb, chromosome = "Y", start = 20400000, end = 21400000)
## Define a genome axis track
gat <- GenomeAxisTrack()

## We have to change the ucscChromosomeNames option to FALSE to enable Gviz usage
## with non-UCSC chromosome names.
options(ucscChromosomeNames = FALSE)

plotTracks(list(gat, GeneRegionTrack(gr)))

options(ucscChromosomeNames = TRUE)

## ----message=FALSE-------------------------------------------------------
seqlevelsStyle(edb) <- "UCSC"
## Retrieving the GRanges objects with seqnames corresponding to UCSC chromosome names.
gr <- getGeneRegionTrackForGviz(edb, chromosome = "chr1", start = 20400000, end = 21400000)
seqnames(gr)
## Define a genome axis track
gat <- GenomeAxisTrack()
plotTracks(list(gat, GeneRegionTrack(gr)))

## ----gviz-separate-tracks, message=FALSE, warning=FALSE, fig.align='center', fig.width=7.5, fig.height=2.25----
protCod <- getGeneRegionTrackForGviz(edb, chromosome = "chr1",
				     start = 20400000, end = 21400000,
				     filter = GenebiotypeFilter("protein_coding"))
lincs <- getGeneRegionTrackForGviz(edb, chromosome = "chr1",
				   start = 1000000, end = 21400000,
				   filter = GenebiotypeFilter("lincRNA"))

plotTracks(list(gat, GeneRegionTrack(protCod, name = "protein coding"),
		GeneRegionTrack(lincs, name = "lincRNAs")), transcriptAnnotation = "symbol")

## At last we change the seqlevels style again to Ensembl
seqlevelsStyle <- "Ensembl"
## ------------------------------------------------------------------------
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
columns(edb)
listColumns(edb)
keytypes(edb)

## Get all gene ids from the database.
gids <- keys(edb, keytype = "GENEID")
length(gids)

## Get all gene names for genes encoded on chromosome Y.
gnames <- keys(edb, keytype = "GENENAME") #filter = SeqnameFilter("Y"))
head(gnames)
length(gnames)

## ----warning=FALSE-------------------------------------------------------
## Use the /standard/ way to fetch data.
gnames1<-select(edb, keys = c(gnames), keytype = "GENENAME",
       columns = c("GENEID", "GENENAME", "TXID", "TXBIOTYPE"))

## Use the filtering system of ensembldb
gnames2<-select(edb, keys = list(GenenameFilter(c(gnames)),
			TxbiotypeFilter("protein_coding")),
       columns = c("GENEID", "GENENAME", "TXID", "TXBIOTYPE"))
head(gnames2)
head(gnames1)
## ------------------------------------------------------------------------
## Use the default method, which just returns the first value for multi mappings.
siggie<-mapIds(edb, keys = c(sigGenes@ids), column = "TXID", keytype = "GENENAME")

siggies<-## Alternatively, specify multiVals="list" to return all mappings.
mapIds(edb, keys = c(sigGenes@ids), column = "TXID", keytype = "GENENAME",
       multiVals = "list")

## And, just like before, we can use filters to map only to protein coding transcripts.
mapIds(edb, keys = list(GenenameFilter(c(sigGenes@ids)),
			TxbiotypeFilter("protein_coding")), column = "TXID",
       multiVals = "list")

## ----eval=FALSE----------------------------------------------------------
  library(ensembldb)
#get all human gene/transcript/exon annotations from Ensembl (75)
#  ## the resulting tables will be stored by default to the current working
fetchTablesFromEnsembl(75, species = "human")
#  ## containing the annotations (again, the function assumes the required
#  ## txt files to be present in the current working directory)
  DBFile <- makeEnsemblSQLiteFromTables()
#  ## and finally we can generate the package
  makeEnsembldbPackage(ensdb = DBFile, version = "0.99.12",
  		     maintainer = "Johannes Rainer <johannes.rainer@eurac.edu>",
  		     author = "J Rainer")
#  ## Load the AnnotationHub data.
  library(AnnotationHub)
  ah <- AnnotationHub()
	#  ## Query all available files for Ensembl release GRCh38.v77 for
 query(ah, c("Homo sapien", "release-77"))
#  ## Get the resource for the gtf file with the gene/transcript definitions.
  Gtf <- ah["AH28812"]
hg19Bed <- ah[9606]  
ucsc <- ah["AH22650"]
ah[GM12892]
head(Gtf);head(hg19Bed)
  ## Create a EnsDb database file from this.
  DbFile <- ensDbFromAH(Gtf)
#  ## We can either generate a database package, or directly load the data
  edb <- EnsDb(DbFile)
# Identify and get the FaFile object with the genomic DNA sequence matching
#  ## the EnsDb annotation.
 Dna <- getGenomeFaFile(edb)
  library(Rsamtools)
#  ## We next retrieve the sequence of all exons on chromosome Y.
  exons <- exons(edb, filter = SeqnameFilter("Y"))
  exonSeq <- getSeq(Dna, exons)
#  
#  ## Alternatively, look up and retrieve the toplevel DNA sequence manually.
  Dna <- ah[["AH22042"]]

## ----message=FALSE-------------------------------------------------------
## Generate a sqlite database from a GRanges object specifying
# genes encoded on chromosome Y
load(system.file("YGRanges.RData", package = "ensembldb"))
Y

DB <- ensDbFromGRanges(Y, path=tempdir(), version = 75,
		       organism = "Homo_sapiens")

edb <- EnsDb(DB)
edb


## ----eval=FALSE----------------------------------------------------------
  library(ensembldb)
#  
#  ## the GTF file can be downloaded from
#  ## ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/
  gtffile <- "Homo_sapiens.GRCh37.75.gtf.gz"
#  ## generate the SQLite database file
  DB <- ensDbFromGtf(gtf = gtffile)
#  
#  ## load the DB file directly
  EDB <- EnsDb(DB)
#  
#  ## alternatively, build the annotation package
#  ## and finally we can generate the package
  makeEnsembldbPackage(ensdb = DB, version = "0.99.12",
  		     maintainer = "Johannes Rainer <johannes.rainer@eurac.edu>",
  		     author = "J Rainer")
###################################################
### code chunk number 11: Contextual results
###################################################
# Show the first four rows and first two columns 
# of the contextual association from the
# demonstration
results$conn_p_value[1:4, 1:2]

# Show the first four rows and first two columns 
# of the pathway overlap scores from the
# demonstration
results$pathway_overlap[1:4, 1:2]

library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

######   genes
#### get all transcripts of a gene
Tx <- transcripts(edb, filter = GenenameFilter(gene.list@ids), 
				  order.by="tx_biotype")
				  
Tx
Txw <- transcripts(edb, filter = GenenameFilter(genes_exp.diff$gene_short_name),
				   order.by=abs(as.numeric("seqnames")))
 Txw
## get all genes endcoded on chromosome Y

					  ## include all transcripts of the gene and their chromosomal
## coordinates, sort by chrom start of transcripts and return as
## GRanges.
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
columns(bsg)
## Get the FaFile with the genomic sequence matching the Ensembl version
listColumns(edb)
keytypes(edb)

gids <- keys(edb, keytype = "GENEID")
gids <- keys(mySigGenes     , keytype = "GENENAME")
length(gids)
## using the AnnotationHub package.
Dna <- getGenomeFaFile(edb)
## Subset to all genes that are encoded on chromosomes for which
## we do have DNA sequence available.
genes <- genes[seqnames(genes) %in% seqnames(seqinfo(Dna))]
gnames <- keys(edb, keytype = "GENENAME", filter = SeqnameFilter("Y"))
gnames <- keys(edb, keytype = "GENENAME", filter = "SEQLENGTH")   #SeqnameFilter("Y"))
head(gnames)
## Use the /standard/ way to fetch data.
siggenes<-select(edb, keys = c(sig_genes_annot$gene_id), keytype = "GENENAME",
       columns = c("GENEID", "GENENAME", "TXID", "TXBIOTYPE"))
allgenes<-select(edb, keys = names(gene_annotation_data$gene_short_name),
            TxbiotypeFilter("protein_coding"),
       columns = c("GENEID", "GENENAME", "TXID", "TXBIOTYPE"))
gnames <- keys(edb, keytype = "GENENAME", filter = SeqnameFilter("Y"))
gnames <- keys(edb, keytype = "GENENAME", filter = SeqnameFilter((sig_genes_annot$gene_id)))
head(gnames)
## Use the /standard/ way to fetch data.
siggenes<-select(edb, keys = c(sig_genes_annot$gene_id), keytype = "GENENAME",
       columns = c("GENEID", "GENENAME", "TXID", "TXBIOTYPE"))
allgenes<-select(edb, keys = names(gene_annotation_data$gene_short_name),
            TxbiotypeFilter("protein_coding"),
       columns = c("GENEID", "GENENAME", "TXID", "TXBIOTYPE"))

## Note: ordering of the results might not match ordering of keys!
## Get the gene sequences, i.e. the sequence including the sequence of
## all of the gene's exons and introns.
geneSeqs <- getSeq(Dna, genes)
geneSeqs <- getSeq(allgenes)

## Extract the sequences of all transcripts encoded on chromosome Y.
yTx <- extractTranscriptSeqs(Dna, edb ) #E filter = SeqnameFilter("Y"))
yTxSeqs <- extractTranscriptSeqs(Dna, yTx)
yTxSeqs

## Along these lines, we could use the method also to retrieve the coding sequence
## of all transcripts on the Y chromosome.
cdsY <- cdsBy(edb, filter = SeqnameFilter(allgenes$GENENAME))
extractTranscriptSeqs(cdsY, edb)

## Retrieving a Gviz compatible GRanges object with all genes
## encoded on chromosome Y.
gr <- getGeneRegionTrackForGviz(edb, chromosome = "Y", start = 20400000, end = 21400000)
## Define a genome axis track
gat <- GenomeAxisTrack()

## We have to change the ucscChromosomeNames option to FALSE to enable Gviz usage
## with non-UCSC chromosome names.
options(ucscChromosomeNames = FALSE)

plotTracks(list(gat, GeneRegionTrack(gr)))
options(ucscChromosomeNames = TRUE)

## Note: ordering of the results might not match ordering of keys!
## Get the gene sequences, i.e. the sequence including the sequence of
## all of the gene's exons and introns.
geneSeqs <- getSeq(Dna, genes)
geneSeqs <- getSeq(allgenes)

## Extract the sequences of all transcripts encoded on chromosome Y.
yTx <- extractTranscriptSeqs(Dna, edb ) #E filter = SeqnameFilter("Y"))
yTxSeqs <- extractTranscriptSeqs(Dna, yTx)
yTxSeqs

## Along these lines, we could use the method also to retrieve the coding sequence
## of all transcripts on the Y chromosome.
cdsY <- cdsBy(edb, filter = SeqnameFilter("Y"))
extractTranscriptSeqs(Dna, cdsY)

## Retrieving a Gviz compatible GRanges object with all genes
## encoded on chromosome Y.
gr <- getGeneRegionTrackForGviz(edb, chromosome = "Y", start = 20400000, end = 21400000)
## Define a genome axis track
gat <- GenomeAxisTrack()

## We have to change the ucscChromosomeNames option to FALSE to enable Gviz usage
## with non-UCSC chromosome names.
options(ucscChromosomeNames = FALSE)

plotTracks(list(gat, GeneRegionTrack(gr)))
options(ucscChromosomeNames = TRUE)

seqlevelsStyle(edb) <- "UCSC"
## Retrieving the GRanges objects with seqnames corresponding to UCSC chromosome names.
gr <- getGeneRegionTrackForGviz(edb, chromosome = "chrY", start = 20400000, end = 21400000)
seqnames(gr)
## Define a genome axis track
gat <- GenomeAxisTrack()
plotTracks(list(gat, GeneRegionTrack(gr)))

protCod <- getGeneRegionTrackForGviz(edb, chromosome = "chrY",
                     start = 20400000, end = 21400000,
                     filter = GenebiotypeFilter("protein_coding"))
lincs <- getGeneRegionTrackForGviz(edb, chromosome = "chrY",
                   start = 20400000, end = 21400000,
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


######   transcripts
##
## get all transcripts of a gene
Txall <- transcripts(edb,
                  filter=GeneidFilter(genes_exp.diff$tracking_id),
                  order(genes_exp.diff["logFC"], DECREASING=T))
Txall
)
## get all transcripts of two genes along with some information on the
## gene and transcript
Tx <- transcripts(edb,
                  filter=GeneidFilter(c("ENSG00000184895",
                      "ENSG00000092377")),
                      columns=c("gene_id", "gene_seq_start",
                          "gene_seq_end", "gene_biotype", "tx_biotype"))
Tx

######   promoters
##
## get the bona-fide promoters (2k up- to 200nt downstream of TSS)
promoters(edb, filter=GeneidFilter(c("ENSG00000184895",
                                     "ENSG00000092377")))

######   exons
##
## get all exons of the provided genes
Exon <- exons(edb,
              filter=GeneidFilter(c("ENSG00000184895",
                  "ENSG00000092377")),
              order.by="exon_seq_start",
              columns=c( "gene_id", "gene_seq_start",
                  "gene_seq_end", "gene_biotype"))
Exon



#####    exonsBy
##
## get all exons for transcripts encoded on chromosomes X and Y.
ETx <- exonsBy(edb, by="tx",
               filter=SeqnameFilter(c("X", "Y")))
ETx
## get all exons for genes encoded on chromosome 1 to 22, X and Y and
## include additional annotation columns in the result
EGenes <- exonsBy(edb, by="gene",
                  filter=SeqnameFilter(c("X", "Y")),
                  columns=c("gene_biotype", "gene_name"))
EGenes

## Note that this might also contain "LRG" genes.
length(grep(names(EGenes), pattern="LRG"))

## to fetch just Ensemblgenes, use an GeneidFilter with value
## "ENS%" and condition "like"


#####    transcriptsBy
##
TGenes <- transcriptsBy(edb, by="gene",
                        filter=SeqnameFilter(c("X", "Y")))
TGenes

## convert this to a SAF formatted data.frame that can be used by the
## featureCounts function from the Rsubreader package.
head(toSAF(TGenes))


#####   transcriptsByOverlaps
##
ir <- IRanges(start=c(2654890, 2709520, 28111770),
              end=c(2654900, 2709550, 28111790))
gr <- GRanges(rep("Y", length(ir)), ir)

## Retrieve all transcripts overlapping any of the regions.
txs <- transcriptsByOverlaps(edb, gr)
txs

## Alternatively, use a GRangesFilter
grf <- GRangesFilter(gr, condition="overlapping")
txs <- transcripts(edb, filter=grf)
txs


####    cdsBy
## Get the coding region for all transcripts on chromosome Y.
## Specifying also additional annotation columns (in addition to the default
## exon_id and exon_rank).
cds <- cdsBy(edb, by="tx", filter=SeqnameFilter("Y"),
             columns=c("tx_biotype", "gene_name"))

####    the 5' untranslated regions:
fUTRs <- fiveUTRsByTranscript(edb, filter=SeqnameFilter("Y"))

####    the 3' untranslated regions with additional column gene_name.
tUTRs <- threeUTRsByTranscript(edb, filter=SeqnameFilter("Y"),
                               columns="gene_name")

