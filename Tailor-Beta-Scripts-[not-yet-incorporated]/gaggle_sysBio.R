source("http://bioconductor.org/workflows.R")
workflowInstall("ExpressionNormalizationWorkflow")
library("GenomicFeatures")

We indicate that none of our sequences (chromosomes) are circular using a 0-length character vector.

gtffile <- file.path(indir,"Homo_sapiens.GRCh37.75_subset.gtf")
txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())
txdb

## TxDb object:
## # Db type: TxDb
## # Supporting package: GenomicFeatures
## # Data source: /var/lib/jenkins/R/x86_64-pc-linux-gnu-library/3.4/airway/extdata/Homo_sapiens.GRCh37.75_subset.gtf
## # Organism: NA
## # Taxonomy ID: NA
## # miRBase build ID: NA
## # Genome: NA
## # transcript_nrow: 65
## # exon_nrow: 279
## # cds_nrow: 158
## # Db created by: GenomicFeatures package from Bioconductor
## # Creation time: 2017-04-10 14:26:22 -0700 (Mon, 10 Apr 2017)
## # GenomicFeatures version at creation time: 1.27.14
## # RSQLite version at creation time: 1.1-2
## # DBSCHEMAVERSION: 1.1
source("http://gaggle.systemsbiology.net/R/genome.browser.support.R")
gaggleInit()
ds <- getDatasetDescription()

getSequenceNames(ds)

len <- ds$sequences$chr$length
starts <- seq(1,len,100)
track.fwd <- data.frame(sequence='chr',
                        strand='+',
                        start=starts,
                        end=starts+99,
                        value=sin(starts/1000.0))
track.rev <- data.frame(sequence='chr',
                        strand='-',
                        start=starts+49,
                        end=starts+148,
                        value=sin(starts/900.0))
track <- rbind(track.fwd, track.rev)

head(track)
# data.frame in the format required by addTrack().

attr <- list(color='0x804B0082',source='Finklestein, et al. 2009', top=0.20,
             height=0.15, viewer='Scaling', group='bogus data')

addTrack(ds, track, name='waves', attributes=attr)


source('source_https.r')
## You can now source R scripts from GitHub. The RAW URL is needed.
source_https('https://raw.github.com/bobthecat/codebox/master/GO_over.r')
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## Define the universe
library(mouse4302.db)
uniqueId <- unique(as.vector(unlist(as.list(mouse4302ENTREZID))))
entrezUniverse <- uniqueId[!is.na(uniqueId)]
length(entrezUniverse)
[1] 20877

## ORA with conditional hypergeometric test
mfhyper <- GO_over(entrezUniverse, glist, annot='mouse4302.db')

## Information on the Directed Acyclic Graph (DAG)
goDag(mfhyper)
A graphNEL graph with directed edges
Number of Nodes = 2723
Number of Edges = 5643

## How many gene were mapped in the end?
geneMappedCount(mfhyper)
[1] 257

## Write out the results
source('https://raw.github.com/bobthecat/codebox/master/write.GOhyper.r')

mrnaGO <- write.GOhyper(mfhyper, filename="BP_mRNA_significant.xls")
dim(mrnaGO <- mrnaGO[mrnaGO$adjPvalue <= 0.05,])
[1] 59  8

head(mrnaGO)
