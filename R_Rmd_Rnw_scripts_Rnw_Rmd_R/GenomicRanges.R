source("https://bioconductor.org/biocLite.R")
biocLite("alternativeSplicingEvents.hg19")
biocLite("alternativeSplicingEvents.hg38")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
biocLite("BSgenome.Hsapiens.UCSC.hg38")
biocLite("EnsDb.Hsapiens.v86")
biocLite("FDb.UCSC.snp137common.hg19")
biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene")
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")


library(EnsDb.Hsapiens.v86)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(alternativeSplicingEvents.hg19)
library(AnnotationHub)
hub <- AnnotationHub()

## Load the alternative splicing events for Human (hg19)
events <- query(hub, "alternativeSplicingEvents.hg19")[[1]]



for (i in 1:length(args)) {
  eval(parse(text=args[[i]]))
}

library('Rsubread')

readfile1 = paste0(inFile, "_R1.fastq.gz")
readfile2 = paste0(inFile, "_R2.fastq.gz")

print(paste0("readfile1=",readfile1))
print(paste0("readfile2=",readfile2))
# ----selectOrg2, eval=FALSE----------------------------------------------
library(HSAUR)
library(hackerschool)
library(R500)
play500()
library(devtools)

install_github("alyssafrazee/hackerschool",subdir="R500",force=T)
install_github("hackerschool", "alyssafrazee", subdir="R500")
for each KEGG pathway{
  
  ks_F <- ks.test( x = PVALUES OF GENES IN A KEGG PATHWAY),
  > y = PVALUES OF GENES NOT IN A KEGG PATHWAY ),
alternative = "greater" );
install.packages("refGenome")
library(refGenome)

# change to directory where you downloaded GTF file
setwd("~/muse/refgenome")

# create ensemblGenome object for storing Ensembl genomic annotation data
ens <- ensemblGenome()

# read GTF file into ensemblGenome object
read.gtf(ens, "Arabidopsis_thaliana.TAIR10.36.gtf")
[read.gtf.refGenome] Reading file 'Arabidopsis_thaliana.TAIR10.36.gtf'.
[GTF]   888100 lines processed.
[read.gtf.refGenome] Extracting genes table.
[read.gtf.refGenome] Found 32,833 gene records.
[read.gtf.refGenome] Finished.
# counts all annotations on each seqname
tableSeqids(ens)

1      2      3      4      5     Mt     Pt 
229818 132448 165053 131565 195121    690    567

# returns all annotations on mitochondria
extractSeqids(ens, 'Mt')
Object of class 'ensemblGenome' with 690 rows and 27 columns.
id seqid    source     feature start   end score strand frame   2"   4"   5"   3"   8"  H3"
886555 886555    Mt araport11  transcript   273   734     .      -     . <NA> <NA> <NA> <NA> <NA> <NA>
  886556 886556    Mt araport11        exon   273   734     .      -     . <NA> <NA> <NA> <NA> <NA> <NA>
  886557 886557    Mt araport11         CDS   276   734     .      -     0 <NA> <NA> <NA> <NA> <NA> <NA>
  886558 886558    Mt araport11 start_codon   732   734     .      -     0 <NA> <NA> <NA> <NA> <NA> <NA>
  886559 886559    Mt araport11  stop_codon   273   275     .      -     0 <NA> <NA> <NA> <NA> <NA> <NA>
  886561 886561    Mt araport11  transcript  8848 11415     .      -     . <NA> <NA> <NA> <NA> <NA> <NA>
  protein_version  protein_id           exon_id gene_name   gene_biotype exon_number   1
886555            <NA>        <NA>              <NA>   ORF153A protein_coding        <NA> <NA>
  886556            <NA>        <NA> ATMG00010.1.exon1   ORF153A protein_coding           1 <NA>
  886557               1 ATMG00010.1              <NA>   ORF153A protein_coding           1 <NA>
  886558            <NA>        <NA>              <NA>   ORF153A protein_coding           1 <NA>
  886559            <NA>        <NA>              <NA>   ORF153A protein_coding           1 <NA>
  886561            <NA>        <NA>              <NA>     RRN26           rRNA        <NA> <NA>
  transcript_biotype transcript_source transcript_id gene_source   gene_id
886555     protein_coding         araport11   ATMG00010.1   araport11 ATMG00010
886556     protein_coding         araport11   ATMG00010.1   araport11 ATMG00010
886557     protein_coding         araport11   ATMG00010.1   araport11 ATMG00010
886558     protein_coding         araport11   ATMG00010.1   araport11 ATMG00010
886559     protein_coding         araport11   ATMG00010.1   araport11 ATMG00010
886561               rRNA         araport11   ATMG00020.1   araport11 ATMG00020

# summarise features in GTF file
tableFeatures(ens)

CDS            exon  five_prime_utr     start_codon      stop_codon three_prime_utr 
285977          313952           56384           48315           48313           48308 
transcript 
54013

# create table of genes
my_gene <- getGenePositions(ens)
dim(my_gene)
[1] 32833    26

# gene IDs are unique
length(my_gene$gene_id)
[1] 32833
length(unique(my_gene$gene_id))
[1] 32833

# use dplyr to create more summaries
# number of genes on each seqname
my_gene %>% group_by(seqid) %>% summarise(n())
# A tibble: 7 x 2
seqid `n()`
<chr> <int>
  1     1  8656
2     2  5174
3     3  6435
4     4  4921
5     5  7362
6    Mt   152
7    Pt   133

my_gene %>% filter(seqid == "Mt") %>% select(gene_id) %>% head()
gene_id
1 ATMG00010
2 ATMG00020
3 ATMG00030
4 ATMG00040
5 ATMG00050
6 ATMG00060

my_gene %>% filter(seqid == 2) %>% select(gene_id) %>% head()
gene_id
1 AT2G01008
# length of genes
my_gene_length <- gt$end - gt$start
my_density <- density(my_gene_length)
plot(my_density, main = 'Distribution of gene lengths')

# what's the first peak?
which.max(my_density$y)
[1] 15
my_density$x[which.max(my_density$y)]
[1] 285.3376

# find peaks in a sequence of numbers
find_peak <- function (x){
  which(diff(sign(diff(x))) < 0) + 1
}
find_peak(my_density$y)
[1]  15  43 213 236 252 262 279 300 312 331 340 366 374 436 487 503

my_density$x[find_peak(my_density$y)[1:2]]
[1]  285.3376 1832.8927

plot(my_density, main = 'Distribution of gene lengths')
abline(v = c(my_density$x[find_peak(my_density$y)[1]],
             my_density$x[find_peak(my_density$y)[2]]),
       col = 'red', lty = 2)

text(x = my_density$x[find_peak(my_density$y)[1]] + 1500,
     y = my_density$y[find_peak(my_density$y)[1]],
     labels = round(my_density$x[find_peak(my_density$y)[1]], 2))

text(x = my_density$x[find_peak(my_density$y)[2]] + 1500,
     y = my_density$y[find_peak(my_density$y)[2]],
     labels = round(my_density$x[find_peak(my_density$y)[2]], 2)
     
     table(my_gene$gene_biotype)
     
     atlncRNA              atRNA             lncRNA              miRNA nontranslating_CDS 
     1037                 78               2444                325                 27 
     otherRNA     protein_coding               rRNA             snoRNA              snRNA 
     221              27628                 15                287                 82 
     tRNA 
     689
     
     my_gene %>% group_by(seqid, gene_biotype) %>% summarise(count = n()) -> my_tally
     
     ggplot(my_tally, aes(x = seqid, y = log2(count))) +
       geom_bar(aes(fill = gene_biotype), stat = 'identity', position = 'dodge')     
     
     
     library(GenomicRanges)
     my_gr <- with(my_gene, GRanges(seqid, IRanges(start, end), strand, id = gene_id))
     
     library(Gviz)
     ref <- GRanges()
     ref_track <- GenomeAxisTrack()
     options(ucscChromosomeNames=FALSE)
     data_track <- AnnotationTrack(my_gr, name = "Genes", showFeatureId = TRUE)
     plotTracks(c(ref_track, data_track),
                from = 1, to = 10000)
     library('rtracklayer')
     my_file <- "Arabidopsis_thaliana.TAIR10.36.gtf"
     my_obj <- import(my_file)
     
     class(my_obj)
     [1] "GRanges"
     attr(,"package")
     [1] "GenomicRanges"
     
     my_obj  
     GRanges object with 888095 ranges and 15 metadata columns:
       seqnames           ranges strand |    source        type     score     phase     gene_id
     <Rle>        <IRanges>  <Rle> |  <factor>    <factor> <numeric> <integer> <character>
       [1]        1     [3631, 5899]      + | araport11        gene      <NA>      <NA>   AT1G01010
     [2]        1     [3631, 5899]      + | araport11  transcript      <NA>      <NA>   AT1G01010
     [3]        1     [3631, 3913]      + | araport11        exon      <NA>      <NA>   AT1G01010
     [4]        1     [3760, 3913]      + | araport11         CDS      <NA>         0   AT1G01010
     [5]        1     [3760, 3762]      + | araport11 start_codon      <NA>         0   AT1G01010
     ...      ...              ...    ... .       ...         ...       ...       ...         ...
     [888091]       Pt [152806, 153195]      + | araport11         CDS      <NA>         0   ATCG01310
     [888092]       Pt [152806, 152808]      + | araport11 start_codon      <NA>         0   ATCG01310
     [888093]       Pt [153878, 154312]      + | araport11        exon      <NA>      <NA>   ATCG01310
     [888094]       Pt [153878, 154309]      + | araport11         CDS      <NA>         0   ATCG01310
     [888095]       Pt [154310, 154312]      + | araport11  stop_codon      <NA>         0   ATCG01310
     gene_name gene_source   gene_biotype transcript_id transcript_source transcript_biotype
     <character> <character>    <character>   <character>       <character>        <character>
       [1]      NAC001   araport11 protein_coding          <NA>              <NA>               <NA>
       [2]      NAC001   araport11 protein_coding   AT1G01010.1         araport11     protein_coding
     [3]      NAC001   araport11 protein_coding   AT1G01010.1         araport11     protein_coding
     [4]      NAC001   araport11 protein_coding   AT1G01010.1         araport11     protein_coding
     [5]      NAC001   araport11 protein_coding   AT1G01010.1         araport11     protein_coding
     ...         ...         ...            ...           ...               ...                ...
     [888091]      rpl2-B   araport11 protein_coding   ATCG01310.1         araport11     protein_coding
     [888092]      rpl2-B   araport11 protein_coding   ATCG01310.1         araport11     protein_coding
     [888093]      rpl2-B   araport11 protein_coding   ATCG01310.1         araport11     protein_coding
     [888094]      rpl2-B   araport11 protein_coding   ATCG01310.1         araport11     protein_coding
     [888095]      rpl2-B   araport11 protein_coding   ATCG01310.1         araport11     protein_coding
     exon_number           exon_id  protein_id protein_version
     <character>       <character> <character>     <character>
       [1]        <NA>              <NA>        <NA>            <NA>
       [2]        <NA>              <NA>        <NA>            <NA>
       [3]           1 AT1G01010.1.exon1        <NA>            <NA>
       [4]           1              <NA> AT1G01010.1               1
     [5]           1              <NA>        <NA>            <NA>
       ...         ...               ...         ...             ...
     [888091]           1              <NA> ATCG01310.1               1
     [888092]           1              <NA>        <NA>            <NA>
       [888093]           2 ATCG01310.1.exon2        <NA>            <NA>
       [888094]           2              <NA> ATCG01310.1               1
     [888095]           2              <NA>        <NA>            <NA>
       -------
       seqinfo: 7 sequences from an unspecified genome; no seqlengths   
     class(ens)
     [1] "ensemblGenome"
     attr(,"package")
     [1] "refGenome"