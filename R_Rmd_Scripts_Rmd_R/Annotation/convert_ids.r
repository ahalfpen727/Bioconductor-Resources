
##### use bioconductor annotation packages #######

source("http://Bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
biocLite(c("GenomicFeatures", "AnnotationDbi"))

library("org.Hs.eg.db")
library("AnnotationDbi")
library("GenomicFeatures")

# all the possible mappings 
ls("package:org.Hs.eg.db")

# convert Entrez_ids to gene_symbols
myEntrez_ids <- c("1","10","100","1000","37690")
mySymbols<- mget(myEntrez_ids, org.Hs.egSYMBOL, ifnotfound=NA)
mySymbols
unlist(mySymbols)

# convert gene_symbols to Entrez_ids

mySymbols_2 <- c("VEGFA","CTCF", "SNAI1","KDM1A")
myEntrez_ids_2<- mget(mySymbols_2, org.Hs.egSYMBOL2EG, ifnotfound=NA)
unlist(myEntrez_ids_2)


?AnnotationDbi::mget  # get help

# or use the select function 
?AnnotationDbi::select
head(keys(org.Hs.eg.db))
keytypes(org.Hs.eg.db)
select(org.Hs.eg.db, keys = mySymbols_2, columns=c("SYMBOL","REFSEQ","GENENAME","ENTREZID"),keytype="SYMBOL")

select(org.Hs.eg.db, keys = myEntrez_ids, columns=c("SYMBOL","REFSEQ","GENENAME","ENTREZID"),keytype="ENTREZID")

# How many gene symbols 
symbol <- keys(org.Hs.eg.db, "SYMBOL")
length(symbol)

############### use biomart ###################

library(biomaRt)
mart<- useMart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')

# get sequences

seq <- getSequence(id = 'BRCA1', type='hgnc_symbol',seqType="3utr", mart = mart)  # pretty slow...
show(seq)

seq2 <-getSequence(id="ENST00000520540",type='ensembl_transcript_id',seqType='gene_flank', upstream =30, mart=mart)
show(seq2)


# convert gene ids  gene symbol to refseq

geneList<- c("VEGFA","CTCF", "SNAI1","KDM1A")
results<- getBM(attributes = c("refseq_mrna","hgnc_symbol"), filters="hgnc_symbol", values=geneList, mart=mart)
results

?getBM