## Define expressed regions for study SRP045638, only for chromosome 21
library('recount')
regions <- expressed_regions("SRP045638", "chr21", cutoff = 5L,
    maxClusterGap = 3000L)

## Import the Gencode v25 hg38 gene annotation
library("rtracklayer")
gencode_v25_hg38 <- import(paste0(
    "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human",
    "/release_25/gencode.v25.annotation.gff3.gz"))
            
## Keep only the chr21 info
gencode_v25_hg38 <- keepSeqlevels(gencode_v25_hg38, "chr21",
    pruning.mode = "coarse")

## Get the chromosome information for hg38
library("GenomicFeatures")
chrInfo <- getChromInfoFromUCSC("hg38")
chrInfo$chrom <- as.character(chrInfo$chrom)
chrInfo <- chrInfo[chrInfo$chrom %in% seqlevels(regions), ]
chrInfo$isCircular <- FALSE

## Assign the chromosome information to the object we will use to
## create the txdb object
si <- with(chrInfo, Seqinfo(as.character(chrom), length, isCircular,
    genome = "hg38"))
seqinfo(gencode_v25_hg38) <- si

## Switch from Gencode gene IDs to Ensembl gene IDs
gencode_v25_hg38$gene_id <- gsub("\\..*", "", gencode_v25_hg38$gene_id)

## Create the TxDb object
gencode_v25_hg38_txdb <- makeTxDbFromGRanges(gencode_v25_hg38)

## Explore the TxDb object
gencode_v25_hg38_txdb


library("bumphunter")
## Annotate all transcripts for gencode v25 based on the TxDb object
## we built previously.
ann_gencode_v25_hg38 <- annotateTranscripts(gencode_v25_hg38_txdb,
    annotationPackage = "org.Hs.eg.db",
    mappingInfo = list("column" = "ENTREZID", "keytype" = "ENSEMBL",
    "multiVals" = "first"))
    
# Error in .testForValidKeys(x, keys, keytype, fks) :
#   None of the keys entered are valid keys for 'ENSEMBL'. Please use the keys method to see a listing of valid arguments.
# > traceback()
# 10: stop(msg)
# 9: .testForValidKeys(x, keys, keytype, fks)
# 8: testSelectArgs(x, keys = keys, cols = cols, keytype = keytype,
#        fks = fks, skipValidKeysTest = skipValidKeysTest)
# 7: .select(x, keys, columns, keytype, jointype = jointype, ...)
# 6: select(x, keys = unique(keys), columns = column, keytype = keytype)
# 5: select(x, keys = unique(keys), columns = column, keytype = keytype)
# 4: mapIds_base(x, keys, column, keytype, ..., multiVals = multiVals)
# 3: mapIds(get(annotationPackage), keys = geneid, column = mappingInfo$column,
#        keytype = mappingInfo$keytype, multiVals = mappingInfo$multiVals)
# 2: mapIds(get(annotationPackage), keys = geneid, column = mappingInfo$column,
#        keytype = mappingInfo$keytype, multiVals = mappingInfo$multiVals)
# 1: annotateTranscripts(gencode_v25_hg38_txdb, annotationPackage = "org.Hs.eg.db",
#        mappingInfo = list(column = "ENTREZID", keytype = "ENSEMBL",
#            multiVals = "first"))


txdb = gencode_v25_hg38_txdb
tt <- transcriptsBy(txdb, by="gene")
RR <- ranges(tt)
geneid <- names(RR)
# > head(geneid)
# [1] "ENSG00000141956.13" "ENSG00000141959.16" "ENSG00000142149.8"  "ENSG00000142156.14" "ENSG00000142166.12"
# [6] "ENSG00000142168.14"

annotationPackage = "org.Hs.eg.db"
mappingInfo = list("column" = "ENTREZID", "keytype" = "ENSEMBL",
"multiVals" = "first")


geneid <- mapIds(get(annotationPackage), keys = geneid, column= mappingInfo$column, keytype = mappingInfo$keytype, multiVals = mappingInfo$multiVals)
# > traceback()
# 9: stop(msg)
# 8: .testForValidKeys(x, keys, keytype, fks)
# 7: testSelectArgs(x, keys = keys, cols = cols, keytype = keytype,
#        fks = fks, skipValidKeysTest = skipValidKeysTest)
# 6: .select(x, keys, columns, keytype, jointype = jointype, ...)
# 5: select(x, keys = unique(keys), columns = column, keytype = keytype)
# 4: select(x, keys = unique(keys), columns = column, keytype = keytype)
# 3: mapIds_base(x, keys, column, keytype, ..., multiVals = multiVals)
# 2: mapIds(get(annotationPackage), keys = geneid, column = mappingInfo$column,
#        keytype = mappingInfo$keytype, multiVals = mappingInfo$multiVals)
# 1: mapIds(get(annotationPackage), keys = geneid, column = mappingInfo$column,
#        keytype = mappingInfo$keytype, multiVals = mappingInfo$multiVals)


geneid2 <- gsub('\\..*', '', head(geneid))
mapIds(get(annotationPackage), keys = geneid2, column= mappingInfo$column, keytype = mappingInfo$keytype, multiVals = mappingInfo$multiVals)
# 'select()' returned 1:1 mapping between keys and columns
# ENSG00000141956 ENSG00000141959 ENSG00000142149 ENSG00000142156 ENSG00000142166 ENSG00000142168
#         "63977"          "5211"         "30811"          "1291"          "3454"          "6647"

small_reprex <- c("ENSG00000141956.13", "ENSG00000141959.16", "ENSG00000142149.8", "ENSG00000142156.14", "ENSG00000142166.12", "ENSG00000142168.14")
mappingInfo <- list("column" = "ENTREZID", "keytype" = "ENSEMBL", "multiVals" = "first")
mapIds(get('org.Hs.eg.db'), keys = small_reprex, column= mappingInfo$column, keytype = mappingInfo$keytype, multiVals = mappingInfo$multiVals)
# Error in .testForValidKeys(x, keys, keytype, fks) :
#   None of the keys entered are valid keys for 'ENSEMBL'. Please use the keys method to see a listing of valid arguments.
  
small_reprex2 <- gsub('\\..*', '', head(geneid))
mapIds(get('org.Hs.eg.db'), keys = small_reprex2, column= mappingInfo$column, keytype = mappingInfo$keytype, multiVals = mappingInfo$multiVals)
# 'select()' returned 1:1 mapping between keys and columns
# ENSG00000141956 ENSG00000141959 ENSG00000142149 ENSG00000142156 ENSG00000142166 ENSG00000142168
#         "63977"          "5211"         "30811"          "1291"          "3454"          "6647"


#### AHHHHH!!! Using the argument simplifyGeneID fixes this issue!
ann_gencode_v25_hg38 <- annotateTranscripts(gencode_v25_hg38_txdb,
    annotationPackage = "org.Hs.eg.db",
    mappingInfo = list("column" = "ENTREZID", "keytype" = "ENSEMBL",
    "multiVals" = "first"), simplifyGeneID = TRUE)
# Getting TSS and TSE.
# Getting CSS and CSE.
# Getting exons.
# Annotating genes.
# 'select()' returned 1:many mapping between keys and columns
        
library('sessioninfo')
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.5.2 Patched (2019-02-17 r76113)
#  os       macOS Mojave 10.14.3
#  system   x86_64, darwin15.6.0
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2019-02-27
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  acepack                1.4.1     2016-10-29 [1] CRAN (R 3.5.0)
#  AnnotationDbi        * 1.44.0    2018-10-30 [1] Bioconductor
#  assertthat             0.2.0     2017-04-11 [1] CRAN (R 3.5.0)
#  backports              1.1.3     2018-12-14 [1] CRAN (R 3.5.0)
#  base64enc              0.1-3     2015-07-28 [1] CRAN (R 3.5.0)
#  bibtex                 0.4.2     2017-06-30 [1] CRAN (R 3.5.0)
#  Biobase              * 2.42.0    2018-10-30 [1] Bioconductor
#  BiocGenerics         * 0.28.0    2018-10-30 [1] Bioconductor
#  BiocParallel         * 1.16.6    2019-02-10 [1] Bioconductor
#  biomaRt                2.38.0    2018-10-30 [1] Bioconductor
#  Biostrings             2.50.2    2019-01-03 [1] Bioconductor
#  bit                    1.1-14    2018-05-29 [1] CRAN (R 3.5.0)
#  bit64                  0.9-7     2017-05-08 [1] CRAN (R 3.5.0)
#  bitops                 1.0-6     2013-08-17 [1] CRAN (R 3.5.0)
#  blob                   1.1.1     2018-03-25 [1] CRAN (R 3.5.0)
#  BSgenome               1.50.0    2018-10-30 [1] Bioconductor
#  bumphunter           * 1.24.5    2018-12-01 [1] Bioconductor
#  checkmate              1.9.1     2019-01-15 [1] CRAN (R 3.5.2)
#  cli                    1.0.1     2018-09-25 [1] CRAN (R 3.5.0)
#  cluster                2.0.7-1   2018-04-13 [1] CRAN (R 3.5.2)
#  codetools              0.2-16    2018-12-24 [1] CRAN (R 3.5.2)
#  colorout             * 1.2-0     2019-02-18 [1] Github (jalvesaq/colorout@cc5fbfa)
#  colorspace             1.4-0     2019-01-13 [1] CRAN (R 3.5.2)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  data.table             1.12.0    2019-01-13 [1] CRAN (R 3.5.2)
#  DBI                    1.0.0     2018-05-02 [1] CRAN (R 3.5.0)
#  DelayedArray         * 0.8.0     2018-10-30 [1] Bioconductor
#  derfinder              1.16.1    2018-12-03 [1] Bioconductor
#  derfinderHelper        1.16.1    2018-12-03 [1] Bioconductor
#  digest                 0.6.18    2018-10-10 [1] CRAN (R 3.5.0)
#  doRNG                  1.7.1     2018-06-22 [1] CRAN (R 3.5.0)
#  downloader             0.4       2015-07-09 [1] CRAN (R 3.5.0)
#  dplyr                  0.8.0.1   2019-02-15 [1] CRAN (R 3.5.2)
#  foreach              * 1.4.4     2017-12-12 [1] CRAN (R 3.5.0)
#  foreign                0.8-71    2018-07-20 [1] CRAN (R 3.5.2)
#  Formula                1.2-3     2018-05-03 [1] CRAN (R 3.5.0)
#  GenomeInfoDb         * 1.18.2    2019-02-12 [1] Bioconductor
#  GenomeInfoDbData       1.2.0     2019-02-18 [1] Bioconductor
#  GenomicAlignments      1.18.1    2019-01-04 [1] Bioconductor
#  GenomicFeatures      * 1.34.3    2019-01-28 [1] Bioconductor
#  GenomicFiles           1.18.0    2018-10-30 [1] Bioconductor
#  GenomicRanges        * 1.34.0    2018-10-30 [1] Bioconductor
#  GEOquery               2.50.5    2018-12-22 [1] Bioconductor
#  ggplot2                3.1.0     2018-10-25 [1] CRAN (R 3.5.0)
#  glue                   1.3.0     2018-07-17 [1] CRAN (R 3.5.0)
#  gridExtra              2.3       2017-09-09 [1] CRAN (R 3.5.0)
#  gtable                 0.2.0     2016-02-26 [1] CRAN (R 3.5.0)
#  Hmisc                  4.2-0     2019-01-26 [1] CRAN (R 3.5.2)
#  hms                    0.4.2     2018-03-10 [1] CRAN (R 3.5.0)
#  htmlTable              1.13.1    2019-01-07 [1] CRAN (R 3.5.2)
#  htmltools              0.3.6     2017-04-28 [1] CRAN (R 3.5.0)
#  htmlwidgets            1.3       2018-09-30 [1] CRAN (R 3.5.0)
#  httr                   1.4.0     2018-12-11 [1] CRAN (R 3.5.0)
#  IRanges              * 2.16.0    2018-10-30 [1] Bioconductor
#  iterators            * 1.0.10    2018-07-13 [1] CRAN (R 3.5.0)
#  jsonlite               1.6       2018-12-07 [1] CRAN (R 3.5.0)
#  knitr                  1.21      2018-12-10 [1] CRAN (R 3.5.2)
#  lattice                0.20-38   2018-11-04 [1] CRAN (R 3.5.2)
#  latticeExtra           0.6-28    2016-02-09 [1] CRAN (R 3.5.0)
#  lazyeval               0.2.1     2017-10-29 [1] CRAN (R 3.5.0)
#  limma                  3.38.3    2018-12-02 [1] Bioconductor
#  locfit               * 1.5-9.1   2013-04-20 [1] CRAN (R 3.5.0)
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  Matrix                 1.2-15    2018-11-01 [1] CRAN (R 3.5.2)
#  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.5.0)
#  memoise                1.1.0     2017-04-21 [1] CRAN (R 3.5.0)
#  munsell                0.5.0     2018-06-12 [1] CRAN (R 3.5.0)
#  nnet                   7.3-12    2016-02-02 [1] CRAN (R 3.5.2)
#  org.Hs.eg.db         * 3.7.0     2019-02-18 [1] Bioconductor
#  pillar                 1.3.1     2018-12-15 [1] CRAN (R 3.5.0)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.0)
#  pkgmaker               0.27      2018-05-25 [1] CRAN (R 3.5.0)
#  plyr                   1.8.4     2016-06-08 [1] CRAN (R 3.5.0)
#  prettyunits            1.0.2     2015-07-13 [1] CRAN (R 3.5.0)
#  progress               1.2.0     2018-06-14 [1] CRAN (R 3.5.0)
#  purrr                  0.3.0     2019-01-27 [1] CRAN (R 3.5.2)
#  qvalue                 2.14.1    2019-01-10 [1] Bioconductor
#  R6                     2.4.0     2019-02-14 [1] CRAN (R 3.5.2)
#  RColorBrewer           1.1-2     2014-12-07 [1] CRAN (R 3.5.0)
#  Rcpp                   1.0.0     2018-11-07 [1] CRAN (R 3.5.0)
#  RCurl                  1.95-4.11 2018-07-15 [1] CRAN (R 3.5.0)
#  readr                  1.3.1     2018-12-21 [1] CRAN (R 3.5.0)
#  recount              * 1.8.2     2019-02-16 [1] Bioconductor
#  registry               0.5       2017-12-03 [1] CRAN (R 3.5.0)
#  rentrez                1.2.1     2018-03-05 [1] CRAN (R 3.5.0)
#  reshape2               1.4.3     2017-12-11 [1] CRAN (R 3.5.0)
#  rlang                  0.3.1     2019-01-08 [1] CRAN (R 3.5.2)
#  rngtools               1.3.1     2018-05-15 [1] CRAN (R 3.5.0)
#  rpart                  4.1-13    2018-02-23 [1] CRAN (R 3.5.2)
#  Rsamtools              1.34.1    2019-01-31 [1] Bioconductor
#  RSQLite                2.1.1     2018-05-06 [1] CRAN (R 3.5.0)
#  rstudioapi             0.9.0     2019-01-09 [1] CRAN (R 3.5.2)
#  rtracklayer          * 1.42.1    2018-11-21 [1] Bioconductor
#  S4Vectors            * 0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [1] CRAN (R 3.5.0)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.5.0)
#  stringi                1.3.1     2019-02-13 [1] CRAN (R 3.5.2)
#  stringr                1.4.0     2019-02-10 [1] CRAN (R 3.5.2)
#  SummarizedExperiment * 1.12.0    2018-10-30 [1] Bioconductor
#  survival               2.43-3    2018-11-26 [1] CRAN (R 3.5.2)
#  tibble                 2.0.1     2019-01-12 [1] CRAN (R 3.5.2)
#  tidyr                  0.8.2     2018-10-28 [1] CRAN (R 3.5.0)
#  tidyselect             0.2.5     2018-10-11 [1] CRAN (R 3.5.0)
#  VariantAnnotation      1.28.11   2019-02-18 [1] Bioconductor
#  withr                  2.1.2     2018-03-15 [1] CRAN (R 3.5.0)
#  xfun                   0.5       2019-02-20 [1] CRAN (R 3.5.2)
#  XML                    3.98-1.17 2019-02-08 [1] CRAN (R 3.5.2)
#  xml2                   1.2.0     2018-01-24 [1] CRAN (R 3.5.0)
#  xtable                 1.8-3     2018-08-29 [1] CRAN (R 3.5.0)
#  XVector                0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc               1.28.0    2018-10-30 [1] Bioconductor
#
# [1] /Library/Frameworks/R.framework/Versions/3.5/Resources/library
# >