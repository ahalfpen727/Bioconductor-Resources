# makeSynMap.R
#
# Purpose:  Prepare a list of synonyms for HUGO recognized symbols of genes
#           with a protein product, rRNA, tRNA, and vault RNA, T-Cell receptors,
#           immunoglobulin genes, and protocadherin
#           cf. https://www.genenames.org/cgi-bin/statistics
#
#           Result: inst/extdata/synMap.RData is a character vector
#           containing synonyms for the above HUGO symbols that can be
#           uniquely mapped. These are either old names, or recognized
#           aliases.
#
# Version:  1.0
# Date:     2018 01 30
# Author:   Boris Steipe <boris.steipe@utoronto.ca>
#
# Version history:
#           1.0  production dataset, January 2018
#
# ToDo:
# Notes:
#
# ==============================================================================


# ====  PARAMETERS  ============================================================

BASEURL <- "ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_types/"

myHUGOfiles <- data.frame(tags = c("ORF", "rRNA", "tRNA", "vRNA",
                                   "TCR", "IG", "PCDH"),
                          fNames = c("gene_with_protein_product.txt",
                                     "RNA_ribosomal.txt",
                                     "RNA_transfer.txt",
                                     "RNA_vault.txt",
                                     "T_cell_receptor_gene.txt",
                                     "immunoglobulin_gene.txt",
                                     "protocadherin.txt"),
                          stringsAsFactors = FALSE)


# ====  PACKAGES  ==============================================================

if (! require(readr, quietly=TRUE)) {
  install.packages("readr")
  library(readr)
}

# ====  PROCESS  ===============================================================
#


rawData <- data.frame(symbols = character(),
                      alias = character(),
                      previous = character())

# ftp all required datasets and keep symbols, aliases and previous names
for (fName in myHUGOfiles$fNames) {
  tmp <- read_tsv(paste0(BASEURL, fName))
  rawData <- rbind(rawData, data.frame(symbols = tmp$symbol,
                                       alias = tmp$alias_symbol,
                                       previous = tmp$prev_symbol,
                                       stringsAsFactors = FALSE))
}
rm(tmp)


# We don't distinguish between alias and previous symbol, so we
# stack aliases and previous names into one column, order by symbol
rawData <- data.frame(symbols = c(rawData$symbols, rawData$symbols),
                      synonyms = c(rawData$alias, rawData$previous),
                      stringsAsFactors = FALSE)
rawData <- rawData[order(rawData$symbols), ]   # 40696 rows

# preprocess to have one row per symbol, containing the symbol
# and all its synonyms

# synSrc holds this data
synSrc <- data.frame(symbols = unique(rawData$symbols),
                     synonyms = NA,
                     stringsAsFactors = FALSE)

# fill synSrc symbol by symbol, synonyms go in string separated by "|"
l <- length(synSrc$symbols)
for (i in 1:l) {
  pBar(i, l)  # progress bar. See ./R/pBar.R
  sym <- synSrc$symbols[i]
  synSrc$synonyms[i] <- paste0(rawData$synonyms[rawData$symbols == sym],
                           collapse = "|")
}

head(synSrc)   # 20348 rows

# remove all NAs and all synonyms that are empty

synSrc$synonyms <- gsub("\\|NA\\|", "|", synSrc$synonyms)
synSrc$synonyms <- gsub("^NA\\|", "", synSrc$synonyms)
synSrc$synonyms <- gsub("\\|NA$", "", synSrc$synonyms)

sel <- synSrc$synonym != "NA" &
  synSrc$synonym != "" &
  ! is.na(synSrc$synonym)
synSrc <- synSrc[sel, ]

head(synSrc)   # 17445 rows

# remove all synonyms that are outright duplicated
x <- unique(synSrc$synonyms[duplicated(synSrc$synonyms)])
#    sanity check
for (i in seq_along(x)) {
  print(synSrc[synSrc$synonyms == x[i], ])
}
#    remove
synSrc <- synSrc[!(synSrc$synonyms %in% x) , ]
head(synSrc)   # 17391 rows

# build the synonym map
synMap <- data.frame(synonyms = character(),
                     symbols = character(),
                     stringsAsFactors = FALSE)


# for each symbol, strsplit() the synonyms, fill a dataframe with symbol and
# synonyms, and rbind() it to synMap
l <- length(synSrc$symbols)
for (i in 1:l) {
  pBar(i, l)
  sym <- synSrc$symbols[i]
  synMap <- rbind(synMap,
                  data.frame(symbols = sym,
                             synonyms = strsplit(synSrc$synonyms[i],
                                                 "\\|")[[1]],
                             stringsAsFactors = FALSE))
}

head(synMap, 20)   # 43252 rows

# Duplicated synonyms that map to different symbols must be removed.
# First: are there any duplicated symbols that map to the same symbol? That
# would be oK:
x <- paste(synMap$symbols, synMap$synonyms, sep = "|")
any(duplicated(x)) # FALSE ... no. Therefore all remaining duplicated
                   # synonyms map to different symbols. They need to be removed.
                   # How many?

sum(duplicated(synMap$synonyms))  #1437

# fetch them from the synonyms column
x <- unique(synMap$synonyms[duplicated(synMap$synonyms)])

# confirm and sanity check - choose N duplicates at random
N <- 10
for (i in sample(1:length(x), N)) {
  print(synMap[synMap$synonyms == x[i], ])
} # indeed: smae synonym, different symbol

# remove them all
synMap <- synMap[!(synMap$synonyms %in% x) , ]
head(synMap)   # 40647 rows

# final sanity check: are there symbols in the synonym list?
sum(unique(synMap$symbols) %in% synMap$synonyms)  #  461 Ooops

# fetch them from the synonyms column
x <- unique(synMap$symbols[synMap$symbols %in% synMap$synonyms])  # 461

# inspect some of them and assess on https://www.genenames.org/
synMap[synMap$symbols == x[1], ]
synMap[synMap$synonyms == x[1], ]

synMap[synMap$symbols == x[20], ]
synMap[synMap$synonyms == x[20], ]

synMap[synMap$symbols == x[30], ]
synMap[synMap$synonyms == x[30], ]

# Messy. It seems that these are really identical identifiers for different
# genes. Let's define that synonyms that are identical to a current gene symbol
# shoud not appear in our map.
synMap <- synMap[!(synMap$synonyms %in% x) , ]
head(synMap)   # 40186 rows

#Final confirmation
any(duplicated(synMap$synonyms))              # Must be FALSE
any((synMap$symbols %in% synMap$synonyms))    # Must be FALSE
str(synMap)                                   # Must be character vectors
nrow(synMap)                                  # 40186
length(unique(synMap$symbols))                # 17090

# save
save(synMap, file = "inst/extdata/synMap.RData")


# [END]


