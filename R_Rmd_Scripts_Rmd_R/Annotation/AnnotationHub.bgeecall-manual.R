## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE-------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("BgeeCall")

## ---- message = FALSE, warning = FALSE-----------------------------------
library(BgeeCall)

## ---- eval=FALSE---------------------------------------------------------
#  library("ShortRead")
#  # keep 48.000 reads
#  sampler <- FastqSampler(file.path("absolute_path","/SRX099901/SRR350955.fastq.gz"), 48000)
#  set.seed(1); SRR350955 <- yield(sampler)
#  writeFastq(object = SRR350955, file =file.path( "absolute_path","SRX099901_subset", "SRR350955_subset.fastq.gz"), mode = "w", full = FALSE, compress = TRUE)

## ---- message = FALSE, warning = FALSE-----------------------------------
ah <- AnnotationHub::AnnotationHub()
ah.Ensembl <- AnnotationHub::query(ah, c("Ensembl", "Homo sapiens"))
ah.UCSC <- AnnotationHub::query(ah, c("UCSC", "Homo sapiens"))
ah.UCSC <- AnnotationHub::query(ah, c("UCSC", "Homo sapiens", "GRanges"))
Encode.annotation.obj <- ah[[" AH5016"]]
transcriptome_object <- rtracklayer::import.2bit(ah_resources[["AH50453"]])

## ---- message = FALSE, warning = FALSE-----------------------------------
# create an object of class UserMetadata and specify the species ID
user_BgeeCall <- new("UserMetadata", species_id = "6239")
# import annotation and transcriptome in the user_BgeeCall object
# it is possible to import them using an S4 object (GRanges, DNAStringSet) or a file (gtf, fasta)
user_BgeeCall <- setAnnotationFromObject(user_BgeeCall, annotation_object, "WBcel235_84")
user_BgeeCall <- setTranscriptomeFromObject(user_BgeeCall, transcriptome_object, "WBcel235")
# provide path to the directory of your RNA-Seq library
user_BgeeCall <- setRNASeqLibPath(user_BgeeCall,
                                  system.file("extdata",
                                              "SRX099901_subset",
                                              package = "BgeeCall"))

## ---- eval = FALSE-------------------------------------------------------
#  calls_output <- generate_calls_workflow(userMetadata = user_BgeeCall)

## ---- echo=FALSE---------------------------------------------------------
user_BgeeCall <- setWorkingPath(user_BgeeCall, system.file("extdata", package = "BgeeCall"))
user_BgeeCall<- setSimpleArborescence(user_BgeeCall, TRUE)
calls_output <- generate_presence_absence(myUserMetadata = user_BgeeCall)

## ---- message = FALSE, warning = FALSE-----------------------------------
head.DataTable(x = read.table(calls_output$calls_tsv_path, header = TRUE), n = 5)

## ---- message = FALSE, warning = FALSE-----------------------------------
read.table(calls_output$cutoff_info_file_path)

## ---- message = FALSE, warning = FALSE-----------------------------------
head.DataTable(x = read.table(calls_output$abundance_tsv, header = TRUE), n = 5)
calls_output$TPM_distribution_path
calls_output$abundance_tsv

## ---- eval = FALSE-------------------------------------------------------
#  Biobase::openPDF(calls_output$TPM_distribution_path)

## ---- eval=FALSE---------------------------------------------------------
#  generate_calls_workflow(userFile = "path_to_your_file.tsv")

## ------------------------------------------------------------------------
list_bgee_species()

## ------------------------------------------------------------------------
list_intergenic_release()

## ------------------------------------------------------------------------
bgee <- new("BgeeMetadata", intergenic_release = "0.1")

## ---- eval=FALSE---------------------------------------------------------
#  kallisto <- new("KallistoMetadata", txOut = TRUE)
#  calls_output <- generate_calls_workflow(myAbundanceMetadata = kallisto, userMetadata = user_BgeeCall)

## ---- eval=FALSE---------------------------------------------------------
#  kallisto <- new("KallistoMetadata", install_kallisto = TRUE)
#  calls_output <- generate_calls_workflow(myAbundanceMetadata = kallisto, userMetadata = user_BgeeCall)

## ---- eval=FALSE---------------------------------------------------------
#  kallisto <- new("KallistoMetadata", single_end_parameters = "-t 3 --single -l 150 -s 30", pair_end_parameters = "-t 2 -b --seed 36")
#  calls_output <- generate_calls_workflow(myAbundanceMetadata = kallisto, userMetadata = user_BgeeCall)

## ---- eval=FALSE---------------------------------------------------------
#  # libraries with reads smaller than 70nt will use the index with kmer size = 21
#  kallisto <- new("KallistoMetadata", read_size_kmer_threshold = 70)
#  calls_output <- generate_calls_workflow(myAbundanceMetadata = kallisto, userMetadata = user_BgeeCall)

## ---- eval=FALSE---------------------------------------------------------
#  # RNA-Seq run SRR350955_subsetof from the RNA-Seq library will be used to generate the calls
#  user_BgeeCall <- setRunIds(user_BgeeCall, c("SRR350955_subset"))
#  calls_output <- run_from_object(myUserMetadata = user_BgeeCall)

## ------------------------------------------------------------------------
kallisto <- new("KallistoMetadata", cutoff = 0.1)

## ---- eval=FALSE---------------------------------------------------------
#  user_BgeeCall <- setRunIds(user_BgeeCall, "")
#  user_BgeeCall <- setSimpleArborescence(user_BgeeCall, TRUE)
#  calls_output <- run_from_object(myUserMetadata = user_BgeeCall)

## ----sessioninfo---------------------------------------------------------
sessionInfo()

## ----cleanup_after, echo=FALSE, message=FALSE, warning=FALSE-------------
unlink(BgeeCall:::get_kallisto_dir_path(kallisto, user_BgeeCall), recursive = TRUE)
unlink(file.path(getWorkingPath(user_BgeeCall), paste0(getIntergenicPrefix(bgee), "*")), recursive = TRUE)

