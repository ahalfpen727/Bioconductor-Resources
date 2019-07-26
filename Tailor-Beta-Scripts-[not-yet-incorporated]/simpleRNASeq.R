## ----style, echo=FALSE, results='asis'-----------------------------------
BiocStyle::markdown()

## ----message=FALSE,warning=FALSE,results='hide',echo=FALSE---------------
options(digits=2)

## ----library-------------------------------------------------------------
library(easyRNASeq)

## ----P trichocarpa annotation--------------------------------------------
library(curl)
curl_download(url=paste0("ftp://ftp.plantgenie.org/Data/PopGenIE/",
                         "Populus_trichocarpa/v3.0/v10.1/GFF3/",
                         "Ptrichocarpa_210_v3.0_gene_exons.gff3.gz"),
                  destfile=,"./Ptrichocarpa_210_v3.0_gene_exons.gff3.gz")

## ----AnnotParam----------------------------------------------------------
    annotParam <- AnnotParam(
        datasource="./Ptrichocarpa_210_v3.0_gene_exons.gff3.gz")

## ----create synthetic transcripts----------------------------------------
annotParam <- createSyntheticTranscripts(annotParam,verbose=FALSE)

## ----save the object-----------------------------------------------------
save(annotParam,
file="./Ptrichocarpa_210_v3.0_gene_exons_synthetic-transcripts_annotParam.rda")

## ----create synthetic transcripts as gI----------------------------------
gI <- createSyntheticTranscripts(
    "./Ptrichocarpa_210_v3.0_gene_exons.gff3.gz",
    verbose=FALSE)

## ----export the file-----------------------------------------------------
writeGff3(gI,
          file="./Ptrichocarpa_210_v3.0_gene_exons_synthetic-transcripts.gff3.gz")

## ----bam files-----------------------------------------------------------
library(curl)
curl_download(url=paste0("ftp://ftp.plantgenie.org/Tutorials/RnaSeqTutorial/",
                         "data/star/md5.txt"),
                  destfile="md5.txt")

## ----data----------------------------------------------------------------
data(RobinsonDelhomme2014)
lapply(RobinsonDelhomme2014[1:6,"Filename"],function(f){
    curl_download(url=paste0("ftp://ftp.plantgenie.org/Tutorials/",
                             "RnaSeqTutorial/data/star/",f),
                  destfile=f)
})

## ----bamParam------------------------------------------------------------
bamParam <- BamParam(paired = TRUE,
                     stranded = FALSE)

## ----rnaSeqParam---------------------------------------------------------
rnaSeqParam <- RnaSeqParam(annotParam = annotParam,
                           bamParam = bamParam,
                           countBy = "genes",
                           precision = "read")

## ----session info, echo=FALSE--------------------------------------------
sessionInfo()

## ----cleanup, echo=FALSE-------------------------------------------------
    file.remove(c(
        "./Ptrichocarpa_210_v3.0_gene_exons.gff3.gz",
        "./Ptrichocarpa_210_v3.0_gene_exons_synthetic-transcripts_annotParam.rda",
        "./Ptrichocarpa_210_v3.0_gene_exons_synthetic-transcripts.gff3.gz"))

