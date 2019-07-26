### R code from vignette source 'SplicingGraphs.Rnw'
#source("https://bioconductor.org/biocLite.R")
#biocLite("SplicingGraphs")
vignette("SplicingGraphs")
###################################################
### code chunk number 1: settings
###################################################
library(SplicingGraphs)
# sg <- SplicingGraphs(txdb)
# sg
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(RNAseqData.HNRNPC.bam.chr14)
isActiveSeq(txdb)[1:25]
txdb
save.image("~/Rscripts/splicegraphs_env.RData")
# splicing patterns of hg19 genome
# 
.precomputed_results_path <- "precomputed_results"
.loadPrecomputed <- function(objname)
{
    filename <- paste0(objname, ".rda")
    path <- file.path(.precomputed_results_path, filename)
    tempenv <- new.env(parent=emptyenv())
    load(path, envir=tempenv)
    get(objname, envir=tempenv)
}

.loadEnv<-function(objname)
{
    #filename <- paste0(objname, ".rda")
	filename <- paste0(objname, ".RData")
    path <- file.path(filename)
    tempenv <- new.env(parent=emptyenv())
    load(path, envir=tempenv)
    get(objname, envir=tempenv)
}
###################################################
### code chunk number 4: keep_chr14_only
###################################################
isActiveSeq(txdb)[-match("chr14", names(isActiveSeq(txdb)))] <- FALSE
names(which(isActiveSeq(txdb)))

###################################################
### code chunk number 6: load_sg_with_bubbles
###################################################
## We load a precomputed 'sg' that contains all the bubbles in
sg@.bubbles_cache
sg_with_bubbles <- function(x) (x = .loadPrecomputed(sg@.bubbles_cache))
sg_with_bubbles <- .loadPrecomputed("sg_with_bubbles")
sg_env<-.loadEnv("~/Rscripts/splicegraphs_env")
sg@.bubbles_cache <- sg_with_bubbles@.bubbles_cache
## Replace NAs with FALSE in circularity flag (because having the flag set
## to NA instead of FALSE (or vice-versa) is not considered a significant
## difference between the 2 objects).
isCircular(sg) <- isCircular(sg) %in% TRUE
isCircular(sg_with_bubbles) <- isCircular(sg_with_bubbles) %in% TRUE
if (!identical(sg, sg_with_bubbles))
    stop("'sg' is not identical to precomputed version")


###################################################
### code chunk number 7: names_sg
###################################################
names(sg)[1:20]
seqnames(sg)[1:20]
strand(sg)[1:20]
table(strand(sg))
sg[1:100][[1]]
elementNROWS(sg)[1:20]
sg[["3183"]]
mcols(sg[["3183"]])
mcols(sg[["3183"]])$txpath
txpath(sg[["3183"]])
plotTranscripts(sg[["3183"]],)
# pdf("3183-transcripts.pdf", width=6, height=3)
plotTranscripts(sg[["3183"]])
###################################################
### code chunk number 16: unlist_sg
###################################################
ex_by_tx <- unlist(sg)
head(names(ex_by_tx))
ex_by_tx[names(ex_by_tx) %in% c("10001", "100129075")]
sg[strand(sg) == "-"]
sg[1:20]
tail(sg)  # equivalent to 'sg[tail(seq_along(sg))]'
sgTry<-sg["3183"]
sgedges(sg["3183"])
sgnodes(sg["3183"])
edges_by_gene <- sgedgesByGene(sg)
# sgraph(x = ,keep.dup.edges = ,tx_id.as.edge.label = ,as.igraph = ,
edges_by_gene[["3183"]]
plot(sg["3183"])
plot(sgraph(sg["3183"], tx_id.as.edge.label=TRUE))
# pdf("3183-sgraph.pdf", width=3, height=5)
plot(sgraph(sg["3183"]))
plot(sgraph(sg["100309464"], tx_id.as.edge.label=TRUE))
dev.off()
###################################################
### code chunk number 26: bubbles_sg_100309464
###################################################
bubbles(sg["100309464"])
codes <- bubbles(sg["10202"])$AScode
code.df<-data.frame(AScode=codes, description=ASCODE2DESC[codes], row.names=NULL)
head(code.df)

AScode_list <- lapply(seq_along(sgTry), function(i) bubbles(sgTry[i])$AScode)
names(AScode_list) <- names(sgTry)
AScode_table <- table(unlist(AScode_list))
AScode_table <- sort(AScode_table, decreasing=TRUE)
AScode_summary <- data.frame(AScode=names(AScode_table),
                             NbOfEvents=as.vector(AScode_table),
                             Desciption=ASCODE2DESC[names(AScode_table)])
nrow(AScode_summary)
head(AScode_summary, n=10)

###################################################
### code chunk number 30: nb_bubbles_per_gene
###################################################
nb_bubbles_per_gene <- elementNROWS(AScode_list)
head(sort(nb_bubbles_per_gene, decreasing=TRUE))
nb_unique_bubbles_per_gene<-elementNROWS(unique(CharacterList(AScode_list)))
head(sort(nb_unique_bubbles_per_gene, decreasing=TRUE))

###################################################
### code chunk number 34: bam_files
###################################################
# library(RNAseqData.HNRNPC.bam.chr14)
bam_files <- RNAseqData.HNRNPC.bam.chr14_BAMFILES
names(bam_files)  # the names of the runs
quickBamFlagSummary(bam_files[1], main.groups.only=TRUE)
###################################################
### code chunk number 36: readGAlignmentPairs
###################################################
param <- ScanBamParam(which=GRanges("chr14", IRanges(1, 20000000)))
galp <- readGAlignmentPairs(bam_files[1], param=param)
length(galp)  # nb of alignment pairs
galp
###################################################
### code chunk number 37: ScanBamParam
###################################################
flag0 <- scanBamFlag(isSecondaryAlignment=FALSE,
                     isNotPassingQualityControls=FALSE,
                     isDuplicate=FALSE)
param0 <- ScanBamParam(flag=flag0)
###################################################
### code chunk number 38: assignReads (eval = FALSE)
###################################################
## ## The following loop takes about 7 minutes on a modern laptop/desktop...
## for (i in seq_along(bam_files)) {
##     bam_file <- bam_files[i]
##     cat("Processing run ", names(bam_file), " ... ", sep="")
##     galp <- readGAlignmentPairs(bam_file, use.names=TRUE, param=param0)
##     sg <- assignReads(sg, galp, sample.name=names(bam_file))
##     cat("OK\n")
## }
###################################################
### code chunk number 39: load_sg_with_reads
###################################################
sg_with_reads <- .loadPrecomputed("sg_with_reads")
## Remove the reads from 'sg_with_reads' and compare with 'sg'.
if (!isTRUE(all.equal(removeReads(sg_with_reads), sg)))
    stop("after removal of the hits, precomputed version of 'sg_with_reads' ",
         "is not identical to 'sg'")
sg <- sg_with_reads

###################################################
### code chunk number 40: sg_3183_assignments
###################################################
edges_by_tx <- sgedgesByTranscript(sg["3183"], with.hits.mcols=TRUE)
edge_data <- mcols(unlist(edges_by_tx))
colnames(edge_data)
head(edge_data[ , c("sgedge_id", "ERR127306.hits")])

###################################################
### code chunk number 42: load_reads_in_3183_region
###################################################
param <- ScanBamParam(flag=flag0, which=range(unlist(sg[["3183"]])))
reads <- readGAlignmentPairs(bam_files[1], use.names=TRUE, param=param)
junction_reads <- reads[njunc(first(reads)) + njunc(last(reads)) != 0L]
plotTranscripts(sg[["3183"]], reads=junction_reads, from=21675000, to=21702000)
#pdf("3183-transcripts-and-reads.pdf", width=12, height=12)
plotTranscripts(sg[["3183"]], reads=junction_reads, from=21698400, to=21698600)
#pdf("3183-transcripts-and-reads-zoom.pdf", width=12, height=12)
plotTranscripts(sg[["3183"]], reads=junction_reads, from=21698400, to=21698600)
dev.off()

###################################################
### code chunk number 47: countReads_sg
###################################################
sg_counts <- countReads(sg)
dim(sg_counts)
head(sg_counts[1:5])
sapply(sg_counts[-(1:2)], sum)


