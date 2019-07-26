## ---- echo=FALSE---------------------------------------------------------
library(knitr)
opts_chunk$set(cache=FALSE,
               error=FALSE)
library(alpineData)
dir <- system.file("extdata",package="alpineData")
metadata <- read.csv(file.path(dir,".alpineData/metadata.csv"),
                     stringsAsFactors=FALSE)
metadata[,c("Title","Performer","Date","Population")]
## ----message=FALSE-------------------------------------------------------
dir <- system.file("extdata",package="alpineData")
metadata <- read.csv(file.path(dir,"metadata.csv"),
                     stringsAsFactors=FALSE)
metadata[,c("Title","Performer","Date","Population")]
source("http://bioconductor.org/biocLite.R")
biocLite("polyester")
biocLite("alpine")
library(alpine)
library(polyester)
library(GenomicAlignments)
ERR188297()

## ----message=FALSE-------------------------------------------------------
library(rtracklayer)
dir <- system.file(package="alpineData", "extdata")
for (sample.name in metadata$Title) {
  # the reads are accessed with functions named
  # after the sample name. the following line calls
  # the function with the sample name and saves 
  # the reads to `gap`
  gap <- match.fun(sample.name)()
  file.name <- file.path(dir,paste0(sample.name,".bam"))
  export(gap, con=file.name)
}
bam.files <- file.path(dir, paste0(metadata$Title, ".bam"))
names(bam.files) <- metadata$Title
stopifnot(all(file.exists(bam.files)))

## ---- eval=FALSE---------------------------------------------------------
## library(ensembldb)
## gtf.file <- "Homo_sapiens.GRCh38.84.gtf"
## txdb <- EnsDb(gtf.file) # already an EnsDb
## txdf <- transcripts(txdb, return.type="DataFrame")
## tab <- table(txdf$gene_id)
## one.iso.genes <- names(tab)[tab == 1]
## # pre-selected genes based on medium to high counts
## # calculated using Rsubread::featureCounts
## selected.genes <- scan("selected.genes.txt", what="char")
## one.iso.txs <- txdf$tx_id[txdf$gene_id %in%
##                           intersect(one.iso.genes, selected.genes)]
## ebt0 <- exonsBy(txdb, by="tx")
## ebt.fit <- ebt0[one.iso.txs]

## ----message=FALSE-------------------------------------------------------
library(GenomicRanges)

## ------------------------------------------------------------------------
library(alpine)
data(preprocessedData)
# filter small genes and long genes
min.bp <- 600
max.bp <- 7000 
gene.lengths <- sum(width(ebt.fit))
summary(gene.lengths)
ebt.fit <- ebt.fit[gene.lengths > min.bp & gene.lengths < max.bp]
length(ebt.fit)
set.seed(1)
# better to use ~100 genes
ebt.fit <- ebt.fit[sample(length(ebt.fit),10)] 

## ------------------------------------------------------------------------
w <- getFragmentWidths(bam.files[1], ebt.fit[[1]])
c(summary(w), Number=length(w))
quantile(w, c(.025, .975))

## ------------------------------------------------------------------------
getReadLength(bam.files)

## ----message=FALSE-------------------------------------------------------
library(alpine)
library(BSgenome.Hsapiens.NCBI.GRCh38)
readlength <- 75 
minsize <- 125 # better 80 for this data
maxsize <- 175 # better 350 for this data
gene.names <- names(ebt.fit)
names(gene.names) <- gene.names

## ----buildFragtype-------------------------------------------------------
system.time({
fragtypes <- lapply(gene.names, function(gene.name) {
                      buildFragtypes(exons=ebt.fit[[gene.name]],
                                     genome=Hsapiens,
                                     readlength=readlength,
                                     minsize=minsize,
                                     maxsize=maxsize,
                                     gc.str=FALSE)
                    })
})
print(object.size(fragtypes), units="auto")

## ------------------------------------------------------------------------
head(fragtypes[[1]], 3)

## ------------------------------------------------------------------------
models <- list(
  "GC" = list(
    formula = "count ~ ns(gc,knots=gc.knots,Boundary.knots=gc.bk) + ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk) + gene",
    offset=c("fraglen")
  ),
  "all" = list(
    formula = "count ~ ns(gc,knots=gc.knots,Boundary.knots=gc.bk) + ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk) + gene",
    offset=c("fraglen","vlmm")
  )
)

## ----fitBiasModels-------------------------------------------------------
system.time({
fitpar <- lapply(bam.files, function(bf) {
                   fitBiasModels(genes=ebt.fit,
                                 bam.file=bf,
                                 fragtypes=fragtypes,
                                 genome=Hsapiens,
                                 models=models,
                                 readlength=readlength,
                                 minsize=minsize,
                                 maxsize=maxsize)
                 })
})
# this object saved was 'fitpar.small' for examples in alpine man pages
# fitpar.small <- fitpar 

## ------------------------------------------------------------------------
library(RColorBrewer)
palette(brewer.pal(8,"Dark2"))

## ----fraglen-------------------------------------------------------------
perf <- as.integer(factor(metadata$Performer))
plotFragLen(fitpar, col=perf)

## ----gccurve-------------------------------------------------------------
plotGC(fitpar, model="all", col=perf)

## ----relpos--------------------------------------------------------------
plotRelPos(fitpar, model="all", col=perf)

## ----vlmm----------------------------------------------------------------
plotOrder0(fitpar[["ERR188297"]][["vlmm.fivep"]][["order0"]])
plotOrder0(fitpar[["ERR188297"]][["vlmm.threep"]][["order0"]])

## ------------------------------------------------------------------------
print(head(fitpar[["ERR188297"]][["summary"]][["all"]]), row.names=FALSE)

## ---- eval=FALSE---------------------------------------------------------
## one.iso.genes <- intersect(names(tab)[tab == 1], selected.genes)
## two.iso.genes <- intersect(names(tab)[tab == 2], selected.genes)
## three.iso.genes <- intersect(names(tab)[tab == 3], selected.genes)
## set.seed(1)
## genes.theta <- c(sample(one.iso.genes, 2),
##                  sample(two.iso.genes, 2),
##                  sample(three.iso.genes, 2))
## txdf.theta <- txdf[txdf$gene_id %in% genes.theta,]
## ebt.theta <- ebt0[txdf.theta$tx_id]

## ------------------------------------------------------------------------
model.names <- c("null","fraglen.vlmm","GC")

## ----estimateAbundance---------------------------------------------------
system.time({
res <- lapply(genes.theta, function(gene.name) {
         txs <- txdf.theta$tx_id[txdf.theta$gene_id == gene.name]
         estimateAbundance(transcripts=ebt.theta[txs],
                           bam.files=bam.files,
                           fitpar=fitpar,
                           genome=Hsapiens,
                           model.names=model.names)
       })
})

## ------------------------------------------------------------------------
res[[1]][["ERR188297"]][["GC"]]
res[[6]][["ERR188297"]][["GC"]]

## ------------------------------------------------------------------------
mat <- extractAlpine(res, model="GC")
mat

## ------------------------------------------------------------------------
se <- extractAlpine(res, model="GC", transcripts=ebt.theta)
se

## ---- eval=FALSE---------------------------------------------------------
## norm.mat <- normalizeDESeq(mat, cutoff=0.1)

## ------------------------------------------------------------------------
data(preprocessedData)
prob.mat <- plotGC(fitpar, "all", return.type=2)
head(prob.mat)

## ------------------------------------------------------------------------
model.names <- c("fraglen","fraglen.vlmm","GC","all")

## ------------------------------------------------------------------------
fitpar[[1]][["model.params"]][c("minsize","maxsize")]

## ----predictCoverage-----------------------------------------------------
system.time({
  pred.cov <- predictCoverage(gene=ebt.fit[["ENST00000245479"]],
                              bam.files=bam.files["ERR188204"],
                              fitpar=fitpar,
                              genome=Hsapiens,
                              model.names=model.names)
})

## ------------------------------------------------------------------------
palette(brewer.pal(9, "Set1"))
frag.cov <- pred.cov[["ERR188204"]][["frag.cov"]]
plot(frag.cov, type="l", lwd=3, ylim=c(0,max(frag.cov)*1.5))
for (i in seq_along(model.names)) {
  m <- model.names[i]
  pred <- pred.cov[["ERR188204"]][["pred.cov"]][[m]]
  lines(pred, col=i, lwd=3)
}
legend("topright", legend=c("observed",model.names),
       col=c("black",seq_along(model.names)), lwd=3)

## ------------------------------------------------------------------------
sessionInfo()

