###################################################
### chunk number 2: initialize
###################################################
library(GenomicRanges)
library(Gviz)
library(GenomicRanges)
> data(cpgIslands)
> class(cpgIslands)
##################################################
### chunk number 1: options
###################################################
library(Rsamtools)
testFile <- system.file("bam", "isowt5_13e.bam", package="leeBamViews")
aligns <- readBamGappedAlignments(testFile)


###################################################
### chunk number 20: GetAnnoations
###################################################
library(GenomicFeatures)
txdb <- makeTranscriptDbFromUCSC(genome="sacCer2", tablename="sgdGene")


###################################################
### chunk number 21: exonsBy
###################################################
exonRanges <- exonsBy(txdb, "tx")
length(exonRanges)
exonRanges[1]


###################################################
### chunk number 22: ammendData
###################################################
levels(rname(aligns))
levels(rname(aligns)) <-
   c(paste("chr", as.roman(1:16), sep=""), "chrM")


###################################################
### chunk number 23: count
###################################################
counts <- countOverlaps(exonRanges, aligns)


###################################################
### chunk number 24: numBases
###################################################
numBases <- sum(width(exonRanges))
geneLengthsInKB <- numBases / 1000


###################################################
### chunk number 25: millionsMapped
###################################################
millionsMapped <- sum(counts) / 1000000


###################################################
### chunk number 26: RPM
###################################################
# counted reads / total reads in millions
rpm <- counts / millionsMapped


###################################################
### chunk number 27: RPKM
###################################################
# reads per million per geneLength in Kb
rpkm <- rpm / geneLengthsInKB


###################################################
### chunk number 28: sort
###################################################
sortedRPKM <- sort(rpkm)
highScoreGenes <- tail(sortedRPKM)


###################################################
### chunk number 29: annotate1
###################################################
txs <- transcripts(txdb,
                   vals=list(tx_id=names(highScoreGenes)),
                   columns=c("tx_id","gene_id"))
systNames <- as.vector(unlist(elementMetadata(txs)["gene_id"]))


###################################################
### chunk number 30: annnotate2
###################################################
library(org.Sc.sgd.db)
toTable(org.Sc.sgdGENENAME[systNames])


###################################################
### chunk number 31: filterKnowns
###################################################
filtData <- subsetByOverlaps(aligns, exonRanges)
length(filtData)


###################################################
### chunk number 32: filterKnowns2
###################################################
filtData2 <- subsetByOverlaps(aligns, transcriptsBy(txdb, "gene"))
length(filtData2)


###################################################
### chunk number 33: coverage
###################################################
cov <- coverage(filtData)


###################################################
### chunk number 34: coverageSubset
###################################################
cov <- cov[13]


###################################################
### chunk number 35: islands
###################################################
islands <- slice(cov, lower = 1)


###################################################
### chunk number 36: continous >1000
###################################################
transcribedRegions <- islands[width(islands) > 1000]
txr <- islands[width(islands) > 1000]


###################################################
### chunk number 37: getSequence
###################################################
library(BSgenome.Scerevisiae.UCSC.sacCer2)
getYeastSequence <- function(data) {
   chr <- rep("chrXIII",length(start(data[[1]])))
   starts <- start(data[[1]])
   ends <- end(data[[1]])
   strands <- rep("+",length(start(data[[1]])))
   getSeq(Scerevisiae, names=chr, start=starts, end=ends, strand=strands)
}
DNASet <- DNAStringSet(getYeastSequence(transcribedRegions))

###################################################
### chunk number 3: EatonEtAlChIPseq
###################################################
library(EatonEtAlChIPseq)
data(orcAlignsRep1)
orcAlignsRep1


###################################################
### chunk number 4: EatonAlignmentFiltering-1
###################################################
subsetRep1 <-
   orcAlignsRep1[alignData(orcAlignsRep1)[["nMismatchBestHit"]] <= 3L]
length(subsetRep1) / length(orcAlignsRep1)
subsetRep1 <- subsetRep1[occurrenceFilter(withSread=FALSE)(subsetRep1)]
length(subsetRep1) / length(orcAlignsRep1)
subsetRep1



###################################################
### chunk number 5: EatonOrcRanges
###################################################
rangesRep1 <- as(subsetRep1, "GRanges")
head(rangesRep1, 3)


###################################################
### chunk number 6: EatonAddSeqlengths
###################################################
seqlengths(rangesRep1) <- 784333


###################################################
### chunk number 7: EatonNegStarts
###################################################
negRangesRep1 <- rangesRep1[strand(rangesRep1) == "-"]
negStartsRep1 <- resize(negRangesRep1, 1)


###################################################
### chunk number 8: EatonPossibleEnds
###################################################
posRangesRep1 <- rangesRep1[strand(rangesRep1) == "+"]
posEndsRep1 <- shift(posRangesRep1, 99)
posEndsRep1 <- resize(posEndsRep1, 100)
strand(posEndsRep1) <- "-"


###################################################
### chunk number 9: EatonMatchingPairs
###################################################
strandMatching <- findOverlaps(negStartsRep1, posEndsRep1)
posKeep <- unique(subjectHits(strandMatching))
negKeep <- unique(queryHits(strandMatching))
length(posKeep) / length(posEndsRep1)
length(negKeep) / length(negStartsRep1)
(length(posKeep) + length(negKeep)) / length(orcAlignsRep1)


###################################################
### chunk number 10: EatonAlignmentFiltering-2
###################################################
posFilteredRangesRep1 <- posRangesRep1[posKeep]
negFilteredRangesRep1 <- negRangesRep1[negKeep]


###################################################
### chunk number 11: CoverageWeights
###################################################
posWeights <- c(seq(0.01, 1, length = 100), rep(c(1, 0), c(101, 200)))
negWeights <- rev(posWeights)
plot(-200:200, posWeights, xlab = "Relative Position",
     ylab = "Coverage Weight", type = "l")


###################################################
### chunk number 12: EatonStartCoverage
###################################################
posStartsCoverRep1 <- coverage(resize(posFilteredRangesRep1, 1))
negStartsCoverRep1 <- coverage(resize(negFilteredRangesRep1, 1))


###################################################
### chunk number 13: EatonExtendedCoverage
###################################################
posExtCoverRep1 <-
  round(runwtsum(posStartsCoverRep1, k = 401, wt = posWeights,
                 endrule = "constant"))
negExtCoverRep1 <-
  round(runwtsum(negStartsCoverRep1, k = 401, wt = negWeights,
                 endrule = "constant"))


###################################################
### chunk number 14: plotCoverage
###################################################
plotCoverage <-
function(x, xlab = "Position", ylab = "Coverage",...)
{
    plot(c(start(x), length(x)), c(runValue(x), tail(runValue(x), 1)),
         type = "s", col = "blue", xlab = xlab, ylab = ylab, ...)
}
plotStrandedCoverage <-
function(positive, negative, xlab = "Position", ylab = "Coverage",...)
{
    ylim <- min(max(positive), max(negative)) * c(-1, 1)
    plotCoverage(positive, ylim = ylim, ...)
    lines(c(start(negative), length(negative)),
          - c(runValue(negative), tail(runValue(negative), 1)),
          type = "s", col = "red")
    abline(h = 0, col = "dimgray")
}


###################################################
### chunk number 15: EatonPlotStrands
###################################################
plotStrandedCoverage(posExtCoverRep1[[1]], negExtCoverRep1[[1]])


###################################################
### chunk number 16: EatonExtendedCoverage
###################################################
combExtCoverRep1 <- pmin(posExtCoverRep1, negExtCoverRep1)
quantile(combExtCoverRep1, c(0.5, 0.9, 0.95))


###################################################
### chunk number 17: EatonPeaks
###################################################
peaksRep1 <- slice(combExtCoverRep1, lower = 5)
peakMaxsRep1 <- viewMaxs(peaksRep1)
tail(sort(peakMaxsRep1[[1]]), 30)
peaksRep1 <- peaksRep1[peakMaxsRep1 >= 28]
peakRangesRep1 <-
  GRanges("chrXIV", as(peaksRep1[[1]], "IRanges"),
          seqlengths = seqlengths(rangesRep1))
length(peakRangesRep1)


###################################################
### chunk number 18: EatonExtendedReads
###################################################
data(orcPeaksRep1)
countOverlaps(orcPeaksRep1, peakRangesRep1)
countOverlaps(peakRangesRep1, orcPeaksRep1)


###################################################
### chunk number 19: YeastData
###################################################

###################################################
### chunk number 38: SessionInfo
###################################################
sessionInfo()


