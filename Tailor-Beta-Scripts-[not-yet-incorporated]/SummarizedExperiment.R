## ----style, echo=FALSE, results='asis'-----------------------------------
BiocStyle::markdown()
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("SummarizedExperiment")
## ----include = FALSE-----------------------------------------------------
# download current version of SE diagram
#download.file("https://docs.google.com/feeds/download/drawings/Export?id=18OcDb80FpvSGRYnFl-8vUqwNNLaNHrG1I9SWKHCselo&exportFormat=svg", "SE.svg")
download.file("https://docs.google.com/feeds/download/drawings/Export?id=1kiC8Qlo1mhSnLDqkGiRNPSo6GWn3C2duBszCFbJCB-g&exportFormat=svg", "SE.svg")
source("https://bioconductor.org/biocLite.R")
biocLite("airway")
library("airway")
## ---- echo=FALSE---------------------------------------------------------
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(data(airway, package="airway"))

## ------------------------------------------------------------------------
library(SummarizedExperiment)
data(airway, package="airway")
se <- airway
se
assays(se)$counts

knitr::kable(assays(se)$counts[1:10,])

## ----rowRanges-----------------------------------------------------------
rowRanges(se)
colData(se)
se[, se$dex == "trt"]
metadata(se)
metadata(se)$formula <- counts ~ dex + albut
metadata(se)
nrows <- 200
ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, TRUE),
                     feature_id=sprintf("ID%03d", 1:200))
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])

SummarizedExperiment(assays=list(counts=counts),
                     rowRanges=rowRanges, colData=colData)
SummarizedExperiment(assays=list(counts=counts), colData=colData)

## ----2d------------------------------------------------------------------
# subset the first five transcripts and first three samples
se[1:5, 1:3]
se[, se$cell == "N61311"]
counts <- matrix(1:15, 5, 3, dimnames=list(LETTERS[1:5], LETTERS[1:3]))
dates <- SummarizedExperiment(assays=list(counts=counts),
                              rowData=DataFrame(month=month.name[1:5], day=1:5))
dates[rowData(dates)$month == "January", ]
assays(se)
assays(se)[[1]][1:5, 1:5]
assay(se)[1:5, 1:5]
assay(se, 1)[1:5, 1:5]

roi <- GRanges(seqnames="1", ranges=100000:1100000)
subsetByOverlaps(se, roi)

## ----rseSubclass---------------------------------------------------------
setClass("MyRSESubclass",
    contains="RangedSummarizedExperiment",
    representation=representation(
        slot1="integer",
        slot2="function"
        ## ... maybe more ...))
library(survival)
###################################################
### code chunk number 4: surv
###################################################
data(leukemia)
head(leukemia)
help(Surv)
mysurv <- Surv(leukemia$time, leukemia$status)
head(mysurv)


###################################################
### code chunk number 5: survival
###################################################
leuk.km <- survfit(Surv(time, status) ~ x, data=leukemia)
plot(leuk.km, lty=1, col=c("darkblue","darkred"))
legend("topright", legend=c('Maintain', 'Non-main'), lty=1:2, col=c("darkblue","darkred"))
leuk.km2 <- survfit(Surv(time, status) ~ x, data=leukemia, conf.type='log-log')
summary(leuk.km2)
plot(leuk.km2, mark.time=FALSE, conf.int=TRUE, lty=1, col=c("darkblue","darkred"))
legend("topright", legend=c('Maintain', 'Non-main'), lty=1, 
	   col=c("darkblue","darkred"))
survdiff(Surv(time, status) ~ x,
                   data=leukemia)


###################################################
### code chunk number 8: survivalkm (eval = FALSE)
###################################################
colon.km <- survfit(Surv(time, status) ~ rx, data=colon)
plot(colon.km, lty=1, col=c("darkblue", "darkgreen", "darkred"))
legend("topright", legend=c("Obs", "Lev", "Lev+5FU"), lty=1, col=c("darkblue", "darkgreen", "darkred"))
## logrank test
survdiff(formula = Surv(time, status) ~ rx, data = colon)
leuk.ph <- coxph(Surv(time, status) ~ x, data=leukemia)
summary(leuk.ph)
summary(coxph(Surv(time, status) ~ ., data=colon))


source("http://www.bioconductor.org/biocLite.R")
biocLite("survcomp")
library(survcomp)
data(breastCancerData)


###################################################
### code chunk number 14: bcataset
###################################################
print(vdx7g)
print(transbig7g)


###################################################
### code chunk number 15: vdx7gsurvdata
###################################################
## clincal information
print(head(phenoData(vdx7g)@data))
dd <- phenoData(vdx7g)@data[ ,c("t.dmfs", "e.dmfs")]
colnames(dd) <- c("time", "event")
## append gene expressions
ge <- t(exprs(vdx7g))
dd <- cbind(dd, ge)


###################################################
### code chunk number 16: vdxfitcox
###################################################
mm <- coxph(Surv(time, event) ~ ., data=data.frame(dd))
print(summary(mm))
predtr <- predict(mm, newdata=data.frame(dd), type="risk")


###################################################
### code chunk number 17: transbigvalidatecox
###################################################
dd2 <- t(exprs(transbig7g))
predts <- predict(mm, newdata=data.frame(dd2), type="risk")


###################################################
### code chunk number 18: transbigsurvd
###################################################
survdd2 <- phenoData(transbig7g)@data[ ,c("t.dmfs", "e.dmfs")]
colnames(survdd2) <- c("time", "event")


###################################################
### code chunk number 19: hr
###################################################
perf <- hazard.ratio(x=predts, surv.time=survdd2[ ,"time"], surv.event=survdd2[ ,"event"], na.rm=TRUE)
print(perf[1:6])


###################################################
### code chunk number 20: dindex
###################################################
perf <- D.index(x=predts, surv.time=survdd2[ ,"time"], surv.event=survdd2[ ,"event"], na.rm=TRUE)
print(perf[1:6])


###################################################
### code chunk number 21: cindex
###################################################
perf <- concordance.index(x=predts, surv.time=survdd2[ ,"time"], surv.event=survdd2[ ,"event"], method="noether", na.rm=TRUE)
print(perf[1:5])


###################################################
### code chunk number 22: tdrocc
###################################################
perf <- tdrocc(x=predts, surv.time=survdd2[ ,"time"], surv.event=survdd2[ ,"event"], time=5*365, na.rm=TRUE)
plot(x=1-perf$spec, y=perf$sens, type="l", xlab="1 - specificity", ylab="sensitivity", xlim=c(0, 1), ylim=c(0, 1), main="Time-dependent ROC curve\nat 5 years")
lines(x=c(0,1), y=c(0,1), lty=3, col="red")


###################################################
### code chunk number 23: sbrier
###################################################
ddtr <- cbind("time"=dd[ ,"time"], "event"=dd[ ,"event"], "score"=predtr)
ddts <- cbind("time"=survdd2[ ,"time"], "event"=survdd2[ ,"event"], "score"=predts)
perf <- sbrier.score2proba(data.tr=data.frame(ddtr), data.ts=data.frame(ddts), method="cox")
plot(x=perf$time, y=perf$bsc, xlab="Time (days)", ylab="Brier score", type="l")
## null model
ddtr <- cbind("time"=dd[ ,"time"], "event"=dd[ ,"event"], "score"=1)
ddts <- cbind("time"=survdd2[ ,"time"], "event"=survdd2[ ,"event"], "score"=1)
perfnull <- sbrier.score2proba(data.tr=data.frame(ddtr), data.ts=data.frame(ddts), method="cox")
lines(x=perfnull$time, y=perfnull$bsc, col="red", lty=2)
## legend
legend("bottomright", title="Integrated Brier score", legend=c(sprintf("Cox model = %.3g", perf$bsc.integrated), sprintf("Null model = %.3g", perfnull$bsc.integrated)), col=c("black", "red"), lty=c(1,2))

stamFile <-system.file("data", "stam.Rda", package="SeattleIntro2012Data")
load(stamFile)
68> 
ridx <- rowSums(assays(stam)[["Tags"]] > 0) == ncol(stam)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
tx <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
tss <- resize(tx, width=1)
peak<-resize(rowData(stam)[ridx], width=1, fix="center")


idx <- nearest(peak, tss)
sgn <- as.numeric(ifelse(strand(tss)[idx] == "+", 1, -1))
> dist <- (start(peak) - start(tss)[idx]) * sgn
Here we summarize the distances as a simple table and density plot, focusing
on binding sites within 1000 bases of a transcription start site; the density plot
is in Figure 7.
>
>
>
>
bound <- 1000
ok <- abs(dist) < bound
dist <- dist[ok]
table(sign(dist))
-1
1262
0
4
1
707
> print(densityplot(dist[ok], plot.points=FALSE,
+
xlab="Distance to Nearest TSS"))
The distance to transcript start site is a useful set of operations, so let’s make
it a re-usable function
> distToTss <-
+
function(peak, tx)
+ {
+
tss <- resize(tx, width=1)
+
idx <- nearest(peak, tss)
+
sgn <- as.numeric(ifelse(strand(tss)[idx] == "+", 1, -1))
+
(start(peak) - start(tss)[idx]) * sgn
+ }

