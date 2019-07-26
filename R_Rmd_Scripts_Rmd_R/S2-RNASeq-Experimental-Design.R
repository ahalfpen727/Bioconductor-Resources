## ----style, echo = FALSE, results = 'asis'--------------------------------------------------------
BiocStyle::markdown()
options(width=100, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")))

## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE--------------------------------------------
suppressPackageStartupMessages({
    library(DESeq2)
    library(limma)
})

## ----configure-test-------------------------------------------------------------------------------
stopifnot(
    getRversion() >= '3.2' && getRversion() < '3.3',
    BiocInstaller::biocVersion() == "3.2"
)

## ----sleep-t.test---------------------------------------------------------------------------------
head(sleep)
plot(extra ~ group, data = sleep)
## Traditional interface
with(sleep, t.test(extra[group == 1], extra[group == 2]))
## Formula interface
t.test(extra ~ group, sleep)
## equal variance between groups
t.test(extra ~ group, sleep, var.equal=TRUE)

## ----sleep-lm-------------------------------------------------------------------------------------
## linear model; compare to t.test(var.equal=TRUE)
fit <- lm(extra ~ group, sleep)
anova(fit)

## ----sleep-model.matrix---------------------------------------------------------------------------
## underlying model, used in `lm.fit()`
model.matrix(extra ~ group, sleep)     # last column indicates group effect
model.matrix(extra ~ 0 + group, sleep) # contrast between columns

## ----sleep-diff-----------------------------------------------------------------------------------
fit0 <- lm(extra ~ ID, sleep)
fit1 <- lm(extra ~ ID + group, sleep)
anova(fit0, fit1)
t.test(extra ~ group, sleep, var.equal=TRUE, paired=TRUE)

