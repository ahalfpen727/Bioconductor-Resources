## ----style, eval=TRUE, echo=FALSE, results='asis'--------------------------
BiocStyle::latex(use.unsrturl=FALSE)
library(BgeeDB)
## ----echo=TRUE,eval=FALSE--------------------------------------------------
#    library(PROPER)
#    sim.opts.Cheung = RNAseq.SimOptions.2grp(ngenes = 20000, p.DE=0.05,
#      lOD="cheung", lBaselineExpr="cheung")

## ----echo=TRUE,eval=FALSE--------------------------------------------------
#    sim.opts.Bottomly = RNAseq.SimOptions.2grp(ngenes = 20000, p.DE=0.05,
#      lOD="bottomly", lBaselineExpr="bottomly")

## ----echo=TRUE,eval=FALSE,result=FALSE-------------------------------------
#    simres = runSims(Nreps = c(3, 5, 7, 10), sim.opts=sim.opts.Cheung,
#      DEmethod="edgeR", nsims=20)

## ----echo=TRUE,eval=FALSE--------------------------------------------------
#    powers = comparePower(simres, alpha.type="fdr", alpha.nominal=0.1,
#      stratify.by="expr", delta=0.5)

## ----echo=TRUE,eval=FALSE--------------------------------------------------
#    summaryPower(powers)
install.packages("circlize")
## ----eval=FALSE,echo=TRUE--------------------------------------------------
#    plotPower(powers)

## ----eval=FALSE,echo=TRUE--------------------------------------------------
#    plotPowerTD(powers)

## ----eval=FALSE,echo=TRUE--------------------------------------------------
#    plotFDcost(powers)
install.packages('BioCircos')
## ----echo=TRUE,eval=FALSE--------------------------------------------------
#    plotAll(powers)
if (!require('devtools')){install.packages('devtools')}
## ----echo=TRUE,eval=FALSE--------------------------------------------------
#    power.seqDepth(simres, powers)

## ----echo=TRUE,eval=FALSE--------------------------------------------------
#    powers = comparePower(simres, alpha.type="fdr", alpha.nominal=0.1,
#      strata = c(0, 10, 2^(1:7)*10, Inf), filter.by="expr",
#      strata.filtered=1, stratify.by="expr", delta=0.5)

## ----echo=TRUE,eval=FALSE--------------------------------------------------
#      powers = comparePower(simres, alpha.type="fdr", alpha.nominal=0.1,
#        stratify.by="dispersion", target.by="effectsize", delta=1)

## ----echo=TRUE,eval=FALSE--------------------------------------------------
#      powers = comparePower(simres, alpha.type="pval", alpha.nominal=0.001,
#        stratify.by="dispersion", target.by="effectsize", delta=1)

## ----echo=TRUE, result=TRUE------------------------------------------------
sessionInfo()

