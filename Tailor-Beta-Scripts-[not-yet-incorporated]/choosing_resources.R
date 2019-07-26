## ----first_load, echo=TRUE, warning=FALSE, results='hide', message=FALSE----
library(RITANdata)
library(RITAN)
library(knitr)
kable( attr(network_list, 'network_data_sources') )
kable( sapply(geneset_list, length), col.names = c('# Genesets') )

## ----geneset_overlap_P, echo=TRUE, eval=FALSE, warning=FALSE, fig.width = 7, fig.height = 7, fig.align='center'----
#  genes <- geneset_list$MSigDB_C5$APICAL_JUNCTION_COMPLEX
#  e <- network_overlap( genes, resources = c('HPRD','CCSB','Biogrid','STRING'), minStringScore = 700 )
#  ## We also strongly encourage use of the BioPlex database, which we do not distribute with RITAN in compliance with their licensing.

## ----geneset_overlap_M, echo=TRUE, eval=FALSE, warning=FALSE, fig.width = 7, fig.height = 7, fig.align='center'----
#  genes2 <- geneset_list$MSigDB_C5$AMINE_METABOLIC_PROCESS
#  e2 <- network_overlap( genes, resources = c('PID','HumanNet','Biogrid') )

## ----geneset_overlap_E, echo=TRUE, eval=FALSE, warning=FALSE, fig.width = 7, fig.height = 7, fig.align='center'----
#  genes2 <- geneset_list$MSigDB_C5$AMINE_METABOLIC_PROCESS
#  e2 <- network_overlap( genes, resources = c('ChEA','HumanNet') )

