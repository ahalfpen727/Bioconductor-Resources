## ----style, echo = FALSE, results = 'asis'--------------------------------------------------------
options(width=100)
knitr::opts_chunk$set(fig.align = "center")
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(cache=TRUE)

## ----setup, echo=FALSE----------------------------------------------------------------------------
suppressMessages(suppressWarnings(suppressPackageStartupMessages({
library(EnrichmentBrowser)
library(BiocStyle)
library(ALL)
library(hgu95av2.db)
library(airway)
library(regioneR)
library(MultiAssayExperiment)
library(mogsa)
})))

## ----deTbl----------------------------------------------------------------------------------------
deTable <-
     matrix(c(28, 142, 501, 12000),
            nrow = 2,
            dimnames = list(c("DE", "Not.DE"),
                            c("In.gene.set", "Not.in.gene.set")))
deTable

## ----fisher---------------------------------------------------------------------------------------
fisher.test(deTable, alternative = "greater")

## ----loadEBrowser---------------------------------------------------------------------------------
library(EnrichmentBrowser)

## ----loadALL--------------------------------------------------------------------------------------
library(ALL)
data(ALL)

## ----subsetALL------------------------------------------------------------------------------------
ind.bs <- grep("^B", ALL$BT)
ind.mut <- which(ALL$mol.biol %in% c("BCR/ABL", "NEG"))
sset <- intersect(ind.bs, ind.mut)
all.eset <- ALL[, sset]

## ----showALL--------------------------------------------------------------------------------------
dim(all.eset)
exprs(all.eset)[1:4,1:4]

## ----probe2gene-----------------------------------------------------------------------------------
all.eset <- probe.2.gene.eset(all.eset) 
head(featureNames(all.eset))

## ----loadAirway-----------------------------------------------------------------------------------
library(airway)
data(airway)

## ----processAirway--------------------------------------------------------------------------------
air.eset <- as(airway, "ExpressionSet")
annotation(air.eset) <- "hsa"

## ----processAirway2-------------------------------------------------------------------------------
air.eset <- air.eset[grep("^ENSG", rownames(air.eset)), ]
dim(air.eset)
exprs(air.eset)[1:4,1:4]

## ----pdataALL-------------------------------------------------------------------------------------
pData(all.eset)$GROUP <- ifelse(all.eset$mol.biol == "BCR/ABL", 1, 0)
table(pData(all.eset)$GROUP)

## ----pdataAirway----------------------------------------------------------------------------------
pData(air.eset)$GROUP <- ifelse(colData(airway)$dex == "trt", 1, 0)
table(pData(air.eset)$GROUP)

## ----pdataAirway2---------------------------------------------------------------------------------
pData(air.eset)$BLOCK <- colData(airway)$cell
table(pData(air.eset)$BLOCK)

## ----deALL----------------------------------------------------------------------------------------
all.eset <- de.ana(all.eset)
head(fData(all.eset), n=4)

## ----deAirway-------------------------------------------------------------------------------------
air.eset <- de.ana(air.eset, de.method="edgeR")
head(fData(air.eset), n=4)

## ----keggGS, eval=FALSE---------------------------------------------------------------------------
#  kegg.gs <- get.kegg.genesets("hsa")

## ----goGS, eval=FALSE-----------------------------------------------------------------------------
#  go.gs <- get.go.genesets(org="hsa", onto="BP", mode="GO.db")

## ----udefGS---------------------------------------------------------------------------------------
data.dir <- system.file("extdata", package="EnrichmentBrowser")
gmt.file <- file.path(data.dir, "hsa_kegg_gs.gmt")
hsa.gs <- parse.genesets.from.GMT(gmt.file)
length(hsa.gs)
hsa.gs[1:2]

## ----oraALL---------------------------------------------------------------------------------------
ora.all <- sbea(method="ora", eset=all.eset, gs=hsa.gs, perm=0, alpha=0.2)
gs.ranking(ora.all)

## ----browseEA, eval=FALSE-------------------------------------------------------------------------
#  ea.browse(ora.all)

## ----mapIDs---------------------------------------------------------------------------------------
air.eset <- map.ids(air.eset, org="hsa", from="ENSEMBL", to="ENTREZID")
ora.air <- sbea(method="ora", eset=air.eset, gs=hsa.gs, perm=0)
gs.ranking(ora.air)

## ----gseaALL--------------------------------------------------------------------------------------
gsea.all <- sbea(method="gsea", eset=all.eset, gs=hsa.gs, perm=1000)  
gs.ranking(gsea.all)

## ----gseaAir, eval=FALSE--------------------------------------------------------------------------
#  gsea.air <- sbea(method="gsea", eset=air.eset, gs=hsa.gs, perm=100)

## ----roastAir-------------------------------------------------------------------------------------
roast.air <- sbea(method="roast", eset=air.eset, gs=hsa.gs)
gs.ranking(roast.air)  

## ----sbeaM----------------------------------------------------------------------------------------
sbea.methods()

## ----dwnldKegg, eval=FALSE------------------------------------------------------------------------
#  pwys <- download.kegg.pathways("hsa")

## ----compGRN--------------------------------------------------------------------------------------
pwys <- file.path(data.dir, "hsa_kegg_pwys.zip")
hsa.grn <- compile.grn.from.kegg(pwys)
head(hsa.grn)

## ----spiaALL, eval=FALSE--------------------------------------------------------------------------
#  spia.all <- nbea(method="spia", eset=all.eset, gs=hsa.gs, grn=hsa.grn, alpha=0.2)
#  gs.ranking(spia.all)

## ----ggeaALL--------------------------------------------------------------------------------------
ggea.all <- nbea(method="ggea", eset=all.eset, gs=hsa.gs, grn=hsa.grn)
gs.ranking(ggea.all)

## ----nbeaM----------------------------------------------------------------------------------------
nbea.methods()

## ----combEA---------------------------------------------------------------------------------------
res.list <- list(ora.all, gsea.all)
comb.res <- comb.ea.results(res.list)
gs.ranking(comb.res)

## ----loadRegioneR---------------------------------------------------------------------------------
library(regioneR)

## ----loadCPGs-------------------------------------------------------------------------------------
cpgHMM <- toGRanges("http://www.haowulab.org/software/makeCGI/model-based-cpg-islands-hg19.txt")
cpgHMM <- filterChromosomes(cpgHMM, chr.type="canonical")
cpgHMM <- sort(cpgHMM)
cpgHMM

## ----loadProms------------------------------------------------------------------------------------
promoters <- toGRanges("http://gattaca.imppc.org/regioner/data/UCSC.promoters.hg19.bed")
promoters <- filterChromosomes(promoters, chr.type="canonical")
promoters <- sort(promoters)
promoters

## ----subsetCPGProms-------------------------------------------------------------------------------
cpg <- cpgHMM[seqnames(cpgHMM) %in% c("chr21", "chr22")]
prom <- promoters[seqnames(promoters) %in% c("chr21", "chr22")]

## ----permTest-------------------------------------------------------------------------------------
pt <- overlapPermTest(cpg, prom, genome="hg19", ntimes=100, per.chromosome=TRUE, count.once=TRUE)
pt
summary(pt[[1]]$permuted)

## ----loadMAE--------------------------------------------------------------------------------------
library(MultiAssayExperiment)

## ----loadOVMAE------------------------------------------------------------------------------------
data.dir <- system.file("extdata", package="enrichOmics")
mae.file <-  file.path(data.dir, "tcga_ov_mae.rds")
mae <- readr::read_rds(mae.file)
mae

## ----coldataOVMAE---------------------------------------------------------------------------------
head(colData(mae))

## ----matchOVMAE-----------------------------------------------------------------------------------
intersectRows(mae)
intersectColumns(mae)

## ----mogsa----------------------------------------------------------------------------------------
library(mogsa)
assayMats <- assays(mae)
annotSup <- prepSupMoa(as.list(assayMats), geneSets = hsa.gs)
ov.moa <- moa(assayMats, proc.row = "center_ssq1", w.data = "inertia")
smoa <- sup.moa(ov.moa, sup=annotSup, nf=5)
score.moa <- slot(smoa, "score")
rownames(score.moa) <- substring(rownames(score.moa),1,8)
colnames(score.moa) <- substring(colnames(score.moa),1,12)
score.moa[,1:5]

## ----limmaMOGSA-----------------------------------------------------------------------------------
library(limma)
group <- colData(mae)[,"BRCA_Mutation"]
design <- model.matrix(~group)
fit <- lmFit(score.moa, design)
fit <- eBayes(fit)
topTable(fit, number=nrow(score.moa), coef="groupPresent")

## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()

