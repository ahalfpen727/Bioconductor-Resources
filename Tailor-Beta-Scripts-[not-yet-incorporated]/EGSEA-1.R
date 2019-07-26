## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"--------------------
BiocStyle::latex()

## ----echo=TRUE,  eval=FALSE,tidy=TRUE,tidy.opts=list(blank=FALSE, width.cutoff=60)----
#  source("http://www.bioconductor.org/biocLite.R")
#  biocLite(c("PADOG", "GSVA", "AnnotationDbi", "topGO", "pathview", "gage",
#  "globaltest", "limma", "edgeR", "safe", "org.Hs.eg.db", "org.Mm.eg.db",
#  "org.Rn.eg.db"))

## ----echo=TRUE,  eval=FALSE,tidy=TRUE--------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("EGSEAdata")

## ----echo=TRUE,  eval=FALSE,tidy=TRUE--------------------------------------
#  library(devtools)
#  install_bitbucket("malhamdoosh/egseadata", ref="Stable_Release")

## ----echo=TRUE,  eval=FALSE,tidy=TRUE--------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("EGSEA")

## ----echo=TRUE,  eval=FALSE,tidy=TRUE--------------------------------------
#  library(BiocInstaller)
#  useDevel()
#  biocLite("EGSEA")

## ----echo=TRUE,  eval=FALSE,tidy=TRUE--------------------------------------
#  library(devtools)
#  install_bitbucket("malhamdoosh/egsea", ref="Devel_Release")

## ----echo=TRUE,  eval=TRUE,tidy=TRUE, message=F, warning=F-----------------
library(EGSEA)

## ----echo=TRUE,  eval=TRUE,tidy=TRUE---------------------------------------
library(EGSEAdata)
data(il13.data)
v = il13.data$voom
names(v)
v$design
contrasts = il13.data$contra
contrasts

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE, width.cutoff=60)----
gs.annots = buildIdx(entrezIDs=rownames(v$E), species="human", msigdb.gsets="c5",
 			kegg.exclude = c("Metabolism"))
names(gs.annots)

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE, width.cutoff=60)----
summary(gs.annots$kegg)
summary(gs.annots$c5)

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE, width.cutoff=60)----
baseMethods = egsea.base()[-c(2, 12)]
baseMethods

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE, width.cutoff=60)----
egsea.sort()

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE, width.cutoff=60)----
# perform the EGSEA analysis
# set report = TRUE to generate HTML report.
# set display.top = 20 to display more gene sets. It takes longer time to run.
gsa = egsea(voom.results=v, contrasts=contrasts,  gs.annots=gs.annots,
 			symbolsMap=v$genes, baseGSEAs=baseMethods,
            report.dir="./il13-egsea-report",
 			 sort.by="avg.rank", num.threads = 4, report=FALSE)

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE, width.cutoff=60)----
summary(gsa)

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE, width.cutoff=60)----
gs.annots = buildIdx(entrezIDs=rownames(v$E), species="human", gsdb.gsets="all")
names(gs.annots)

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE,width.cutoff=60)----
show(gsa)

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE,width.cutoff=60)----
topSets(gsa, contrast=1, gs.label="kegg", number = 10)

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE,width.cutoff=60)----
t = topSets(gsa, contrast=1, gs.label="c5", sort.by="ora", number = 10, names.only=FALSE)
t

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE,width.cutoff=60)----
showSetByName(gsa, "c5", rownames(t)[1])

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE,width.cutoff=60)----
t = topSets(gsa, contrast="comparison", gs.label="kegg", number = 10)
t

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE,width.cutoff=60)----

showSetByName(gsa, "kegg",  rownames(t)[1])


## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE,width.cutoff=60)----
plotMethods(gsa, gs.label = "kegg", contrast=1, file.name="X24IL13-X24-kegg-methods")

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE,width.cutoff=60)----
plotSummary(gsa, gs.label="kegg", contrast = 1, file.name="X24IL13-X24-kegg-summary")

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE,width.cutoff=60)----
showSetByID(gsa, gs.label="kegg", c("hsa04060", "hsa04640"))

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE,width.cutoff=60)----
plotGOGraph(gsa, gs.label="c5", file.name="X24IL13-X24-c5-top-", sort.by="avg.rank")

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE,width.cutoff=60)----
plotHeatmap(gsa, "Asthma", gs.label="kegg", contrast=1,
        file.name="asthma-hm")

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE,width.cutoff=60)----
plotPathway(gsa, "Asthma", gs.label="kegg", file.name="asthma-pathway")

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE,width.cutoff=60)----
plotSummary(gsa, gs.label="kegg", contrast=c(1,2), file.name="kegg-summary-cmp")

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE,width.cutoff=60)----
plotSummaryHeatmap(gsa, gs.label="kegg", show.vals = "p.adj",
        file.name="il13-sum-heatmap")

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE,width.cutoff=60)----
plotHeatmap(gsa, "Asthma", gs.label="kegg", contrast="comparison",
        file.name="asthma-hm-cmp")

## ----echo=TRUE, eval=TRUE,tidy=TRUE,tidy.opts=list(blank=FALSE,width.cutoff=60)----
plotPathway(gsa, "Asthma", gs.label="kegg", contrast=0,  file.name="asthma-pathway-cmp")

## ----echo=TRUE,  eval=FALSE,tidy=TRUE,tidy.opts=list(blank=FALSE, width.cutoff=60)----
#  # load the mammary dataset
#  library(EGSEA)
#  library(EGSEAdata)
#  data(mam.data)
#  v = mam.data$voom
#  names(v)
#  v$design
#  contrasts = mam.data$contra
#  contrasts
#  # build the gene set collections
#  gs.annots = buildIdx(entrezIDs=rownames(v$E), species="mouse",
#          msigdb.gsets = "c2",
#          kegg.exclude = "all")
#  names(gs.annots)
#  # create Entrez IDs - Symbols map
#  symbolsMap = v$genes[,c(1,3)]
#  colnames(symbolsMap) = c("FeatureID", "Symbols")
#  symbolsMap[, "Symbols"] = as.character(symbolsMap[, "Symbols"])
#  # replace NA Symbols with IDs
#  na.sym = is.na(symbolsMap[, "Symbols"])
#  na.sym
#  symbolsMap[na.sym, "Symbols"] = symbolsMap[na.sym, "FeatureID"]
#  # perform the EGSEA analysis
#  # set report = TRUE to generate the EGSEA interactive report
#  baseMethods = c("camera", "safe", "gage", "padog", "zscore",
#          "gsva", "globaltest", "ora")
#  gsa = egsea(voom.results=v, contrasts=contrasts, gs.annots=gs.annots,
#          symbolsMap=symbolsMap, baseGSEAs=baseMethods,
#          sort.by="med.rank",
#          num.threads=4, report=FALSE)
#  # show top 20 comparative gene sets in C2 collection
#  summary(gsa)
#  topSets(gsa, gs.label="c2", contrast="comparison", number = 20)

## ----echo=TRUE,  eval=FALSE,tidy=TRUE,tidy.opts=list(blank=FALSE, width.cutoff=60)----
#  # load the count matrix and other relevant data
#  library(EGSEAdata)
#  data(il13.data.cnt)
#  cnt = il13.data.cnt$counts
#  group = il13.data.cnt$group
#  group
#  design = il13.data.cnt$design
#  contrasts = il13.data.cnt$contra
#  genes = il13.data.cnt$genes
#  # build the gene set collections
#  gs.annots = buildIdx(entrezIDs=rownames(cnt), species="human", msigdb.gsets="none",
#   			 kegg.exclude = c("Metabolism"))
#  # perform the EGSEA analysis
#  # set report = TRUE to generate the EGSEA interactive report
#  gsa = egsea.cnt(counts=cnt, group=group, design=design, contrasts=contrasts,
#          gs.annots=gs.annots, symbolsMap=genes, baseGSEAs=egsea.base()[-c(2,12)],
#   	    sort.by="avg.rank", num.threads = 4, report=FALSE)

## ----echo=TRUE, eval=TRUE, tidy=TRUE, tidy.opts=list(blank=FALSE, width.cutoff=60)----
# load IL-13 dataset
library(EGSEAdata)
data(il13.data)
voom.results = il13.data$voom
contrast = il13.data$contra
# find Differentially Expressed genes
library(limma)
vfit = lmFit(voom.results, voom.results$design)
vfit = contrasts.fit(vfit, contrast)
vfit = eBayes(vfit)
# select DE genes (Entrez IDs and logFC) at p-value <= 0.05 and |logFC| >= 1
top.Table = topTable(vfit, coef=1, number=Inf, p.value=0.05, lfc=1)
deGenes = as.character(top.Table$FeatureID)
logFC = top.Table$logFC
names(logFC) = deGenes
# build the gene set collection index
gs.annots = buildIdx(entrezIDs=deGenes, species="human", msigdb.gsets="none",
 			kegg.exclude = c("Metabolism"))
# perform the ORA analysis
# set report = TRUE to generate the EGSEA interactive report
gsa = egsea.ora(geneIDs=deGenes, universe= as.character(voom.results$genes[,1]),
				logFC =logFC, title="X24IL13-X24",  gs.annots=gs.annots,
  			    symbolsMap=top.Table[, c(1,2)], display.top = 5,
  			    report.dir="./il13-egsea-ora-report", num.threads = 4, report=FALSE)


## ----echo=TRUE, eval=TRUE, tidy=TRUE, tidy.opts=list(blank=FALSE, width.cutoff=60)----
library(EGSEAdata)
data(il13.data)
v = il13.data$voom
# load KEGG pathways
data(kegg.pathways)
# select 50 pathways
gsets = kegg.pathways$human$kg.sets[1:50]
gsets[1]
# build custom gene set collection using these 50 pathways
gs.annot = buildCustomIdx(geneIDs=rownames(v$E), gsets= gsets, species="human")
class(gs.annot)
show(gs.annot)

## ----echo=TRUE, eval=TRUE, tidy=TRUE, tidy.opts=list(blank=FALSE, width.cutoff=60)----
sessionInfo()

