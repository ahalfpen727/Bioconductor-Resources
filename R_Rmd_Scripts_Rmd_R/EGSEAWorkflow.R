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

## ----setup, message=FALSE, echo = FALSE------------------------------------
library(limma)
library(edgeR)
url = "http://bioinf.wehi.edu.au/EGSEA/mam.rnaseq.rdata"
utils::download.file(url, destfile="mam.rnaseq.rdata", mode="wb") 
load("mam.rnaseq.rdata")
names(mam.rnaseq.data)
dim(mam.rnaseq.data)
x = calcNormFactors(mam.rnaseq.data, method = "TMM")
design = model.matrix(~0+x$samples$group+x$samples$lane)
colnames(design) = gsub("x\\$samples\\$group", "", colnames(design))
colnames(design) = gsub("x\\$samples\\$lane", "", colnames(design))
head(design)
contr.matrix = makeContrasts(
         BasalvsLP = Basal-LP,
         BasalvsML = Basal - ML,
         LPvsML = LP - ML,
         levels = colnames(design))
head(contr.matrix)

## ----setup2, eval=TRUE-----------------------------------------------------
v = voom(x, design, plot=FALSE)
names(v)

## ----collections-----------------------------------------------------------
library(EGSEAdata)
egsea.data("mouse")

## ----collectionlist--------------------------------------------------------
info = egsea.data("mouse", returnInfo = TRUE)
names(info)
info$msigdb$info$collections

## ----loadegsea, message=FALSE, warning=FALSE-------------------------------
library(EGSEA)

## ----indexing--------------------------------------------------------------
gs.annots = buildIdx(entrezIDs=v$genes$ENTREZID, species="mouse", 
           msigdb.gsets=c("c2", "c5"), go.part = TRUE)
names(gs.annots)

## ----exploresets-----------------------------------------------------------
class(gs.annots$c2)
summary(gs.annots$c2)
show(gs.annots$c2)
s = getSetByName(gs.annots$c2, "SMID_BREAST_CANCER_LUMINAL_A_DN")
class(s)
names(s)
names(s$SMID_BREAST_CANCER_LUMINAL_A_DN)

## ----indexclass------------------------------------------------------------
slotNames(gs.annots$c2)

## ----symbolmap-------------------------------------------------------------
colnames(v$genes)
symbolsMap = v$genes[, c(1, 2)]
colnames(symbolsMap) = c("FeatureID", "Symbols")
symbolsMap[, "Symbols"] = as.character(symbolsMap[, "Symbols"])

## ----base------------------------------------------------------------------
egsea.base()

## ----selectbasemethods-----------------------------------------------------
baseMethods = egsea.base()[-2]
baseMethods

## ----combine---------------------------------------------------------------
egsea.combine()

## ----sort------------------------------------------------------------------
egsea.sort()

## ----egseatest-------------------------------------------------------------
gsa = egsea(voom.results=v, contrasts=contr.matrix,  
         gs.annots=gs.annots, symbolsMap=symbolsMap,
         baseGSEAs=baseMethods, sort.by="med.rank",
         num.threads = 16, report = FALSE)

## ----showegsea-------------------------------------------------------------
show(gsa)

## ----summariseegsea--------------------------------------------------------
summary(gsa)

## ----topsets---------------------------------------------------------------
topSets(gsa, gs.label="c2", contrast = "comparison", names.only=TRUE)

## ----topsetslim------------------------------------------------------------
t = topSets(gsa, contrast = "comparison",
             names.only=FALSE, number = Inf, verbose = FALSE)
t[grep("LIM_", rownames(t)), c("p.adj", "Rank", "med.rank", "vote.rank")]

## ----topsets2--------------------------------------------------------------
topSets(gsa, gs.label="kegg", contrast="BasalvsLP", sort.by="med.rank")
topSets(gsa, gs.label="kegg", contrast="comparison", sort.by="med.rank")

## ----heatmaps--------------------------------------------------------------
plotHeatmap(gsa, gene.set="LIM_MAMMARY_STEM_CELL_UP", gs.label="c2",
         contrast = "comparison", file.name = "hm_cmp_LIM_MAMMARY_STEM_CELL_UP", format="png")
plotHeatmap(gsa, gene.set="LIM_MAMMARY_STEM_CELL_DN", gs.label="c2",
         contrast = "comparison", file.name = "hm_cmp_LIM_MAMMARY_STEM_CELL_DN", format="png")

## ----pathwayplot1, eval=FALSE----------------------------------------------
#  plotPathway(gsa, gene.set = "Vascular smooth muscle contraction",
#               contrast = "BasalvsLP", gs.label = "kegg",
#               file.name = "Vascular_smooth_muscle_contraction")

## ----pathwayplot2, eval=FALSE----------------------------------------------
#  plotPathway(gsa, gene.set = "Vascular smooth muscle contraction",
#               contrast = "comparison", gs.label = "kegg",
#               file.name = "Vascular_smooth_muscle_contraction_cmp")

## ----mdsplot---------------------------------------------------------------
plotMethods(gsa, gs.label = "c2", contrast = "BasalvsLP", 
         file.name = "mds_c2_BasalvsLP", format="png")
plotMethods(gsa, gs.label = "c5BP", contrast = "BasalvsLP", 
         file.name = "mds_c5_BasalvsLP", format="png")

## ----keggsummaryplot1------------------------------------------------------
plotSummary(gsa, gs.label = 3, contrast = 3, 
         file.name = "summary_kegg_LPvsML", format="png")

## ----c2summaryplot2--------------------------------------------------------
plotSummary(gsa, gs.label = 1, contrast = 3, 
         file.name = "summary_c2_LPvsML", 
         x.axis = "med.rank", format="png")

## ----c2summaryplot3--------------------------------------------------------
plotSummary(gsa, gs.label = 1, contrast = 3, 
         file.name = "summary_sig_c2_LPvsML", 
         x.axis = "med.rank", x.cutoff=300, format="png")

## ----summaryplotkegg1and2--------------------------------------------------
plotSummary(gsa, gs.label = "kegg", contrast = c(1,2), 
         file.name = "summary_kegg_1vs2", format="png")

## ----gographs--------------------------------------------------------------
plotGOGraph(gsa, gs.label="c5BP", contrast = 1, file.name="BasalvsLP-c5BP-top-", format="png")
plotGOGraph(gsa, gs.label="c5CC", contrast = 1, file.name="BasalvsLP-c5CC-top-", format="png")

## ----summarybarplot--------------------------------------------------------
plotBars(gsa, gs.label = "c2", contrast="comparison", file.name="comparison-c2-bars", format="png")

## ----summaryheatmap--------------------------------------------------------
plotSummaryHeatmap(gsa, gs.label="c2", hm.vals = "avg.logfc.dir",
         file.name="summary_heatmaps_c2", format="png")
plotSummaryHeatmap(gsa, gs.label="kegg", hm.vals = "avg.logfc.dir",
         file.name="summary_heatmaps_kegg", format="png")

## ----toptable--------------------------------------------------------------
t = limmaTopTable(gsa, contrast=1)
head(t)

## ----htmlreport, warning=FALSE, eval=FALSE---------------------------------
#  generateReport(gsa, number = 20, report.dir="./mam-rnaseq-egsea-report")

## ----mareadidats, eval=TRUE------------------------------------------------
library(limma)
url = "http://bioinf.wehi.edu.au/EGSEA/arraydata.zip"
utils::download.file(url, destfile="arraydata.zip", mode="wb")
utils::unzip("arraydata.zip", exdir = ".")
targets = read.delim("targets.txt", header=TRUE, sep=" ")
data = read.idat(as.character(targets$File), 
                   bgxfile="GPL6887_MouseWG-6_V2_0_R0_11278593_A.bgx",
                   annotation=c("Entrez_Gene_ID","Symbol", "Chromosome"))
data$other$Detection = detectionPValues(data)
data$targets = targets
colnames(data) = targets$Sample

## ----manormalize, eval=TRUE------------------------------------------------
data = neqc(data)

## ----mafilter, eval=TRUE---------------------------------------------------
table(targets$Celltype)
keep.exprs = rowSums(data$other$Detection<0.05)>=5
table(keep.exprs)
data = data[keep.exprs,]
dim(data)
head(data$genes)
sum(is.na(data$genes$Entrez_Gene_ID))
data1 = data[!is.na(data$genes$Entrez_Gene_ID), ]
dim(data1)
ord = order(lmFit(data1)$Amean, decreasing=TRUE)
ids2keep = data1$genes$Array_Address_Id[ord][!duplicated(data1$genes$Entrez_Gene_ID[ord])]
data1 = data1[match(ids2keep, data1$genes$Array_Address_Id),]
dim(data1)

expr = data1$E
group = as.factor(data1$targets$Celltype)
probe.annot = data1$genes[, 2:4]
head(probe.annot)

## ----malinearmodel, eval=TRUE----------------------------------------------
head(data1$targets)
experiment = as.character(data1$targets$Experiment)
design = model.matrix(~0 + group + experiment)
colnames(design) = gsub("group", "", colnames(design))
design
contr.matrix = makeContrasts(
         BasalvsLP = Basal-LP,
         BasalvsML = Basal-ML,
         LPvsML = LP-ML,
         levels = colnames(design))
contr.matrix

## ----maindex, eval=FALSE---------------------------------------------------
#  library(EGSEA)
#  library(EGSEAdata)
#  gs.annots = buildIdx(entrezIDs=probe.annot[, 2],
#               species="mouse",
#               msigdb.gsets=c("c2", "c5"), go.part = TRUE)
#  names(gs.annots)

## ----maegsea, eval=FALSE---------------------------------------------------
#  baseMethods = egsea.base()[-2]
#  baseMethods
#  
#  gsam = egsea.ma(expr=expr, group=group,
#   		  	 probe.annot = probe.annot,
#   		  	 design = design,
#           contrasts=contr.matrix,
#           gs.annots=gs.annots,
#           baseGSEAs=baseMethods, sort.by="med.rank",
#           num.threads = 8, report = FALSE)

## ----mareport, warning=FALSE, eval=FALSE-----------------------------------
#  generateReport(gsam, number = 20, report.dir="./mam-ma-egsea-report")

## ----matopsets, eval=FALSE-------------------------------------------------
#  topSets(gsam, gs.label="c2", contrast="comparison", names.only=TRUE, number=5)

## ----softwareinfo----------------------------------------------------------
sessionInfo()

