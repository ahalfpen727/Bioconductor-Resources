#############################################################################
#        R CODE FOR ILLUMINA BEADCHIPS ANALYSIS USING LUMI AND LIMMA        #
#############################################################################

## YOU NEED 3 FILES TO RUN THIS SCRIPT ###

# 1. Metadata consisting of sample names and array barcodes  (see instructions on how to create this document in help on format).

# 2. Annotation file (see example in the help on format document)

# 3. Illumina microarray raw data (see example in the help on format document)

# see tips on how to run this code in the help on format document

#########################################################################
######### PART OF THE CODE (VARIABLES) THAT YOU NEED TO MODIFY ##########
#########################################################################

### NOTE: if you followed exactly the format of the 3 input files, you don't need to modify anything else in the code, just run the code line by line. 

#The following variables need to be assigned inside the code, e.g.,
#- workDirectory <- "C:/Users/changjiang/Documents/ut/LA02"
#- metadataFileName <- "SampleBarcode_minus_65055.txt"     
#- IlluminaDataFileName <- "sample_probe_nonnorm_FinalReport.txt"
#- annotationFileName <- "HumanHT-12_V4_0_R2_15002873_B.txt"

workDirectory <- "/Users/veroniquevoisin/Documents/Laurie Ailles/LA03/DEBUG_CODE"
metadataFileName <- "SampleBarcode.txt"
annotationFileName <- "HumanHT-12_V4_0_R2_15002873_B.txt"
IlluminaDataFileName <- "sample_probe_nonnorm_FinalReport.txt"
group1 <- "T"  # e.g "T"
group2 <- "E" # e.g "C"
designType = "paired"  #choice between "unpaired" or "paired"

### NOT IMPLEMENTED YET: do not change
chiptype <- 'illumina_humanht_12_v4' ##### list number of choice here
species <- 'human'  ### choices are 'human' or 'mouse'

### it is possible to add 1 extra confounding factor e.g tumor location, technical batch, tumor grade, male or female: the additional factor is an optional column in metadata called "variableX" and is incorporated in the model fitting. 

 ########################################################################
###############     R PACKAGES AND WORK DIRECTORY       ################
########################################################################

###  note : you just need to install the packages (also called library) once but YOU NEED to load the libraries each time you opened your workspace

sessionInfo() # give you information about package versions and R versions that you used, information that you may need for publication purpose.

#(1) install packages
#(2) load packages
#(3) set work directory

##Install packages
#source("http://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
#biocLite()   #This installs a number of base packages
#biocLite("limma")
#biocLite("lumi")
#####biocLite("affy")
######biocLite("preprocessCore") #In BioC software repository
#biocLite("gplots")
#biocLite("scatterplot3d")
####biocLite("lumiHumanAll.db")
#biocLite("ggplot2")
#biocLite("biomaRt")

#update.packages(repos = biocinstallRepos())
#install.packages("mgcv")

##Load library
library(lumi)
library(limma)
library(scatterplot3d)
#library(affy)
#library(preprocessCore)
library(gplots)
library(gplots2)


##work directory

###
setwd(workDirectory)
#date <- date()
#date <- gsub(" ", "", date)
#temp = paste0(group1, "vs", group2, "_", date)
subDir <- "temp"
dir.create(file.path(workDirectory, subDir), showWarnings = FALSE)
getwd()
dir()

########################################################################
###############    METADATA: SAMPLE AND ARRAY BARCODES     #############
########################################################################


##Read metadata file

sample_barcode <- read.table(metadataFileName, header=TRUE, stringsAsFactors=FALSE)
dim(sample_barcode) #36  3
sample_barcode
head(sample_barcode)
tail(sample_barcode)

##unique sample ID: sample_barcode$SampleNumber
if (nrow(sample_barcode) != length(unique(sample_barcode$SampleNumber)))#36
 { print( "warning: your metadata contains non unique sample name")}
 
 
##add patient, marker, and chip
patient <- unlist(lapply(strsplit(sample_barcode$sampleName, split="_"), function(s) s[1]))
sample_barcode$patient <- patient
table(sample_barcode$patient)

#populations (markers)
#marker <- unlist(lapply(strsplit(sample_barcode$sampleName, split="_"), function(s) s[2]))
marker <- gsub(".*_", "", sample_barcode$sampleName)
sample_barcode$marker <- marker
table(sample_barcode$marker)

#chips
chip <- unlist(lapply(strsplit(sample_barcode$ArrayBarcode, split="_"), function(s) s[1]))
chip <- factor(chip)
table(chip)
chip
table(chip)
sample_barcode$chip <- chip

sample_barcode 

########################################################################
###############     PRE_PROCESSING USING LUMI PACKAGE       ############
########################################################################

library(lumi)
#(1) Read raw intensity data
#(2) Normalization (In practice, this step should be after quality control.)
#(3) Quality control (before and after normalization)

####################
#(1) Read raw data
####################

##Illumina data file name
datalumi <- lumiR(IlluminaDataFileName) # read the data
pData(datalumi) # check if the samples have been loaded corectly

##change column names using sampleName in sample_barcode
if (all(colnames(datalumi) %in% sample_barcode$ArrayBarcode) ) {
	indx <- match(colnames(datalumi), sample_barcode$ArrayBarcode)
    newSampleNames <- sample_barcode$sampleName
    colnames(datalumi)[!is.na(indx)] <- newSampleNames[indx[!is.na(indx)]]
    colnames(datalumi)
} else {
  	print("ERROR at change column names, check array barcode")
}   

##summary of the raw data
pData(datalumi)         # check the output to see whether the column names have been changed
datalumi@QC             # display the quality control information of the LumiBatch object

####################
#(2) Normalization
####################

##Variance Stabilization & Quantile normalization
##From raw Illumina probe intensities to expression values 
##return a processed LumiBatch object

datalumiQN <- lumiExpresso(datalumi, QC.evaluation=TRUE, varianceStabilize.param=list(method='log2'), normalize.param = list(method='quantile'))
summary(datalumiQN, 'QC')
pData(datalumiQN)


###################################################################################
###############    QUALITY CONTROL PLOTS  ###########################
#####################################################################################


##boxplot and density plot of both raw and normalized intensities on log2 scale
pdf("temp/density_boxplot.pdf", width=11, height=11) # to save as the plots as pdf
#png("temp/density_boxplot.png", width = 1.5*480, height = 1.5*480) # to save the plots as png
par(mfrow=c(2,2), cex=0.8) # to create 4 panels
plot(datalumi, what='density', main="Raw intensity on log2 scale", addLegend=F)  
plot(datalumi, what='boxplot', main="Raw intensity on log2 scale")
plot(datalumiQN, what='density', main="Quantile normalized data", addLegend=F)
plot(datalumiQN, what='boxplot', main="Quantile normalized data") 
dev.off() # to save an image on the local computer

##Pairs plot of sample intensities in an ExpressionSet object 
#pairs(log2(ave.sig[sample(1:nrow(ave.sig), 5000), c(1,10)] + 1), pch=16, cex=0.5, lower.panel = NULL)
#pairs(datalumi[, indx]) 
indx <- c(1, 3) #a pair of two samples
plot(datalumi[, indx], what='pair') # before normalization
plot(datalumiQN[, indx], what='pair') # after normalization

##MA Plots
indx <- c(1:2, 4:5) #any 4 samples
MAplot(datalumi[, indx], subset=NULL) # before normalization
MAplot(datalumiQN[, indx], subset=NULL) # after normalization


##############################################
#############  CLUSTERING  ###################
##############################################

##1. Based on large coefficient of variance (mean / standard variance)
##Estimate the sample relations based on selected probes (the probes are selected if their coefficient of variance (mean / standard variance) are greater a threshold). 
##Two methods can be used: MDS (Multi-Dimensional Scaling) or hierarchical clustering methods (same as hclust). 

#pdf("temp/sampleRelation.pdf")          
#png("temp/sampleRelation.png")           
cvTh <- 0.1  #a threshold of coefficient of variance
plotSampleRelation(datalumiQN, cv.Th = cvTh, main=paste0("Clustering based on probes with sd/mean > ", cvTh))
dev.off()

##2. Detect the outlier 
##The current outlier detection is based on the distance from the sample to the center (average of all samples after removing 10 percent samples farthest away from the center).

pdf("temp/detectOutlier.pdf")
temp <- detectOutlier(datalumiQN, ifPlot=TRUE) 
dev.off()
temp <- detectOutlier(datalumiQN, ifPlot=FALSE)
any(temp) ##FALSE - none of the samples is outlier

##3. Clustering for the rows, in this case samples clustering, Euclidean Distance used here
dataExprs <- exprs(datalumiQN) # transform the lumi object into a data.frame (table)
tmp <- hclust(dist(t(dataExprs)),method="average")
pdf("temp/hclust.pdf")
#png("temp/hclustwithchip.png") #816 pixels = 8.5 inches
plot(tmp, xlab=paste0(ncol(dataExprs), " samples"))
dev.off()

##4. PCA (principal component analysis)
#dataExprs <- exprs(datalumiQN) # transform the lumi object into a data.frame (table)
clusteringPCA <- prcomp(t(dataExprs), scale=T)
summary(clusteringPCA)

marker <- factor(sample_barcode$marker)
levels(marker) <- 1:length(levels(marker))
markercolors <- palette(rainbow(length(levels(marker))))

pdf("temp/pca2comps.pdf")
plot(clusteringPCA$x[,1:2],col=markercolors[marker],pch=19,main="PCA")
legend("topright", levels(factor(sample_barcode$marker)), col=markercolors, pch=19)
dev.off()

library(scatterplot3d)
pdf("temp/pca3comps.pdf")
scatterplot3d(clusteringPCA$x[,1:3],color=markercolors[marker],pch=19,main="PCA")
legend("topleft", levels(factor(sample_barcode$marker)), col=markercolors, pch=19)
dev.off()


####################
#(3.3) Present counts
####################
## - Estimate or choose a threshold of detection p-values
## - A probe is present if its detection pvalue is less than or equal to the threshold. 
## - The probes that are absent in each sample should be removed for the analysis of differential expression.

##1. Retrieving the detection pvalues
#identical(detection(datalumiQN), detection(datalumi))
pvalues <- detection(datalumi)
class(pvalues)
pvalues[1:4,]

##boxplot
pdf("temp/pvaluesBoxplot.pdf")
par(mar=c(9,4,4,2))
boxplot(pvalues, las=2, cex.axis = 0.75) 
med_m <- apply(pvalues, 2, median) #median
h <- round(max(med_m), digits=4)   #max median
abline(h=h, lty=2, col=2)
axis(2, at=h, h, cex.axis = 0.75, col.axis = "red", las=1)
mtext(" max", side=4, cex = 0.75, col = "red", las=1, at=h)
h <- round(quantile(c(as.matrix(pvalues)), 0.50), digits=4) #50%-quantile of all p-values
abline(h=h, lty=2, col=2)
axis(2, at=h, h, cex.axis = 0.75, col.axis = "red", las=1)
mtext(" 50thQ", side=4, cex = 0.75, col = "red", las=1, at=h)
title("Detection p-values")
dev.off()


##2. Detection p-value threshold

##(1) estimate FDR using hisgram counts of the p-values
##Formula: FDR(x) = FP/(FP+TP) = n0*sum(p <= x)/sum(n[p <= x])
##  p - a vector of all detection p-values
##  x - a given p-value
##  n - counts of hisgram of p-values (i.e., p)
##  n0 - estimated average count of histogram of the absent probes' p-values (under H0), e.g., n0 = mean(n[p > 0.5]), the trimmed mean.

nbreaks <- 1e3 #or length(c(pvalues))/1000
f <- hist(c(pvalues), breaks = nbreaks, plot = F)

n <- f$counts
n0 <- mean(n[f$mids > 0.5])
FDR <- sapply(f$mids, function(p) n0*sum(f$mids <= p)/sum(n[f$mids <= p]))
max(f$mids[FDR <= 0.05]) #0.0435
detection.p.th <- max(f$mids[FDR <= 0.1]) #FDR <= 0.1
detection.p.th #is the detection threshold that has been calculated and that will be used to remove absent probes

plot(FDR, f$mids, type="l") 
abline(0, 1, col="gray")


##3. Estimate the present count of each gene (probe)
##the number of the samples that a probe presents 
##A probe is present if its detection pvalue is less than or equal to detection.p.th, a given threshold. 
##same as, rowSums(pvalues <= detection.p.th) 
presentCount <- detectionCall(datalumiQN, Th = detection.p.th) # number of probes in each samples that are called present.

##histogram of the present counts: 
##The histogram represents a frequency of the probes that are present in the same number of samples.
##The lower end represents the number of probes that have detection-pvalues > detection.p.th in all samples (PresentCount=0)
##The upper end represents the probes that have detection p-values <= detection.p.th in all samples
##To speed up the processing and reduce false positives, remove the unexpressed genes, that is, remove all probes for which presentCount=0

pdf("temp/hist_presentCount.pdf")
hist(presentCount, main=paste0("Histogram of PresentCount of Probes at a Threshold of ", detection.p.th), cex.main=0.85, xlab="Number of samples", breaks=c(-1, 0:max(presentCount)))
dev.off()
sum(presentCount == 0) #6167
sum(presentCount == ncol(datalumi)) #10912



##################################################################################################
###############      Differential expression (DE) analysis using limma package        ############
##################################################################################################



#(1) data filtering: unexpressed genes (probes) and sample(s) (such as technical replicates, irrelevant patients, ...)
#(2) model fitting for unpaired modelDesign (tests)
#(3) model fitting for paired modelDesign (tests)
#model fitting: 
#fit  - lmFit, modelDesign <- model.matrix(...)
#fit1 - contrasts.fit, contrast.matrix <- makeContrasts(...)
#fit2 - eBayes
#Question:
#In the data filtering, which should we do first, filtering genes or samples? 
#But the order does not matter if we give a fixed detection pvalue threshold.


###################
#Annotation file 
###################


##Read annotation file from the analysis folder
annotOriginal <- read.delim(annotationFileName, header = TRUE, sep = "\t", quote="\"", dec=".",fill = TRUE, comment.char="",skip=8, stringsAsFactors = FALSE)
dim(annotOriginal) #48240    28
head(annotOriginal)

annotOriginal$Definition <- gsub("Homo sapiens", "",  annotOriginal$Definition)

##The row names of data must be in the annotation file, that is, the column Array_Address_Id!!!
if (!all(rownames(datalumi) %in% annotOriginal$Array_Address_Id)){ print("some probes are not in the annotation file, please check that you are using the right annotation file")} # TRUE 

##extract the annotation file of the probes in the data
anno <- annotOriginal[annotOriginal$Array_Address_Id %in% rownames(datalumi), ]
dim(anno) #47323    28
any(is.na(anno$Array_Address_Id)) #[1] FALSE
rownames(anno) <- anno$Array_Address_Id
table(anno$Species)
head(anno)
# Homo sapiens ILMN Controls 
#        47231            92 

##combining annotations of interest and normalized expressions
##Notes:
##- The column of Array_Address_Id in the annotation data contains the unique probe IDs in the Illumina expression data.
##- The column of Probe_Id in the annotation data is not the same as the probe IDs in the Illumina expression data.
#normlized expressions
dataExprs <- exprs(datalumiQN)
dim(dataExprs) #47323    36

#annotations of interest
annotationIncluded <- c("RefSeq_ID", "Entrez_Gene_ID", "Symbol", "Probe_Id", "Array_Address_Id", "Probe_Sequence", "Definition")
if (!all(rownames(dataExprs) %in% anno$Array_Address_Id)){print("normalized data is not matching the annotation file")} #[1] TRUE
indx <- NULL
indx <- match(rownames(dataExprs), anno$Array_Address_Id)
anno.exprs <- cbind(anno[indx, annotationIncluded], dataExprs)
dim(anno.exprs) #47323    43
identical(rownames(anno.exprs), as.character(anno.exprs$Array_Address_Id)) #[1] TRUE
identical(rownames(anno.exprs), rownames(anno)[indx]) #[1] TRUE
head(anno.exprs)

filenm <- "temp/annotated_normalized_data"
write.table(anno.exprs, file=paste0(filenm, ".txt"), quote=FALSE, sep="\t", row.names = FALSE)  

##remove all variables that will not be used subsequently.
#Only "sample_barcode", "anno", "datalumi", and "datalumiQN" are needed for the analysis of differetial expression.
#ls()
#rm( list = setdiff(ls(), c("sample_barcode", "anno", "datalumi", "datalumiQN", "detection.p.th")) )
#ls()


################################################################
# Calculation differential expression (DE) between markers
################################################################

###################
#I. Data filtering
###################

####  groups included in the comparison
sample_group1 <- grep (paste(".*_", group1, "$", sep=""), colnames(dataExprs))
sample_group2 <- grep (paste(".*_", group2, "$", sep=""), colnames(dataExprs))
sample_included <- c(sample_group1, sample_group2)

## groups should only include paired samples
if (designType =="paired")
{
 dupPatient <- sample_barcode[ sample_barcode$marker %in% c(group1,group2) , ];
 dupPatient <- dupPatient$patient[duplicated(dupPatient$patient)];
 list <- NULL
 for ( i in 1:length(dupPatient))
 {
	list <- c(list,grep( paste0(dupPatient[i], "*") , colnames(dataExprs)))
 }
 sample_included <- sample_included[sample_included %in% list]
};
sample_included

#filter unexpressed (absent) genes using present counts estimated by detection pvalue threshold

##normalized expressions
dataExprs <- exprs(datalumiQN) # transform the lumi object in a data.frame
dim(dataExprs) #47323    36

#(1) presentCount estimated by detection.p.th 
#### appply the present count only on groups included in the comparison in case normalized data includes different tissue types
identical(colnames(datalumiQN), colnames(dataExprs))
presentCount <- detectionCall(datalumiQN[ , sample_included], Th = detection.p.th) 
head(presentCount)
if (!identical(names(presentCount), rownames(dataExprs))) {print("names present count should have been identical to normalized expression names")} #TRUE


#(3) filtering expressions 
dataf <- dataExprs[presentCount > 0, sample_included ] 
head(dataf)
colnames(dataf)
dim(dataf) #41156     9 #41156    33 #41156    36

samplef <- sample_barcode[sample_included  , ]
setequal(samplef$sampleName, colnames(dataf))  #[1] TRUE
identical(samplef$sampleName, colnames(dataf)) #[1] TRUE. So orders are also identical.
dim(samplef)

samplef

###################
#II. Model fittings
###################

#Testing differential expression (DE) between markers
#- unpaired design: in the case that the effects (mean) of patients are the same.
#model fittings: 
#fit  - lmFit, modelDesign <- model.matrix(...)
#fit1 - contrasts.fit, contrast.matrix <- makeContrasts(...)
#fit2 - eBayes


patient <- factor(samplef$patient)
table(patient)
marker <- factor(samplef$marker)
markerNames <- levels(marker)
table(marker)
chip <- factor(samplef$chip) 
table(chip)
variableX <- 0
if (length(factor(samplef$variableX)) != 0 )
variableX <- factor(samplef$variableX)

####if paired , shoud I remove other samples

if (sum(samplef$variableX ) != nrow(samplef) )
{
  if (designType == "paired") modelDesign <- model.matrix(~ 0 + marker + patient + variableX + chip) 
  if (designType == "unpaired") modelDesign <- model.matrix(~ 0 + marker +  variableX + chip) #not include an intercept.
}
if (sum(samplef$variableX ) == nrow(samplef) )
{
  if (designType == "paired") modelDesign <- model.matrix(~ 0 + marker + patient + chip) 
  if (designType == "unpaired") modelDesign <- model.matrix(~ 0 + marker +  chip) #not include an intercept.	
}

dim(modelDesign)
colnames(modelDesign)
fit <- lmFit(dataf, modelDesign) 

#(2) contrast fitting
marker_group1 <- paste("marker",  group1, sep="")
marker_group2 <- paste("marker",  group2, sep="")
contrastnm <- paste(marker_group1, "-", marker_group2, sep="") #contrast between a pair of markers  #HARD CODED
contrast.matrix <- makeContrasts(contrasts=contrastnm, levels=modelDesign)
contrast.matrix
fit1 <- contrasts.fit(fit, contrast.matrix)

#(3) eBayes fitting
fit2 <- eBayes(fit1)

##Mean Expressions for different marker 
#produce a list of means for each marker (grouped by a marker)
indx <- 1:length(marker)
meangrp_list <- tapply(indx, marker[indx], function(i) rowMeans(dataf[, i])) 
dim(meangrp_list)
meangrp <- sapply(meangrp_list, function(x) x) #changed as a matrix
dim(meangrp) #41156     3
colnames(meangrp) <- paste0("marker", names(meangrp_list))
head(meangrp)
class(meangrp)
logFC2 <- meangrp[,  c(paste0("marker", group1))] - meangrp[,  c(paste0("marker", group2))];
logFC2 <- as.matrix(logFC2);

##Annotations need to be included
annotationIncluded <- c("RefSeq_ID", "Entrez_Gene_ID", "Symbol", "Probe_Id", "Array_Address_Id", "Probe_Sequence", "Definition")

topfit <- topTable(fit2, number=nrow(dataf), adjust="BH")   
colnames(topfit)
head(topfit)

topfit <- topfit[, !(colnames(topfit) %in% c("B"))]
head(topfit)
dim(topfit)

rownames(topfit) <- topfit$ID
topID <- rownames(topfit)
all(rownames(dataf) %in% rownames(meangrp))
all(rownames(dataf) %in% anno$Array_Address_Id)
all(rownames(topfit) %in% anno$Array_Address_Id)
stopifnot( all(topID %in% rownames(meangrp)) ) #TRUE
stopifnot( all(topID %in% anno$Array_Address_Id) ) #TRUE


#combining annotations, mean expressions, and topfit 
indx <- NULL
indx <- match(topID, anno$Array_Address_Id)
anno.topfit <- cbind(anno[indx, annotationIncluded], meangrp[topID, c(marker_group1, marker_group2)], logFC2[topID,], topfit)
head(meangrp)

#add the normalized data for the two markers
indx <- NULL
indx <- match(topID, rownames(dataf)) # rownames(dataf) is array id
anno.topfit.dat <- cbind(anno.topfit, dataf[indx, ])


##########################################################
###### update annotation using the R biomaRt package #####
##########################################################

library(biomaRt)
mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes = getBM(attributes = c('illumina_humanht_12_v4', 'hgnc_symbol', 'description'), filters='illumina_humanht_12_v4', values=anno.topfit.dat$Probe_Id, mart=mart);
genes[genes==""] = NA;
head(genes)
names(genes)
names(genes) <- c("illumina_humanht_12_v4", "hgnc_symbol", "description")

sum(duplicated(genes$hgnc_symbol))
sum(duplicated(genes$description))
sum(duplicated(genes$illumina_humanht_12_v4))
genes$description = gsub("\\[Source.*", "", genes$description);
sum(is.na(genes$hgnc_symbol))

##duplicated probes
##- one probe ID might correspond to more than one gene that could be in different chromosomes!
##- collapse gene names for the duplicated probe ID
dupProbe <- unique( genes$illumina_humanht_12_v4[duplicated(genes$illumina_humanht_12_v4)] )
Gene <-  genes$hgnc_symbol
for (i in 1:length(dupProbe)) {
	indx <- (genes$illumina_humanht_12_v4 == dupProbe[i])
	Gene[indx] <- paste0(genes$hgnc_symbol[indx], collapse = "/")
	}
genes$upd_symbol <- Gene
rm(Gene)


myannotate = function(data)
{m = match(data$Probe_Id, genes$illumina_humanht_12_v4)
  data$description = genes[m,'description']
  data$GeneName = genes[m, 'upd_symbol']
  return(data)
};

anno.topfit.dat_annotated = myannotate(anno.topfit.dat)
head(anno.topfit.dat_annotated)


####################
# save the result table in the temp folder on the local computer
####################
filenm <- contrastnm
write.table(anno.topfit.dat_annotated, file=paste0("temp/", filenm,  "_", designType, ".txt"), quote=FALSE, sep="\t", row.names = FALSE)  
#rm(anno.topfit.dat)



###########################################
##  hierarchical clustering on the 2 groups 
###########################################
tmp <- hclust(dist(t(dataf)),method="average")
pdf(paste("temp/hclust", group1, group2, ".pdf",  sep="_"))
plot(tmp)
dev.off()
############

#############
#### STARTING FROM THIS POINT, TABLES ARE REDUCED AT THE GENE LEVEL BY SELECTION OF ONE PROBE CORRESPONDING TO THE BEST T VALUE 
############
anno.topfit.dat_annotated <- anno.topfit.dat_annotated[order(abs(anno.topfit.dat_annotated$t),decreasing=TRUE), ];
anno.topfit.dat_annotated <- anno.topfit.dat_annotated[!duplicated(anno.topfit.dat_annotated$GeneName) & !is.na(anno.topfit.dat_annotated$GeneName), ];
sum(duplicated(anno.topfit.dat_annotated$GeneName))

################
## Volcano plots
################
library(ggplot2)


head(anno.topfit.dat_annotated)
colnames(anno.topfit.dat_annotated)
##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
threshold = as.factor(abs(anno.topfit.dat_annotated$logFC) > 2 & anno.topfit.dat_annotated$adj.P.Val < 0.05)
 
##Construct the plot object
g = ggplot(data=anno.topfit.dat_annotated, aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position = "none") +
  xlim(c(min(anno.topfit.dat_annotated$logFC-0.5), max(anno.topfit.dat_annotated$logFC+0.5))) + ylim(c(0, max((-log10(anno.topfit.dat_annotated$adj.P.Val))+0.5))) +
  xlab("log2 fold change") + ylab("-log10 p-value")

anno.topfit.dat_annotated <- anno.topfit.dat_annotated[order(abs(anno.topfit.dat_annotated$t), decreasing=TRUE), ]
dd_text2 = data=anno.topfit.dat_annotated[(abs(anno.topfit.dat_annotated$logFC) > 2), ]
dd_text = anno.topfit.dat_annotated[ c(1:40),]

text  = geom_text(data=dd_text, aes(x=logFC, y=-log10(adj.P.Val), label=GeneName), colour="black", size=2, hjust=-0.1, vjust=-0.1)
text2  = geom_text(data=dd_text2, aes(x=logFC, y=-log10(adj.P.Val), label=GeneName), colour="black", size=2, hjust=-0.1, vjust=-0.1)
title <- labs(title=paste(contrastnm, designType))
pdf(paste0("temp/volcano_plot_", contrastnm, "_", designType, ".pdf"))
g + text + text2 + title
dev.off()

###############################################
## Heatmaps of the top 500 and top 50 genes
##############################################
library(gplots)

names(anno.topfit.dat_annotated)
head(anno.topfit.dat_annotated)
toptable500 <-anno.topfit.dat_annotated
#Re-order probes according to adj Pval (FDR) from smallest to largest
toptable500<-toptable500[order(toptable500$adj.P.Val,decreasing=FALSE),]
head(toptable500)

names(toptable500)
sample_group1 <- grep (paste(".*_", group1, "$", sep=""), colnames(toptable500))
sample_group2 <- grep (paste(".*_", group2, "$", sep=""), colnames(toptable500))

toptable500_mx<-as.matrix(toptable500[c(1:500),c(sample_group1, sample_group2)])
toptable500_df<-toptable500[c(1:500),];
names(toptable500_df)

names(toptable500_df)
matrix <- toptable500_mx;
myrownames <-toptable500_df$GeneName[1:500]
mycolnames <-colnames(toptable500_mx)
 
pdf(paste0("temp/heatmap_", group1, "_" ,group2,"_",  designType, "_", "top500.pdf"))
heatmap.2(matrix, col=bluered(69), scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.1,   Rowv=TRUE, Colv=TRUE, dendrogram = c("both"),
          labRow = myrownames, labCol=mycolnames,  cexCol=0.9, sepcolor="lightgray",  
          rowsep=FALSE, sepwidth=c(0.005,0.005), margins=c(9,9), 
          main=paste(group1, group2, designType, "top500", sep="_"));
dev.off()

#To create Excel table of toptable500 heatmap
hm <- heatmap.2(matrix, col=bluered(69), scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.1,   Rowv=TRUE, Colv=TRUE, dendrogram = c("both"),
          labRow = myrownames, labCol=mycolnames,  cexCol=0.9, sepcolor="lightgray",  
          rowsep=FALSE, sepwidth=c(0.005,0.005), margins=c(9,9), 
          main=paste(group1, group2, designType, "top500", sep="_"));

head(toptable500_df)
ordered_hm_df <- cbind(toptable500_df[ , c(1:15) ]  , toptable500_mx[ , hm$colInd])
ordered_hm_df <- ordered_hm_df[rev(hm$rowInd)  , ]

sample_group1 <- grep (paste(".*_", group1, "$", sep=""), colnames(ordered_hm_df))
sample_group2 <- grep (paste(".*_", group2, "$", sep=""), colnames(ordered_hm_df))

ordered_hm_df$mean <- apply(ordered_hm_df[ , c(sample_group1, sample_group2)], 1, mean);
ordered_hm_df$sd <- apply(ordered_hm_df[ , c(sample_group1, sample_group2)], 1, sd);
newdf <- apply(ordered_hm_df[ , c(sample_group1, sample_group2)], 2, function(x) ((x- ordered_hm_df$mean) / ordered_hm_df$sd))

ordered_hm_df_scaled <- cbind(ordered_hm_df, newdf)      
write.csv(ordered_hm_df_scaled, paste0("temp/heatmap_", group1, "_" ,group2,"_",  designType, "_", "top500.csv"))  ;          

sample_group1 <- grep (paste(".*_", group1, "$", sep=""), colnames(toptable500))
sample_group2 <- grep (paste(".*_", group2, "$", sep=""), colnames(toptable500))

## TOP50
toptable50_mx<-as.matrix(toptable500[c(1:50),c(sample_group1, sample_group2)])
toptable50_df<-toptable500[c(1:50),];
matrix <- toptable50_mx;
myrownames <-toptable50_df$GeneName[1:50]
mycolnames <-colnames(toptable50_mx)
pdf(paste0("temp/heatmap_", group1, "_" ,group2,"_",  designType, "_", "top50.pdf"))
heatmap.2(matrix, col=bluered(69), scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5,   Rowv=TRUE, Colv=TRUE, dendrogram = c("both"),
          labRow = myrownames, labCol=mycolnames,  cexCol=0.9, sepcolor="lightgray",  
          rowsep=FALSE, sepwidth=c(0.005,0.005), margins=c(9,9), 
          main=paste(group1, group2, designType, "top50", sep="_"));
dev.off()


#####################################################
## Preparing rank file and expression file for GSEA #
#####################################################

#### prepare rank files for GSEA ####
myrank = function(data)
{data = data[order(abs(data$t),decreasing=TRUE), ];
 data = data[!duplicated(data$GeneName) & !is.na(data$GeneName), ];
 data = data[order(data$t, decreasing=TRUE), ];
 rank = data[ , c("GeneName", "t")];
 return(rank)
}

anno.topfit.dat_rank = myrank(anno.topfit.dat_annotated)

write.table(anno.topfit.dat_rank, paste0("temp/", contrastnm, "_", designType, ".RNK"), sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE);

#### create expression file ####
head(anno.topfit.dat_annotated)

sample_group1 <- grep (paste(".*_", group1, "$", sep=""), colnames(anno.topfit.dat_annotated))
sample_group2 <- grep (paste(".*_", group2, "$", sep=""), colnames(anno.topfit.dat_annotated))

make_expression_file = function(data)
{data = data[order(abs(data$t),decreasing=TRUE), ];
 data = data[!duplicated(data$GeneName) & !is.na(data$GeneName), ];
 data = data[order(data$t, decreasing=TRUE), ];
 data = cbind( data[ , c("GeneName","description")] ,  data[ , c(sample_group1, sample_group2)]    )
 return(data);#edited the Description   
}

anno.topfit.dat_expression = make_expression_file(anno.topfit.dat_annotated) 

# write expression file
write.table(anno.topfit.dat_expression, paste0("temp/", contrastnm, "_", designType, "_expression.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

### save the workspace
save.image(paste0("temp/", contrastnm, ".RData"))

### remove the temp output folder with the name of the comparison
### the output folder should not exist
if (file.exists("temp")) { file.rename("temp",  paste0(contrastnm, designType, "_output"))};






