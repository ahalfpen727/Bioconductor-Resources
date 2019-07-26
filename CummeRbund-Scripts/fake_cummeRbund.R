#Starting from the output of the RNA-seq Tutorial Part 1.

#Install packages and load libraries
#install.packages("ggplot2")
library(ggplot2);library(gplots)

#If X11 not available, open a pdf device for output of all plots
pdf(file="fake_cummeRbund_output.pdf")

#To learn about any command type: ?command_name  OR  help.search("command_name")
#e.g.    ?read.table

#demo(graphics)

#Clean up workspace - i.e. delete variable created by the graphics demo
rm(list = ls(all = TRUE))

#List the variables that exist in your current work space
ls()

#### Import the gene expression data from the Tophat/Cufflinks/Cuffdiff tutorial

#Set working directory where results files exist
working_dir = "~/workspace/rnaseq/de/tophat_cufflinks/ref_only"
setwd(working_dir)

#List the current contents of this directory - it is empty right now so it will be displayed as 'character(0)'
dir()

#Import expression and differential expression results from the Bowtie/Samtools/Tophat/Cufflinks/Cuffdiff pipeline
file1="isoforms.read_group_tracking"
file2="isoform_exp.diff"
file3="isoforms.fpkm_tracking"
file4="cds.read_group_tracking"
file5="cds_exp.diff"
file6="cds.fpkm_tracking"
file7="genes.read_group_tracking"
file8="gene_exp.diff"
file9="genes.fpkm_tracking"
file10="tss_groups.read_group_tracking"
file11="tss_group_exp.diff"
file12="tss_groups.fpkm_tracking"

#Read in tab delimited files and assign the resulting 'dataframe' to a variable
#Use 'as.is' for columns that contain text/character values (i.e. non-numerical values)
all_fpkm = read.table(file1, header=TRUE, sep="\t", as.is=c(1:2,9))
tn_de = read.table(file2, header=TRUE, sep="\t", as.is=c(1:7,14))
tn_fpkm = read.table(file3, header=TRUE, sep="\t", as.is=c(1:9,13,17))
allc_fpkm = read.table(file4, header=TRUE, sep="\t", as.is=c(1:2,9))
tnc_de = read.table(file5, header=TRUE, sep="\t", as.is=c(1:7,14))
tnc_fpkm = read.table(file6, header=TRUE, sep="\t", as.is=c(1:9,13,17))
allg_fpkm = read.table(file7, header=TRUE, sep="\t", as.is=c(1:2,9))
tng_de = read.table(file8, header=TRUE, sep="\t", as.is=c(1:7,14))
tng_fpkm = read.table(file9, header=TRUE, sep="\t", as.is=c(1:9,13,17))
allt_fpkm = read.table(file10, header=TRUE, sep="\t", as.is=c(1:2,9))
tnt_de = read.table(file11, header=TRUE, sep="\t", as.is=c(1:7,14))
tnt_fpkm = read.table(file12, header=TRUE, sep="\t", as.is=c(1:9,13,17))


#### Working with 'dataframes'
#View the first five rows of data (all columns) in one of the dataframes created
#library(biomaRt);library(DOSE);library(org.Hs.eg.db);library(clusterProfiler)
#geneDE<-as.data.frame(cbind(SYMBOL=tng_de$gene,log2foldchange=tng_de$log2.fold_change,qvalue=tng_de$q_value))
#hg38.entrez.de <- bitr(geneDE$SYMBOL, fromType = "SYMBOL",toType =c("ENTREZID"),OrgDb = org.Hs.eg.db)
#hg38_de.df<-merge(hg38.entrez.de, geneDE)
#hg38_DE.df<-as.matrix(hg38_de.df, sep="\t")
#hg38_de.df<-hg38_de.df[,-1]
#hg38.de.df<-file.path("hg38_entrez_DE.df")
#write.table(hg38_de.df, file=hg38.de.df, sep = "\t", row.names = F, col.names = T,quote = F)
#Get the first 3 rows of data and a selection of columns
#tn_de[1:3,c(2:4,7,10,12)]
#Do the same thing, but using the column names instead of numbers
#tn_de[1:3, c("gene_id","locus","value_1","value_2")]

#Rename some of the columns from ugly names to more human readable names
names(all_fpkm) = c("tracking_id", "condition", "replicate", "raw_frags", "internal_scaled_frags", "external_scaled_frags", "FPKM", "effective_length", "status")
names(tn_de) = c("test_id", "gene_id", "gene_name", "locus", "sample_1", "sample_2", "status", "value_1", "value_2", "fold_change", "test_stat", "p_value", "q_value", "significant")
names(tn_fpkm) = c("tracking_id", "class_code", "nearest_ref_id", "gene_id", "gene_name", "tss_id", "locus", "length", "coverage", "LUTS_FPKM", "LUTS_conf_lo", "LUTS_conf_hi", "LUTS_status", "CTRL_FPKM", "CTRL_conf_lo", "CTRL_conf_hi", "CTRL_status")
names(allc_fpkm) = c("tracking_id", "condition", "replicate", "raw_frags", "internal_scaled_frags", "external_scaled_frags", "FPKM", "effective_length", "status")
names(tnc_de) = c("test_id", "gene_id", "gene_name", "locus", "sample_1", "sample_2", "status", "value_1", "value_2", "fold_change", "test_stat", "p_value", "q_value", "significant")
names(tnc_fpkm) = c("tracking_id", "class_code", "nearest_ref_id", "gene_id", "gene_name", "tss_id", "locus", "length", "coverage", "LUTS_FPKM", "LUTS_conf_lo", "LUTS_conf_hi", "LUTS_status", "CTRL_FPKM", "CTRL_conf_lo", "CTRL_conf_hi", "CTRL_status")
names(allg_fpkm) = c("tracking_id", "condition", "replicate", "raw_frags", "internal_scaled_frags", "external_scaled_frags", "FPKM", "effective_length", "status")
names(tng_de) = c("test_id", "gene_id", "gene_name", "locus", "sample_1", "sample_2", "status", "value_1", "value_2", "fold_change", "test_stat", "p_value", "q_value", "significant")
names(tng_fpkm) = c("tracking_id", "class_code", "nearest_ref_id", "gene_id", "gene_name", "tss_id", "locus", "length", "coverage", "LUTS_FPKM", "LUTS_conf_lo", "LUTS_conf_hi", "LUTS_status", "CTRL_FPKM", "CTRL_conf_lo", "CTRL_conf_hi", "CTRL_status")
names(allt_fpkm) = c("tracking_id", "condition", "replicate", "raw_frags", "internal_scaled_frags", "external_scaled_frags", "FPKM", "effective_length", "status")
names(tnt_de) = c("test_id", "gene_id", "gene_name", "locus", "sample_1", "sample_2", "status", "value_1", "value_2", "fold_change", "test_stat", "p_value", "q_value", "significant")
names(tnt_fpkm) = c("tracking_id", "class_code", "nearest_ref_id", "gene_id", "gene_name", "tss_id", "locus", "length", "coverage", "LUTS_FPKM", "LUTS_conf_lo", "LUTS_conf_hi", "LUTS_status", "CTRL_FPKM", "CTRL_conf_lo", "CTRL_conf_hi", "CTRL_status")

#Get ID to gene name mapping
isoform_mapping=tn_fpkm[,"gene_name"]
names(isoform_mapping)=tn_fpkm[,"tracking_id"]

gene_mapping=tng_fpkm[,"gene_name"]
names(gene_mapping)=tng_fpkm[,"tracking_id"]

cds_mapping=tnc_fpkm[,"gene_name"]
names(cds_mapping)=tnc_fpkm[,"tracking_id"]

tss_mapping=tnt_fpkm[,"gene_name"]
names(tss_mapping)=tnt_fpkm[,"tracking_id"]


#Reformat per-replicate gene FPKM data into a standard matrix
LUTS_1=allg_fpkm[allg_fpkm[,"condition"]=="LUTS" & allg_fpkm[,"replicate"]==0,"FPKM"]
LUTS_2=allg_fpkm[allg_fpkm[,"condition"]=="LUTS" & allg_fpkm[,"replicate"]==1,"FPKM"]
LUTS_3=allg_fpkm[allg_fpkm[,"condition"]=="LUTS" & allg_fpkm[,"replicate"]==2,"FPKM"]
LUTS_4=allg_fpkm[allg_fpkm[,"condition"]=="LUTS" & allg_fpkm[,"replicate"]==3,"FPKM"]
LUTS_5=allg_fpkm[allg_fpkm[,"condition"]=="LUTS" & allg_fpkm[,"replicate"]==4,"FPKM"]
LUTS_6=allg_fpkm[allg_fpkm[,"condition"]=="LUTS" & allg_fpkm[,"replicate"]==5,"FPKM"]
LUTS_7=allg_fpkm[allg_fpkm[,"condition"]=="LUTS" & allg_fpkm[,"replicate"]==6,"FPKM"]
LUTS_8=allg_fpkm[allg_fpkm[,"condition"]=="LUTS" & allg_fpkm[,"replicate"]==7,"FPKM"]

CTRL_1=allg_fpkm[allg_fpkm[,"condition"]=="CTRL" & allg_fpkm[,"replicate"]==0,"FPKM"]
CTRL_2=allg_fpkm[allg_fpkm[,"condition"]=="CTRL" & allg_fpkm[,"replicate"]==1,"FPKM"]
CTRL_3=allg_fpkm[allg_fpkm[,"condition"]=="CTRL" & allg_fpkm[,"replicate"]==2,"FPKM"]
CTRL_4=allg_fpkm[allg_fpkm[,"condition"]=="CTRL" & allg_fpkm[,"replicate"]==3,"FPKM"]
CTRL_5=allg_fpkm[allg_fpkm[,"condition"]=="CTRL" & allg_fpkm[,"replicate"]==4,"FPKM"]
CTRL_6=allg_fpkm[allg_fpkm[,"condition"]=="CTRL" & allg_fpkm[,"replicate"]==5,"FPKM"]
CTRL_7=allg_fpkm[allg_fpkm[,"condition"]=="CTRL" & allg_fpkm[,"replicate"]==6,"FPKM"]
CTRL_8=allg_fpkm[allg_fpkm[,"condition"]=="CTRL" & allg_fpkm[,"replicate"]==7,"FPKM"]

ids=unique(allg_fpkm[,"tracking_id"])
gene_names=gene_mapping[ids]
gene.ids=unique(allg_fpkm[,"tracking_id"])
gene_names=gene_mapping[gene.ids]
gene.fpkm_matrix=data.frame(gene.ids, LUTS_1,LUTS_2,LUTS_3,LUTS_4,LUTS_5,LUTS_6,LUTS_7,LUTS_8,CTRL_1,CTRL_2,CTRL_3,CTRL_4,CTRL_5,CTRL_6,CTRL_7,CTRL_8)
row.names(gene.fpkm_matrix)=gene.ids
data_columns=c(2:17)
short_names=c("LUTS_1","LUTS_2","LUTS_3","LUTS_4","LUTS_5","LUTS_6","LUTS_7","LUTS_8","CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5","CTRL_6","CTRL_7","CTRL_8")
geneLUTSgroup<-as.data.frame(gene.fpkm_matrix[,c("LUTS_1","LUTS_2","LUTS_3","LUTS_4","LUTS_5","LUTS_6","LUTS_7","LUTS_8")])
geneCTRLgroup<-as.data.frame(gene.fpkm_matrix[,c("CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5","CTRL_6","CTRL_7","CTRL_8")])

#Reformat per-replicate CDS FPKM data into a standard matrix
LUTS_1=allc_fpkm[allc_fpkm[,"condition"]=="LUTS" & allc_fpkm[,"replicate"]==0,"FPKM"]
LUTS_2=allc_fpkm[allc_fpkm[,"condition"]=="LUTS" & allc_fpkm[,"replicate"]==1,"FPKM"]
LUTS_3=allc_fpkm[allc_fpkm[,"condition"]=="LUTS" & allc_fpkm[,"replicate"]==2,"FPKM"]
LUTS_4=allc_fpkm[allc_fpkm[,"condition"]=="LUTS" & allc_fpkm[,"replicate"]==3,"FPKM"]
LUTS_5=allc_fpkm[allc_fpkm[,"condition"]=="LUTS" & allc_fpkm[,"replicate"]==4,"FPKM"]
LUTS_6=allc_fpkm[allc_fpkm[,"condition"]=="LUTS" & allc_fpkm[,"replicate"]==5,"FPKM"]
LUTS_7=allc_fpkm[allc_fpkm[,"condition"]=="LUTS" & allc_fpkm[,"replicate"]==6,"FPKM"]
LUTS_8=allc_fpkm[allc_fpkm[,"condition"]=="LUTS" & allc_fpkm[,"replicate"]==7,"FPKM"]

CTRL_1=allc_fpkm[allc_fpkm[,"condition"]=="CTRL" & allc_fpkm[,"replicate"]==0,"FPKM"]
CTRL_2=allc_fpkm[allc_fpkm[,"condition"]=="CTRL" & allc_fpkm[,"replicate"]==1,"FPKM"]
CTRL_3=allc_fpkm[allc_fpkm[,"condition"]=="CTRL" & allc_fpkm[,"replicate"]==2,"FPKM"]
CTRL_4=allc_fpkm[allc_fpkm[,"condition"]=="CTRL" & allc_fpkm[,"replicate"]==3,"FPKM"]
CTRL_5=allc_fpkm[allc_fpkm[,"condition"]=="CTRL" & allc_fpkm[,"replicate"]==4,"FPKM"]
CTRL_6=allc_fpkm[allc_fpkm[,"condition"]=="CTRL" & allc_fpkm[,"replicate"]==5,"FPKM"]
CTRL_7=allc_fpkm[allc_fpkm[,"condition"]=="CTRL" & allc_fpkm[,"replicate"]==6,"FPKM"]
CTRL_8=allc_fpkm[allc_fpkm[,"condition"]=="CTRL" & allc_fpkm[,"replicate"]==7,"FPKM"]

cds.ids=unique(allc_fpkm[,"tracking_id"])
cds_names=gene_mapping[cds.ids]
cds.fpkm_matrix=data.frame(cds.ids, LUTS_1,LUTS_2,LUTS_3,LUTS_4,LUTS_5,LUTS_6,LUTS_7,LUTS_8,CTRL_1,CTRL_2,CTRL_3,CTRL_4,CTRL_5,CTRL_6,CTRL_7,CTRL_8)
row.names(cds.fpkm_matrix)=cds.ids
data_columns=c(2:17)
short_names=c("LUTS_1","LUTS_2","LUTS_3","LUTS_4","LUTS_5","LUTS_6","LUTS_7","LUTS_8","CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5","CTRL_6","CTRL_7","CTRL_8")
cdsLUTSgroup<-as.data.frame(cds.fpkm_matrix[,c("LUTS_1","LUTS_2","LUTS_3","LUTS_4","LUTS_5","LUTS_6","LUTS_7","LUTS_8")])
cdsCTRLgroup<-as.data.frame(cds.fpkm_matrix[,c("CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5","CTRL_6","CTRL_7","CTRL_8")])

#Reformat per-replicate tss group FPKM data into a standard matrix
LUTS_1=allt_fpkm[allt_fpkm[,"condition"]=="LUTS" & allt_fpkm[,"replicate"]==0,"FPKM"]
LUTS_2=allt_fpkm[allt_fpkm[,"condition"]=="LUTS" & allt_fpkm[,"replicate"]==1,"FPKM"]
LUTS_3=allt_fpkm[allt_fpkm[,"condition"]=="LUTS" & allt_fpkm[,"replicate"]==2,"FPKM"]
LUTS_4=allt_fpkm[allt_fpkm[,"condition"]=="LUTS" & allt_fpkm[,"replicate"]==3,"FPKM"]
LUTS_5=allt_fpkm[allt_fpkm[,"condition"]=="LUTS" & allt_fpkm[,"replicate"]==4,"FPKM"]
LUTS_6=allt_fpkm[allt_fpkm[,"condition"]=="LUTS" & allt_fpkm[,"replicate"]==5,"FPKM"]
LUTS_7=allt_fpkm[allt_fpkm[,"condition"]=="LUTS" & allt_fpkm[,"replicate"]==6,"FPKM"]
LUTS_8=allt_fpkm[allt_fpkm[,"condition"]=="LUTS" & allt_fpkm[,"replicate"]==7,"FPKM"]

CTRL_1=allt_fpkm[allt_fpkm[,"condition"]=="CTRL" & allt_fpkm[,"replicate"]==0,"FPKM"]
CTRL_2=allt_fpkm[allt_fpkm[,"condition"]=="CTRL" & allt_fpkm[,"replicate"]==1,"FPKM"]
CTRL_3=allt_fpkm[allt_fpkm[,"condition"]=="CTRL" & allt_fpkm[,"replicate"]==2,"FPKM"]
CTRL_4=allt_fpkm[allt_fpkm[,"condition"]=="CTRL" & allt_fpkm[,"replicate"]==3,"FPKM"]
CTRL_5=allt_fpkm[allt_fpkm[,"condition"]=="CTRL" & allt_fpkm[,"replicate"]==4,"FPKM"]
CTRL_6=allt_fpkm[allt_fpkm[,"condition"]=="CTRL" & allt_fpkm[,"replicate"]==5,"FPKM"]
CTRL_7=allt_fpkm[allt_fpkm[,"condition"]=="CTRL" & allt_fpkm[,"replicate"]==6,"FPKM"]
CTRL_8=allt_fpkm[allt_fpkm[,"condition"]=="CTRL" & allt_fpkm[,"replicate"]==7,"FPKM"]

tss.ids=unique(allt_fpkm[,"tracking_id"])
tss_names=gene_mapping[tss.ids]
tss.fpkm_matrix=data.frame(tss.ids, LUTS_1,LUTS_2,LUTS_3,LUTS_4,LUTS_5,LUTS_6,LUTS_7,LUTS_8,CTRL_1,CTRL_2,CTRL_3,CTRL_4,CTRL_5,CTRL_6,CTRL_7,CTRL_8)
row.names(tss.fpkm_matrix)=tss.ids
data_columns=c(2:17)
short_names=c("LUTS_1","LUTS_2","LUTS_3","LUTS_4","LUTS_5","LUTS_6","LUTS_7","LUTS_8","CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5","CTRL_6","CTRL_7","CTRL_8")
tssLUTSgroup<-as.data.frame(tss.fpkm_matrix[,c("LUTS_1","LUTS_2","LUTS_3","LUTS_4","LUTS_5","LUTS_6","LUTS_7","LUTS_8")])
tssCTRLgroup<-as.data.frame(tss.fpkm_matrix[,c("CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5","CTRL_6","CTRL_7","CTRL_8")])

#Reformat per-replicate isoform FPKM data into a standard matrix
LUTS_1=all_fpkm[all_fpkm[,"condition"]=="LUTS" & all_fpkm[,"replicate"]==0,"FPKM"]
LUTS_2=all_fpkm[all_fpkm[,"condition"]=="LUTS" & all_fpkm[,"replicate"]==1,"FPKM"]
LUTS_3=all_fpkm[all_fpkm[,"condition"]=="LUTS" & all_fpkm[,"replicate"]==2,"FPKM"]
LUTS_4=all_fpkm[all_fpkm[,"condition"]=="LUTS" & all_fpkm[,"replicate"]==3,"FPKM"]
LUTS_5=all_fpkm[all_fpkm[,"condition"]=="LUTS" & all_fpkm[,"replicate"]==4,"FPKM"]
LUTS_6=all_fpkm[all_fpkm[,"condition"]=="LUTS" & all_fpkm[,"replicate"]==5,"FPKM"]
LUTS_7=all_fpkm[all_fpkm[,"condition"]=="LUTS" & all_fpkm[,"replicate"]==6,"FPKM"]
LUTS_8=all_fpkm[all_fpkm[,"condition"]=="LUTS" & all_fpkm[,"replicate"]==7,"FPKM"]

CTRL_1=all_fpkm[all_fpkm[,"condition"]=="CTRL" & all_fpkm[,"replicate"]==0,"FPKM"]
CTRL_2=all_fpkm[all_fpkm[,"condition"]=="CTRL" & all_fpkm[,"replicate"]==1,"FPKM"]
CTRL_3=all_fpkm[all_fpkm[,"condition"]=="CTRL" & all_fpkm[,"replicate"]==2,"FPKM"]
CTRL_4=all_fpkm[all_fpkm[,"condition"]=="CTRL" & all_fpkm[,"replicate"]==3,"FPKM"]
CTRL_5=all_fpkm[all_fpkm[,"condition"]=="CTRL" & all_fpkm[,"replicate"]==4,"FPKM"]
CTRL_6=all_fpkm[all_fpkm[,"condition"]=="CTRL" & all_fpkm[,"replicate"]==5,"FPKM"]
CTRL_7=all_fpkm[all_fpkm[,"condition"]=="CTRL" & all_fpkm[,"replicate"]==6,"FPKM"]
CTRL_8=all_fpkm[all_fpkm[,"condition"]=="CTRL" & all_fpkm[,"replicate"]==7,"FPKM"]

ids=unique(all_fpkm[,"tracking_id"])
gene_names=gene_mapping[ids]
gene.ids=unique(all_fpkm[,"tracking_id"])
iso.fpkm_matrix=data.frame(gene.ids, LUTS_1,LUTS_2,LUTS_3,LUTS_4,LUTS_5,LUTS_6,LUTS_7,LUTS_8,CTRL_1,CTRL_2,CTRL_3,CTRL_4,CTRL_5,CTRL_6,CTRL_7,CTRL_8)
row.names(iso.fpkm_matrix)=gene.ids
data_columns=c(2:17)
short_names=c("LUTS_1","LUTS_2","LUTS_3","LUTS_4","LUTS_5","LUTS_6","LUTS_7","LUTS_8","CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5","CTRL_6","CTRL_7","CTRL_8")
isoLUTSgroup<-as.data.frame(iso.fpkm_matrix[,c("LUTS_1","LUTS_2","LUTS_3","LUTS_4","LUTS_5","LUTS_6","LUTS_7","LUTS_8")])
isoCTRLgroup<-as.data.frame(iso.fpkm_matrix[,c("CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5","CTRL_6","CTRL_7","CTRL_8")])

#Assign colors to each.  You can specify color by RGB, Hex code, or name
#To get a list of color names:
colours()
data_colors=c("tomato1","tomato2","tomato3","royalblue1","royalblue2","royalblue3")

#View expression values for the transcripts of a particular gene symbol of chromosome 1.  e.g. 'TST'
i = which(gene.fpkm_matrix[,"gene.ids"] == "CXCL12")
gene.fpkm_matrix[i,]

#What if we want to view values for a list of genes of interest all at once?
genes_of_interest = c("CXCL12", "CXCR4","COL1A1","MMP9","TGFB1","TGFB2")
i = which(gene.fpkm_matrix[,"gene.ids"] %in% genes_of_interest)
gene.fpkm_matrix[i,]
tss.fpkm_matrix[i,]
cds.fpkm_matrix[i,]
iso.fpkm_matrix[i,]
#What if we want to view values for a list of genes of interest all at once?

#### Examine basic features of the differential expression file
#In part 1 of the tutorial, cuffdiff attempted to perform a differential expression test for each row of data (i.e. each gene/transcript)
#However, sometimes this test fails due to insufficient data, etc.  These cases are summarized in the 'status' column
#Summarize the status of all tests
gene_status_counts=table(tng_de[,"status"])
gene_status_counts
#Plot #1 - Make a barplot of these status counts, first using the basic plotting functions of R, and then using the ggplot2 package
barplot(gene_status_counts, col=rainbow(6), xlab="Status", ylab="Transcript count", main="Gene status counts reported by Cuffdiff")

Gene_Status=factor(tng_de[,"status"])
qplot(Gene_Status, data=tng_de, geom="bar", fill=Gene_Status, xlab="Status", ylab="Transcript count", main="Gene status counts reported by Cuffdiff")


CDS_status_counts=table(tnc_de[,"status"])
CDS_status_counts
#Plot #2 - Now the same idea using ggplot2
CDS_Status=factor(tnc_de[,"status"])
qplot(CDS_Status, data=tnc_de, geom="bar", fill=CDS_Status, xlab="Status", ylab="Transcript count", main="CDS status counts reported by Cuffdiff")

isoform_status_counts=table(tn_de[,"status"])
isoform_status_counts
ISO_Status=factor(tn_de[,"status"])
qplot(ISO_Status, data=tn_de, geom="bar", fill=ISO_Status, xlab="Status", ylab="Transcript count", main="Isoform status counts reported by Cuffdiff")


TSS_Status=factor(tnt_de[,"status"])
tss_status_counts=table(tnt_de[,"status"])
qplot(TSS_Status, data=tnt_de, geom="bar", fill=TSS_Status, xlab="Status", ylab="Transcript count", main="TSS-group status counts reported by Cuffdiff")



#Plot #3 - Make a piechart of these status counts, first using the basic plotting functions of R, and then using the ggplot2 package
xvf<-pie(gene_status_counts, col=rainbow(6), main="Status counts reported by Cuffdiff")


#Plot #4 - Make a dotchart of these status counts
dotchart(as.numeric(gene_status_counts), col=rainbow(6), labels=names(gene_status_counts), xlab="Transcript count", main="Gene status transcript counts reported by Cuffdiff", pch=16)

#Each row of data represents a transcript. Many of these transcripts represent the same gene. Determine the numbers of transcripts and unique genes
length(tn_de[,"gene_name"]) #Transcript count
length(unique(tn_de[,"gene_name"])) #Unique Gene count

length(tnc_de[,"gene_name"]) #Transcript count
length(unique(tnc_de[,"gene_name"])) #Unique Gene count

length(tng_de[,"gene_name"]) #Transcript count
length(unique(tng_de[,"gene_name"])) #Unique Gene count


#### Plot #5 - the number of transcripts per gene.
#Many genes will have only 1 transcript, some genes will have several transcripts
#Use the 'table()' command to count the number of times each gene symbol occurs (i.e. the # of transcripts that have each gene symbol)
#Then use the 'hist' command to create a histogram of these counts
#How many genes have 1 transcript?  More than one transcript?  What is the maximum number of transcripts for a single gene?
counts=table(tnt_de[,"gene_name"])
c_one = length(which(counts == 1))
c_more_than_one = length(which(counts > 1))
c_max = max(counts)
hist(counts, breaks=50, col="bisque4", xlab="TSS Transcripts per gene", main="Distribution of TSS transcript count per gene")
legend_text = c(paste("Genes with one TSS transcript =", c_one), paste("Genes with more than one TSS transcript =", c_more_than_one), paste("Max TSS transcripts for single gene = ", c_max))
legend("topright", legend_text, lty=NULL)
plot(counts, col="bisque4", xlab="Transcripts per gene", main="Distribution of TSS transcript count per gene")

counts=table(tnc_de[,"gene_name"])
c_one = length(which(counts == 1))
c_more_than_one = length(which(counts > 1))
c_max = max(counts)
hist(counts, breaks=50, col="bisque4", xlab="CDS transcripts per gene", main="Distribution of CDS transcript count per gene")
legend_text = c(paste("Genes with one CDS transcript =", c_one), paste("Genes with more than CDS one transcript =", c_more_than_one), paste("Max CDS transcripts for single gene = ", c_max))
legend("topright", legend_text, lty=NULL)
plot(counts, col="bisque4", xlab="Transcripts per gene", main="Distribution of CDS transcript count per gene")

counts=table(tn_de[,"gene_name"])
c_one = length(which(counts == 1))
c_more_than_one = length(which(counts > 1))
c_max = max(counts)
hist(counts, breaks=50, col="bisque4", xlab="Isoform Transcripts per gene", main="Distribution of isoform transcript count per gene")
legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
legend("topright", legend_text, lty=NULL)
plot(counts, col="bisque4", xlab="Transcripts per gene", main="Distribution of isoform transcript count per gene")

#### Plot #7 - the distribution of transcript sizes as a histogram
#In this analysis we supplied Cufflinks with transcript models so the lengths will be those of known transcripts
#However, if we had used a de novo transcript discovery mode, this step would give us some idea of how well transcripts were being assembled
#If we had a low coverage library, or other problems, we might get short 'transcripts' that are actually only pieces of real transcripts
hist(tn_fpkm[,"length"], breaks=50, xlab="Transcript length (bp)", main="Distribution of Isoform transcript lengths", col="steelblue")

#hist(tnc_fpkm[,"length"], breaks=50, xlab="Transcript length (bp)", main="Distribution of CDS transcript lengths", col="steelblue")
#hist(tng_fpkm[,"length"], breaks=50, xlab="Transcript length (bp)", main="Distribution of Gene transcript lengths", col="steelblue")
#hist(tnt_fpkm[,"length"], breaks=50, xlab="Transcript length (bp)", main="Distribution of TSS transcript lengths", col="steelblue")


#### Summarize FPKM values for all 6 replicates
#What are the minimum and maximum FPKM values for a particular library?
min(fpkm_matrix[,"LUTS_8"])
max(fpkm_matrix[,"LUTS_8"])

#Set the minimum non-zero FPKM values for use later.
#Do this by grabbing a copy of all data values, coverting 0's to NA, and calculating the minimum or all non NA values
zz = iso.fpkm_matrix[,data_columns]
zz[zz==0] = NA
min_nonzero = min(zz, na.rm=TRUE)
boxplot(log2(iso.fpkm_matrix[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for all Isoform libraries")

zz = tss.fpkm_matrix[,data_columns]
zz[zz==0] = NA
min_nonzero = min(zz, na.rm=TRUE)
boxplot(log2(tss.fpkm_matrix[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for all TSS libraries")

zz = cds.fpkm_matrix[,data_columns]
zz[zz==0] = NA
min_nonzero = min(zz, na.rm=TRUE)
boxplot(log2(cds.fpkm_matrix[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for all CDS libraries")

zz = gene.fpkm_matrix[,data_columns]
zz[zz==0] = NA
min_nonzero = min(zz, na.rm=TRUE)
boxplot(log2(gene.fpkm_matrix[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for all Gene libraries")
#### Plot #8 - View the range of values and general distribution of FPKM values for all libraries
#Create boxplots for this purpose
#Note that the bold horizontal line on each boxplot is the median

#### Plot #9 - plot a pair of replicates to assess reproducibility of technical replicates
#Tranform the data by converting to log2 scale after adding an arbitrary small value to avoid log2(0)
min_nonzero=.01
x = gene.fpkm_matrix[,"LUTS_1"]
y = gene.fpkm_matrix[,"LUTS_2"]
plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (LUTS, Replicate 1)", ylab="FPKM (LUTS, Replicate 2)", main="Comparison of expression values LUTs_1 and LUTS_2")
#Add a straight line of slope 1, and intercept 0
abline(a=0,b=1)
#Calculate the correlation coefficient and display in a legend
rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")

x = iso.fpkm_matrix[,"LUTS_2"]
y = iso.fpkm_matrix[,"LUTS_3"]
plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (LUTS, Replicate 2)", ylab="FPKM (LUTS, Replicate 3)", main="Comparison of expression values LUTS_2 and LUTS_3")
#Add a straight line of slope 1, and intercept 0
abline(a=0,b=1)
#Calculate the correlation coefficient and display in a legend
rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")

x = gene.fpkm_matrix[,"LUTS_3"]
y = gene.fpkm_matrix[,"LUTS_4"]
plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (LUTS, Replicate 3)", ylab="FPKM (LUTS, Replicate 4)", main="Comparison of expression values LUTs_3 and LUTS_4")
#Add a straight line of slope 1, and intercept 0
abline(a=0,b=1)
#Calculate the correlation coefficient and display in a legend
rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")

LUTSgroup<-as.data.frame(gene.fpkm_matrix[,c("LUTS_1","LUTS_2","LUTS_3","LUTS_4","LUTS_5","LUTS_6","LUTS_7","LUTS_8")])
CTRLgroup<-as.data.frame(gene.fpkm_matrix[,c("CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5","CTRL_6","CTRL_7","CTRL_8")])

x = CTRLgroup[,c("CTRL_1","CTRL_2","CTRL_3","CTRL_4")]
y = CTRLgroup[,c("CTRL_5","CTRL_6","CTRL_7","CTRL_8")]
cor(x,y, method ="pearson")

x = LUTSgroup[,c("LUTS_1","LUTS_2","LUTS_3","LUTS_4")]
y = LUTSgroup[,c("LUTS_5","LUTS_6","LUTS_7","LUTS_8")]
cor(x,y, method ="pearson")

#### Plot #10 - Scatter plots with a large number of data points can be misleading ... regenerate this figure as a density scatter plot
colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
smoothScatter(x=log2(x+min_nonzero), y=log2(y+min_nonzero), xlab="FPKM (LUTS, Replicate 1)", ylab="FPKM (LUTS, Replicate 2)", main="Comparison of expression values for a pair of replicates", colramp=colors, nbin=200)


#### Plot all sets of replicates on a single plot
#Create an function that generates an R plot.  This function will take as input the two libraries to be compared and a plot name and color
plotCor = function(lib1, lib2, name, color){
  x=gene.fpkm_matrix[,lib1]
  y=gene.fpkm_matrix[,lib2]
  zero_count = length(which(x==0)) + length(which(y==0))
  plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col=color, cex=0.25, xlab=lib1, ylab=lib2, main=name)
  abline(a=0,b=1)
  rs=cor(x,y, method="pearson")^2
  legend_text = c(paste("R squared = ", round(rs, digits=3), sep=""), paste("Zero count = ", zero_count, sep=""))
  legend("topleft", legend_text, lwd=c(1,NA), col="black", bg="white", cex=0.8)
}
#Open a plotting page with room for two plots on one page
par(mfrow=c(1,2))

#Plot #11 - Now make a call to our custom function created above, once for each library comparison
plotCor("LUTS_1", "CTRL_1", "LUTS_1 vs CTRL_1", "tomato2")
plotCor("LUTS_2", "CTRL_2", "LUTS_2 vs CTRL_2", "royalblue2")


##### One problem with these plots is that there are so many data points on top of each other, that information is being lost
#Regenerate these plots using a density scatter plot
plotCor2 = function(lib1, lib2, name, color){
  x=gene.fpkm_matrix[,lib1]
  y=gene.fpkm_matrix[,lib2]
  zero_count = length(which(x==0 & y==0)) # & length(which(y==0))
  counts_in_one = length(which(x > 0 | y > 0)) # | length(which(y > 0))
  counts_in_both = length(which(x > 0 & y > 0)) # | length(which(y > 0))
  total_count = length(x)
  colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  smoothScatter(x=log2(x+min_nonzero), y=log2(y+min_nonzero), xlab=lib1, ylab=lib2, main=name, colramp=colors, nbin=275)
  abline(a=0,b=1)
  rs=cor(x,y, method="pearson")^2
  legend_text = c(paste("R squared = ", round(rs, digits=3), sep=""), paste("Zero count = ", zero_count, sep=""), paste("Counts in One Group = ", counts_in_one, sep=""), paste("Counts in Both Groups = ", counts_in_both, sep=""), paste("Total Counts = ", total_count, sep=""))
  legend("topleft", legend_text, lwd=c(1,NA), col="black", bg="white", cex=0.8)
}

#### Compare the correlation 'distance' between all replicates
par(mfrow=c(1,2))
plotCor2("LUTS_1", "CTRL_1", "LUTS_1 vs CTRL_1", "tomato2")
plotCor2("LUTS_1", "CTRL_2", "LUTS_1 vs CTRL_2", "tomato2")

par(mfrow=c(1,2))
plotCor2("LUTS_1", "CTRL_3", "LUTS_1 vs CTRL_3", "royalblue2")
plotCor2("LUTS_1", "CTRL_4", "LUTS_1 vs CTRL_4", "royalblue2")

par(mfrow=c(1,2))
plotCor2("LUTS_1", "CTRL_5", "LUTS_1 vs CTRL_5", "royalblue2")
plotCor2("LUTS_1", "CTRL_6", "LUTS_1 vs CTRL_6", "royalblue2")

par(mfrow=c(1,2))
plotCor2("LUTS_1", "CTRL_7", "LUTS_1 vs CTRL_7", "royalblue2")
plotCor2("LUTS_1", "CTRL_8", "LUTS_1 vs CTRL_8", "royalblue2")

par(mfrow=c(1,2))
plotCor2("LUTS_3", "CTRL_4", "LUTS_3 vs CTRL_4", "tomato2")
plotCor2("LUTS_4", "CTRL_5", "LUTS_4 vs CTRL_5", "royalblue2")

par(mfrow=c(1,2))
plotCor2("LUTS_6", "CTRL_7", "LUTS_6 vs CTRL_7", "tomato2")
plotCor2("LUTS_7", "CTRL_8", "LUTS_7 vs CTRL_8", "royalblue2")

#### Compare the correlation 'distance' between all replicates
par(mfrow=c(1,2))
plotCor2("LUTS_1", "LUTS_2", "LUTS_1 vs LUTS_2", "royalblue2")
plotCor2("LUTS_1", "LUTS_3", "LUTS_1 vs LUTS_3", "royalblue2")

par(mfrow=c(1,2))
plotCor2("LUTS_1", "LUTS_3", "LUTS_1 vs LUTS_3", "tomato2")
plotCor2("LUTS_1", "LUTS_4", "LUTS_1 vs LUTS_4", "tomato2")

par(mfrow=c(1,2))
plotCor2("LUTS_1", "LUTS_5", "LUTS_1 vs LUTS_5", "royalblue2")
plotCor2("LUTS_1", "LUTS_6", "LUTS_1 vs LUTS_6", "royalblue2")

par(mfrow=c(1,2))
plotCor2("LUTS_1", "LUTS_7", "LUTS_1 vs LUTS_7", "royalblue2")
plotCor2("LUTS_1", "LUTS_8", "LUTS_1 vs LUTS_8", "royalblue2")

#### Compare the correlation 'distance' between all replicates
#Do we see the expected pattern for all eight libraries (i.e. replicates most similar, then tumor vs. normal)?
plotCors = function(lib1,lib2, name, color){
  x=gene.fpkm_matrix[,lib1]
  y=gene.fpkm_matrix[,lib2]
  counts = length(x)
  zero_count_Luts = length(which(x==0))
  zero_count_Ctrl = length(which(y==0))
  colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  smoothScatter(x=log2(x+min_nonzero), y=log2(y+min_nonzero), xlab=lib1, ylab=lib2, main=name, colramp=colors, nbin=275)
  abline(a=0,b=1)
  rs=cor(x,y, method="pearson")^2
  legend_text = c(paste("R squared = ", round(rs, digits=3), sep=""),paste("Counts = ", counts, sep=""), paste("Zero Counts LUTS= ", zero_count_Luts, sep=""), paste("Zero Counts CTRL= ", zero_count_Ctrl, sep=""))
  legend("topleft", legend_text, lwd=c(1,NA), col="black", bg="white", cex=0.8)
}
#### Compare the correlation 'distance' between all replicates
par(mfrow=c(1,2))
plotCors("LUTS_1", "CTRL_1", "LUTS_1 vs CTRL_1", "tomato2")
plotCors("LUTS_1", "CTRL_2", "LUTS_1 vs CTRL_2", "tomato2")

par(mfrow=c(1,2))
plotCors("LUTS_1", "CTRL_3", "LUTS_1 vs CTRL_3", "royalblue2")
plotCors("LUTS_1", "CTRL_4", "LUTS_1 vs CTRL_4", "royalblue2")

par(mfrow=c(1,2))
plotCors("LUTS_1", "CTRL_5", "LUTS_1 vs CTRL_5", "royalblue2")
plotCors("LUTS_1", "CTRL_6", "LUTS_1 vs CTRL_6", "royalblue2")

par(mfrow=c(1,2))
plotCors("LUTS_1", "CTRL_7", "LUTS_1 vs CTRL_7", "royalblue2")
plotCors("LUTS_1", "CTRL_8", "LUTS_1 vs CTRL_8", "royalblue2")

par(mfrow=c(1,2))
plotCors("LUTS_3", "CTRL_4", "LUTS_3 vs CTRL_4", "tomato2")
plotCors("LUTS_4", "CTRL_5", "LUTS_4 vs CTRL_5", "royalblue2")

par(mfrow=c(1,2))
plotCors("LUTS_6", "CTRL_7", "LUTS_6 vs CTRL_7", "tomato2")
plotCors("LUTS_7", "CTRL_8", "LUTS_7 vs CTRL_8", "royalblue2")

#### Compare the correlation 'distance' between all replicates
par(mfrow=c(1,2))
plotCors("LUTS_1", "LUTS_2", "LUTS_1 vs LUTS_2", "royalblue2")
plotCors("LUTS_1", "LUTS_3", "LUTS_1 vs LUTS_3", "royalblue2")

par(mfrow=c(1,2))
plotCors("LUTS_1", "LUTS_3", "LUTS_1 vs LUTS_3", "tomato2")
plotCors("LUTS_1", "LUTS_4", "LUTS_1 vs LUTS_4", "tomato2")

par(mfrow=c(1,2))
plotCors("LUTS_1", "LUTS_5", "LUTS_1 vs LUTS_5", "royalblue2")
plotCors("LUTS_1", "LUTS_6", "LUTS_1 vs LUTS_6", "royalblue2")

par(mfrow=c(1,2))
plotCors("LUTS_1", "LUTS_7", "LUTS_1 vs LUTS_7", "royalblue2")
plotCors("LUTS_1", "LUTS_8", "LUTS_1 vs LUTS_8", "royalblue2")

#Calculate the FPKM sum for all libraries
gene.fpkm_matrix[,"sum"]=apply(gene.fpkm_matrix[,data_columns], 1, sum)
#Identify the genes with a grand sum FPKM of at least 5 - we will filter out the genes with very low expression across the board
i = which(gene.fpkm_matrix[,"sum"] > 5)
#Calculate the correlation between all pairs of data
r=cor(gene.fpkm_matrix[i,data_columns], use="pairwise.complete.obs", method="pearson")
#Print out these correlation values
r

#Calculate the FPKM sum for all libraries
gene.fpkm_matrix[,"sum"]=apply(gene.fpkm_matrix[,data_columns2], 1, sum)
#Identify the genes with a grand sum FPKM of at least 5 - we will filter out the genes with very low expression across the board
i = which(gene.fpkm_matrix[,"sum"] > 5)
#Calculate the correlation between all pairs of data
r=cor(gene.fpkm_matrix[i,data_columns2], use="pairwise.complete.obs", method="pearson")
#Print out these correlation values
r

#### Plot #13 - Convert correlation to 'distance', and use 'multi-dimensional scaling' to display the relative differences between libraries
#This step calculates 2-dimensional coordinates to plot points for each library
#Libraries with similar expression patterns (highly correlated to each other) should group together
#What pattern do we expect to see, given the types of libraries we have (technical replicates, biologal replicates, tumor/normal)?
d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
par(mfrow=c(1,1))
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes)", xlim=c(-0.12,0.12), ylim=c(-0.12,0.12))
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], short_names, col=data_colors)

#### Plot #14 - View the distribution of differential expression values as a histogram
#Display only those that are significant according to Cuffdiff
sig = which(tn_de[,"p_value"]<0.05)
de = log2(tn_de[sig,"value_1"]+min_nonzero) - log2(tn_de[sig,"value_2"]+min_nonzero)
tn_de[,"de"] = log2(tn_de[,"value_1"]+min_nonzero) - log2(tn_de[,"value_2"]+min_nonzero)
hist(de, breaks=50, col="seagreen", xlab="Log2 difference (LUTS - CTRL)", main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topleft", "Fold-change > 4", lwd=2, lty=2)
r=cor(gene.fpkm_matrix[i,data_columns], use="pairwise.complete.obs", method="pearson")

sig = which(tnc_de[,"p_value"]<0.05)
cde = log2(tnc_de[sig,"value_1"]+min_nonzero) - log2(tnc_de[sig,"value_2"]+min_nonzero)
tnc_de[,"cde"] = log2(tnc_de[,"value_1"]+min_nonzero) - log2(tnc_de[,"value_2"]+min_nonzero)
hist(cde, breaks=50, col="seagreen", xlab="Log2 difference (LUTS - CTRL)", main="Distribution of differential CDS expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topleft", "Fold-change > 4", lwd=2, lty=2)

sig = which(tng_de[,"p_value"]<0.05)
gde = log2(tng_de[sig,"value_1"]+min_nonzero) - log2(tng_de[sig,"value_2"]+min_nonzero)
tng_de[,"gde"] = log2(tng_de[,"value_1"]+min_nonzero) - log2(tng_de[,"value_2"]+min_nonzero)
hist(gde, breaks=50, col="seagreen", xlab="Log2 difference (LUTS - CTRL)", main="Distribution of differential gene expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topleft", "Fold-change > 4", lwd=2, lty=2)


#### Plot #15 - Display the grand expression values from UHR and HBR and mark those that are significantly differentially expressed
x=log2(tn_de[,"value_1"]+min_nonzero)
y=log2(tn_de[,"value_2"]+min_nonzero)
plot(x=x, y=y, pch=16, cex=0.25, xlab="LUTS FPKM (log2)", ylab="CTRL FPKM (log2)", main="LUTS vs CTRL FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

x=log2(tnc_de[,"value_1"]+min_nonzero)
y=log2(tnc_de[,"value_2"]+min_nonzero)
plot(x=x, y=y, pch=16, cex=0.25, xlab="LUTS FPKM (log2)", ylab="CTRL FPKM (log2)", main="LUTS vs CTRL FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

x=log2(tng_de[,"value_1"]+min_nonzero)
y=log2(tng_de[,"value_2"]+min_nonzero)
plot(x=x, y=y, pch=16, cex=0.25, xlab="LUTS FPKM (log2)", ylab="CTRL FPKM (log2)", main="LUTS vs CTRL FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

#Get the gene symbols for the top N (according to corrected p-value) and display them on the plot
topn = order(abs(tn_de[,"fold_change"]), decreasing=TRUE)[1:25]
topn = order(tn_de[,"q_value"])[1:25]
text(x[topn], y[topn], tn_de[topn,"gene_name"], col="black", cex=0.75, srt=45)

topn = order(abs(tnc_de[,"fold_change"]), decreasing=TRUE)[1:25]
topn = order(tnc_de[,"q_value"])[1:25]
text(x[topn], y[topn], tnc_de[topn,"gene_name"], col="black", cex=0.75, srt=45)

topn = order(abs(tng_de[,"fold_change"]), decreasing=TRUE)[1:25]
topn = order(tng_de[,"q_value"])[1:25]
text(x[topn], y[topn], tng_de[topn,"gene_name"], col="black", cex=0.75, srt=45)


#### Write a simple table of differentially expressed transcripts to an output file
#Each should be significant with a log2 fold-change >= 2
sig = which(tn_de[,"p_value"]<0.05 & abs(tn_de[,"de"]) >= 2)
sig_tn_de = tn_de[sig,]

sigc = which(tnc_de[,"p_value"]<0.05 & abs(tnc_de[,"de"]) >= 2)
sigc_tn_de = tnc_de[sigc,]

sigg = which(tng_de[,"p_value"]<0.05 & abs(tng_de[,"de"]) >= 2)
sigg_tn_de = tng_de[sigg,]

#Order the output by or p-value and then break ties using fold-change
o = order(sig_tn_de[,"q_value"], -abs(sig_tn_de[,"de"]), decreasing=FALSE)
output = sig_tn_de[o,c("gene_id","gene_name","locus","value_1","value_2","de","p_value")]
write.table(output, file="SigDE_supplementary_R.txt", sep="\t", row.names=FALSE, quote=FALSE)

o = order(sigc_tn_de[,"q_value"], -abs(sigc_tn_de[,"de"]), decreasing=FALSE)
output = sigc_tn_de[o,c("gene_id","gene_name","locus","value_1","value_2","de","p_value")]
write.table(output, file="SigDE_supplementary_R.txt", sep="\t", row.names=FALSE, quote=FALSE)

o = order(sigg_tn_de[,"q_value"], -abs(sigg_tn_de[,"de"]), decreasing=FALSE)
output = sigg_tn_de[o,c("gene_id","gene_name","locus","value_1","value_2","de","p_value")]
write.table(output, file="SigDE_supplementary_R.txt", sep="\t", row.names=FALSE, quote=FALSE)


#View selected columns of the first 25 lines of output
output[1:25,c(2,4,5,6,7)]

#You can open the file "SigDE.txt" in Excel, Calc, etc.
#It should have been written to the current working directory that you set at the beginning of the R tutorial
dir()


#### Plot #16 - Create a heatmap to vizualize expression differences between the eight samples
#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

main_title="sig DE Transcripts"
par(cex.main=0.8)
sig_genes=tn_de[sig,"test_id"]
sig_gene_names=gene_mapping[sig_genes]
data=log2(as.matrix(gene.fpkm_matrix[sig_genes,data_columns])+1)
heatmap.2(data, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(6,7), Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, cexRow=0.3, cexCol=1, labRow=sig_gene_names,col=rev(heat.colors(75)))

dev.off()

#The output file can be viewed in your browser at the following url:
#Note, you must replace cbw## with your own amazon instance number (e.g., "cbw01"))
#http://__YOUR_IP_ADDRESS__/workspace/rnaseq/de/tophat_cufflinks/ref_only/Tutorial_Part3_Supplementary_R_output.pdf
#To exit R type:
quit(save="no")
