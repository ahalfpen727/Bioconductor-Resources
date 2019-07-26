# This file contains all instructions described in Contreras et al. (2017) to create gene networks from RNA-seq data.
# The instructions are designed by Biolinux version 8 but it can be applied to other unix based systems.
# Please note that data, methods and parameters used at each step are meant to be for demonstration purposes
# and by no means should be taken as the only way or as a general rule to do data analyses in all cases.
# Changes in methods and parameters can have a significant impact on your final results and should be carefully
# evaluated and decided upon depending on your scientific aims as well as experimental design.
#
# Last version June 15, 2017.


# Preparing a working directory before you start:

#To create a working directory, select a location in your disk and create a folder using the following command:

mkdir WorkDir

#To move to your newly created folder, press:

cd WorkDir

#If you want to find the complete path to the working directory, you can use the command:

pwd

# To facilitate further processing, the path of the created directory is created.

export WD=$(pwd)

#Download and install the tools used in this chapter

#The next step is to create a folder for all the tools to be used. To do so, create a directory called "tools" inside your "WorkDir" folder. 

mkdir tools

cd tools

#First, the user must install FastQC tool. Consider that the program version may change in the future. We will use wget to download the program and unzip to decompress the downloaded file.

wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip

unzip fastqc_v0.11.5.zip

#To execute the FastQC tool, the user must change the folder permissions:

chmod 777 FastQC/fastqc

#Then, the user must install Trimmomatic, following similar steps as described above. 

wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip

unzip Trimmomatic-0.36.zip

#Finally, download the Hisat2 align tool

wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.5-Linux_x86_64.zip

unzip hisat2-2.0.5-Linux_x86_64.zip

#To return to the "WorkDir" folder, use:

cd $WD

# RNA-seq data acquisition from public databases

#To download the FASTQ files used in this example, the user can use a graphical web interface or a file generated to download the data by command line available at http://virtualplant.bio.puc.cl/share/CHAPTER/fastq.url.txt. 

# To use the provided script, move to "WorkDir" directory and type the following command: 

cat fastq.url.txt  |xargs -n1 wget $1


#3.5 Data quality assessment

#To run FastQC for every file we will list all the FASTQ symbolic links and run the FastQC for each list element:

# To accelerate data processing, "xargs" command allows to select the number of threads used to process the data in parallel. This number should be selected based on the number of processors available in your computer to perform the analysis. For example, the above command line can be performed in parallel using 2 threads with the options -P: 
 
ls *.fastq.gz | xargs -n1 -P2 $WD/tools/FastQC/fastqc $1

#The FastQC output are html files which will be saved in the working directory and can be visualized with any browser.

# 3.6 Data trimming


# To run Trimmomatic, it is necessary to create a text file "script.trimmo.sh" containing the instructions to be executed. Make sure you are in the "WorkDir" folder and press:

echo 'nohup java -jar $WD/tools/Trimmomatic-0.36/trimmomatic-0.36.jar SE $1 $1.trim.fil.gz LEADING:20 TRAILING:20 AVGQUAL:25 SLIDINGWINDOW:10:30 MINLEN:36  > $1.trim.nohup' > script.trimmo.sh

#Then, we will select all the FASTQ files and execute Trimmomatic on each one. "–P2" indicates the number of threads that will be perform in parallel (2 in this case): 

ls *.fastq.gz |xargs -n1 -P2 sh script.trimmo.sh


#3.7 Aligning trimmed RNA-seq data to reference genome

# Download Assembly and annotation files from the Arabidopsis thaliana genome from Arabidopsis Information Portal (www.araport.org) website. To do this, register in the web and go to Data > Downloads and download TAIR10_Chr.all.fasta.gz file from TAIR10_genome_release/assembly. Then, go to  TAIR10_genome_release/annotation and download Araport11_GFF3_genes_transposons.201606.gff and Araport11_GFF3_genes_transposons.201606.gtf.

#The downloaded files must be decompressed with "gunzip" to continue working with them. The following commands do this:

gunzip TAIR10_Chr.all.fasta.gz

gunzip Araport11_GFF3_genes_transposons.201606.gff.gz

#The following command generates the necessary files to perform the alignment with HISAT2:

$WD/tools/hisat2-2.0.5/hisat2-build TAIR10_Chr.all.fasta TAIR10_Chr.all  

#In the same way as done before, create a simple text file named "script.align.sh" containing the instructions for the single-end reads alignment:

echo 'nohup $WD/tools/hisat2-2.0.5/hisat2 -x TAIR10_Chr.all -U $1 -S $1.sam > $1.align.stat.txt' > script.align.sh

#This instruction creates a sam file for each ".trim.fil.gz" file. The file must be save in the "WorkDir" directory. The following command executes the previous instructions for each trimmed file using 2 cores. 

ls *.trim.fil.gz | xargs -n1 -P10 sh script.align.sh

# The "align.stat.txt" files contains the stats of the alignment. The alignments files with ".sam" extension are the SAM files. Sequence Alignment/Map (SAM) format is a generic alignment text file that describes the genome coordinates were the reads were aligned, the number of possible match sites in the genome and other related information.

# 3.8 Identification of expressed genes from RNA-seq data.
 
#For the identification of expressed genes, we will use only protein coding genes in the Arabidopsis genome. The user can select the genes with the following command:

grep "locus_type=protein_coding" Araport11_GFF3_genes_transposons.201606.gff |cut -f9 |cut -d';' -f1 |cut -d'=' -f2 >protein_coding


#	To assign the reads to a gene we use R software (R Core Team 2015). In this guide we will use the Rsubread package that facilitates the RNA-seq read data analyses, proportionating several metrics: quality assessment of sequence reads, read alignment, read summarization, exon-exon junction among others (Liao, Smyth, & Shi, 2013). In the following instructions, lines that contain commands to execute in R will be indicated with a greater-than sign ">". 

#Start R in the terminal:

R

# Or open Rstudio. Then, make sure that the working directory is "WorkDir". To install Rsubread package, press the following instructions:

source("http://bioconductor.org/biocLite.R")

biocLite("Rsubread")

library(Rsubread)

 
#Then, it is necessary to run the feature counts for all SAM files. First, select the SAM files:

sam.list <- dir(pattern=".sam")

#Since there are different types of library preparation and run configurations, the user should use the appropriate one. In this example, we use the non-stranded (strandSpecific=0) option. For strand-specific library preparations, the user may refer to provided supplementary files.  

#The following command runs the program using 2 threads and store the output in "fc0" object. 

fc0 <- featureCounts( sam.list, annot.ext= "Araport11_GFF3_genes_transposons.201606.gtf",isGTFAnnotationFile=T, allowMultiOverlap=T, isPairedEnd=F, nthreads=2, strandSpecific=0)

#Then, we save the counts and associated stats in tab delimited files:

write.table(fc0$counts, "fc0.counts.txt", sep="\t", col.names=NA, quote=F)

write.table(fc0$stat, "fc0.stat.txt", sep="\t", row.names=F, quote=F)

#Finally, we select protein coding genes with reads counts:

protein_coding<-as.matrix(read.table("protein_coding"))

counts<-fc0$counts[protein_coding,]

counts <- counts[rowSums(counts)>0,]


#3.9 Determination of differentially expressed genes (DEGs) from RNA-seq data using DEseq2.

source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)

#DESeq2 requires the raw count data from the FASTQ files. The user can select the columns corresponding to the wild type plants under the desired treatments, i.e. from the column 20 to 25, named respectively SRR420813, SRR420814, SRR420815, SRR420816, SRR420817 and SRR420818. 

#From these data, generate a table:

select <- counts[,20:25]

#Then, the user must create an object describing the conditions (Control or Treatment) associated with each column in "select" object. In this case, the first three columns of select are "control" (or "c"), and the latter three are "treatments" (or "t"). The position of the column should match the condition position as follow:

condition <- c("c","c","c","t","t","t")

coldata <- data.frame(row.names=colnames(select), condition)

#The next step is to generate an object containing the counts and conditions to be compared as entry for DESeq2:

dds <- DESeqDataSetFromMatrix( countData=select[rowSums(select)>0,],  colData=coldata, design=~condition)

#Run the DESeq2 function to start the statistical analysis. It retrieves a normalized expression level in logarithmic base 2 scale.

results <- results(DESeq(dds))

#Filter the data with the desired thresholds of fold change level and adjusted P-values. In this case fold change is log2FC >1 with adjusted P-values <0.01. Store the data in "filter" object:

res <- na.exclude(as.data.frame(results))

filter <- res[(abs(res$log2FoldChange)>1 & res$padj<0.01),]

#The resulting table can be written as text file for further analyzes:

write.table(filter,"regulated.txt", quote=F,sep="\t",  col.names = NA)

#This table contains all regulated genes with a differential expression level more than 2 times compared to the control and an adjusted p-value < 0.01. In this example we obtained 231 regulated genes.

#3.10 Normalization of gene reads counts from RNA-seq data. 


#The following commands can perform the mentioned normalization. 

source("http://bioconductor.org/biocLite.R")

biocLite("EBSeq")

library(EBSeq)

NormData <- GetNormalizedMat(counts, MedianNorm(counts))

#The resulting gene expression matrix contains unique row identifiers and row counts obtained from different experiments on each column. After normalization, it is recommended to delete very low counts data or sum a unit to all data in order to avoid values equal to zero. Then is useful to generate a logarithmic matrix of the data to standardize the variance.

NormData.log <- log2(NormData+1)

#3.11 Calculating correlation between genes pairs

#Since the goal of this guide is to generate a gene co-expression network, it is necessary to determinate the genes that co-express in our determined conditions. First, the user should select the regulated genes from the normalized table generated previously by the DESeq2 and extract their normalized counts, using the following command:

Norm.interest <- NormData.log[rownames(filter),]

install.packages("https://cran.r-project.org/src/contrib/Archive/psych/psych_1.7.5.tar.gz", repos=NULL, type="source")

library("psych")

Norm.interest.corr <- corr.test( t(Norm.interest), method="pearson", ci=F)

#Among the many results of this function, there are two triangular matrices with the needed data. One matrix contains the correlation values and the other contain the p-values in the upper part and the adjusted p-values in the lower part. To generate a table comprising the data organized properly to build the network, the user should type the following commands: 


Norm.interest.corr$p[lower.tri( Norm.interest.corr$p,diag=TRUE)]=NA

Pval.adj <- as.data.frame(as.table(Norm.interest.corr$p))

Norm.interest.corr$r [lower.tri( Norm.interest.corr$r,diag=TRUE)]=NA

Correlation <- as.data.frame(as.table(Norm.interest.corr$r))

Cor.table <- na.exclude(cbind( Correlation, Pval.adj))[,c(1,2,3,6)]

colnames(Cor.table) <- c("gene1","gene2","cor","p.adj")

#The generated table, can be filer according absolute correlation (0.9) and adjusted p-value (0.01) thresholds:

Cor.table.filt <- Cor.table [(abs(Cor.table[,3])>0.9 & Cor.table[,4] <0.01 ),]

#To facilitate graphs generation in future network, adjusted p-value can be converted to logarithmic base 10 and to avoid problems with log of 0, we propose to sum the minimal adjusted p-value found to the entire column.  
  
Log.p.adj <- log10(Cor.table.filt[,4]+min(Cor.table.filt [Cor.table.filt[,4]!=0,4]))

Cor.table.filt<-cbind(Cor.table.filt,Log.p.adj)

write.table(Cor.table.filt, "Cor.table.filter.txt", sep="\t", row.names=F, quote=F)

#Now we have generated a table containing the statistically significant correlations for every pair of gene differentially regulated.

#3.12 Network Stats

#The basics stats of the network, degree and betweenness, can be calculated using "igraph" R package. 

 install.packages("igraph")
 
 library(igraph)
 
 g <- graph.data.frame( Cor.table.filt[,1:2], directed=FALSE) 
 
 degree <- degree(g)
 
 betweenness <- betweenness(g)
 
 Node_nw_st <- data.frame( degree, betweenness)

#To integrate this two parameters, a combined ranking can be generated and added to the former table:

 Rank_stat <- rowMeans(cbind(rank(Node_nw_st[,1]), rank(Node_nw_st[,2])))

 Node_nw_st <- cbind(Node_nw_st, Rank_stat)

 write.table(Node_nw_st,file="Node_nw_st.txt", sep="\t", col.names = NA,quote=F)

