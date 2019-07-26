#  ## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("mypackage")
## ---- eval=FALSE----------------------------------------------------
suppressPackageStartupMessages(library(VariantAnnotation))
suppressMessages(suppressPackageStartupMessages(library(cgdv17)))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
library(Rsubread)
library("systemPipeR") # Loads the package
library(help="systemPipeR") # Lists package info
#vignette("systemPipeR") # Opens vignetteppress     RsubreadUsersGuide()
library(googleVis)
library(GenomicFeatures)
library(BiocParallel)
library(cummeRbund)
library(Rsamtools)

M<-gvisMotionChart(Fruits, "Fruit", "Year")
plot(M)
cat(M$html$chart, file = "tmp.html")
head(Fruits)

# Get Genome and genes
txdb=makeTranscriptDbFromUCSC(genome='hg19',tablename='ensGene')
ex_by_gene=exonsBy(txdb,'gene')
reads=readBamGappedAlignments("aligned_reads_sorted.bam")

# The annotations have chromosomes called
names(seqlengths(tx_by_gene))

# The reads have chromosomes called
 as.character(unique(rname(reads)))  

 #The annotations have chromosomes called
names(seqlengths(tx_by_gene))
   
# tx19 = TxDb.Hsapiens.UCSC.hg19.knownGene
#The reads have chromosomes called
as.character(unique(rname(reads)))
new_read_chr_names=gsub("(.*)[T]*\\..*","chr\\1",rname(reads))
reads=GRanges(seqnames=new_read_chr_names,ranges=IRanges(start=start(reads),
			  end=end(reads)), strand=rep("*",length(reads)))

 #If we just want to convert chromosome names.

reads=GRanges(seqnames=new_read_chr_names,ranges=IRanges(start=start(reads),
			  end=end(reads)), strand=strand(reads))

# If we just want to make each read be ambiguous with respect to strand.

 reads=GRanges(seqnames=rname(reads),ranges=IRanges(start=start(reads),
 				end=end(reads)),strand=rep("*",length(reads)))
 
# Finally, we get the number of reads that overlap with each gene
#  (or whatever else you're interested in).

counts=countOverlaps(tx_by_gene,reads)
toc=data.frame(condition1_rep1=counts1.1,condition1_rep2=counts1.2,
     condition2_rep1=counts2.1,condition2_rep2=counts2.2,stringsAsFactors=FALSE)
rownames(toc)=names(tx_by_gene)

# Normalizing
library(edgeR)
   
norm_factors=calcNormFactors(as.matrix(toc))

# Statistical testing
DGE=DGEList(toc,lib.size=norm_factors*colSums(toc),
            group=rep(c("Condition1","Condition2"),c(2,2)))
disp=estimateCommonDisp(DGE) 
tested=exactTest(disp)
   
# Gene Set testing (GO)
library(goseq)
genes = as.integer(p.adjust(tested$table$p.value, method = "BH") < 0.05)
names(genes) = row.names(tested$table)

# calculate a probability weighting function, correcting for length bias, 
# a technical bias present in all forms of RNA-seq data.[25]

pwf=nullp(genes,'hg19','ensGene')

#  to correct for total read count bias 
#  we would calculate the pwf using the number of counts from each gene as follows:

pwf=nullp(genes,bias.data=rowsum(counts[match(names(genes),rownames(counts))]))

# p-value for each GO category being over represented amongst DE genes.

GO.pvals=goseq(pwf,'hg19','ensGene')

   
# Processing in R

# load the sorted, indexed, BAM files into R. A
# unstranded RNA-seq, make the strand designator for each read ambiguous
#  (which is done by setting it to "*"). After starting R we run,
library(Rsamtools)
   #Create a list of bam file object containing the short reads
   bamlist=list()
   src_files=list.files(pattern="*_sorted.bam$")
   for(filename in src_files){
     #Since we do not know which strand the reads were originally transcribed,
     #so set the strand to be ambiguous
     tmp=readBamGappedAlignments(filename)
     bamlist[[length(bamlist)+1]]=GRanges(seqnames=rname(tmp),
       ranges=IRanges(start=start(tmp),end=end(tmp)),
       strand=rep("*",length(tmp)))
   }
   names(bamlist)=src_files

# load the files into R, create an annotation object. 
# download one from the UCSC using the GenomicFeatures package
# here the ENSEMBL gene annotation is used

   library(GenomicFeatures)
   txdb=makeTranscriptDbFromUCSC(genome="hg19",tablename="ensGene")

# To compare genes for differential expression, we will summarize by gene 
# we choose to count all reads that fall within the body of the gene 
# (including introns) as counting towards a genes count

   tx_by_gene=transcriptsBy(txdb,"gene")

# Count the number of reads that fall in each gene for each lane 
# record the results in a table of counts.
#Initialize table of counts
   toc=data.frame(rep(NA,length(tx_by_gene)))
   for(i in 1:length(bamlist)){
     toc[,i]=countOverlaps(tx_by_gene,bamlist[[i]])
   }
 #Fix up labels
   rownames(toc)=names(tx_by_gene)
   colnames(toc)=names(bamlist)
  
# Normalization
# calculate appropriate scaling factors for normalization 
# via TMM method with the first lane as the reference.

library(edgeR)
norm_factors=calcNormFactors(as.matrix(toc))
   
# The counts themselves are not changed, 
# instead these scale factors are used as an offset in the negative binomial model 
# This is incorporated in the DGE list object required by edgeR.

 DGE=DGEList(toc,lib.size=norm_factors*colSums(toc),group=gsub("[0-9].*","",
 			colnames(toc)))

# Statistical Test:
# Calculate a common dispersion parameter represents extra Poisson variability 

    disp=estimateCommonDisp(DGE)

# Calculate p-values for genes being differentially expressed.

   tested=exactTest(disp)

   library(goseq)
#Apply benjamini hochberg correction for multiple testing
#choose genes with a p-value less than .05
   genes=as.integer(p.adjust(tested$table$p.value,method="BH") <.05)
   names(genes)=row.names(tested$table)

# Calculate the probability weighting function, 
# which quantifies the length bias effect.

   pwf=nullp(genes,"hg19","ensGene")

# Finally, we calculate the p-values for each GO category over represented 
# amongst DE genes.

   GO.pvals=goseq(pwf,"hg19","ensGene")

  library(DESeq)
 setwd("where the data is")
 c1<-read.table("counts_074284",row.names=1) 
 c2<-read.table("counts_074286",row.names=1)
 c3<-read.table("counts_074262",row.names=1)
 c4<-read.table("counts_074263",row.names=1)
 c5<-read.table("counts_074264",row.names=1)
 c6<-read.table("counts_074285",row.names=1)
 counts<-cbind(c1,c2,c3,c4,c5,c6)
 counts<-counts[-c(32679:32683),]  #remove the more general lines
 colnames(counts)<-c("P1","P2","P3","M1","M2","M3")
 design <- rep (c("P","Mo"),each=3)
 de  <-  newCountDataSet(counts, design)
 de  <-  estimateSizeFactors( de)
 # de  <-  estimateVarianceFunctions(  de  ) # This function has been removed. Use 'estimateDispersions' instead.  
 de <- estimateDispersions( de )
 res  <-  nbinomTest(  de,  "P",  "Mo")

# How many genes are significant?

 sum(na.omit(res$padj<0.05))
  
# Loading Urine1 Data
tx_by_gene=transcriptsBy(txdb,'gene')
geneListString<-c("CXCL12","TGFB","MYC","RAS","CXCR4","IL8","IL6")
 cuff<-readCufflinks(gtfFile=eval(refgtf),genome="genome.fa",rebuild=F)
 refgtf="cuffcmp.combined.gtf"
 genome_path="hg19"
 over="LUTS"
 under="CTRL"
-----
head(speciesMap)
txdb <- makeTxDbFromGFF(file=eval(refgtf), format="gff", organism="Hsapiens") #, available.species())
saveDb(txdb, file= "./data/tair10.sqlite")
## -----------------------------------
## ---- echo=FALSE---------------------------------------------------------


createAnnotationFile(gene_book)
write.Rsubread(gene_book)
genediff<-diffTable(isoforms(cuff))
gene_book = file.path(".GRCh38_gtf_only_Gene_Book.txt")
write.table(genediff,gene_book, sep="\t", row.names = T , col.names = T,quote = F)
genebook<-read.table(gene_book, header=T)
buildindex(basename="hg19_genome", reference="genome.fa")
propmapped(genebook,countFragments=TRUE,properlyPaired=FALSE)
align(index="hg19_genome",readfile1=genebooke,)
features_genes<-featureCounts(genediff, annot.inbuilt="hg19", annot.ext = NULL, 
              isGTFAnnotationFile=F, GTF.featureType="exon",
              GTF.attrType="gene_id", chrAliases=NULL, useMetaFeatures=TRUE,
              allowMultiOverlap=FALSE, minOverlap=1, largestOverlap=FALSE, 
              readExtension5=0, readExtension3=0)
## ---- eval=FALSE---------------------------------------------------------
hg19.tx <- makeTxDbFromUCSC(genome="hg19", tablename="refGene")

gAnnot <- exons(hg19.tx)

colnames(elementMetadata(gAnnot)) <- "exon"

gAnnot <- split(gAnnot,seqnames(gAnnot))

gene_strings<-writeXStringSet(Reduce(append,lapply(seqnames(Hsapiens),function(nam){
+ dss <- DNAStringSet(unmasked(Hsapiens[[nam]]))
+ names(dss)<-nam
+ dss})),file="cuffData.db")

library(DESeq2, quietly=TRUE); library(ape, warn.conflicts=FALSE)
countDF <- as.data.frame(genediff)
colData <- data.frame(row.names=targetsin(args)$SampleName, condition=targetsin(args)$Factor)
dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
d <- cor(assay(rlog(dds)), method="spearman")
hc <- hclust(dist(1-d))
pdf("results/sample_tree.pdf")
plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
dev.off()

system("wget ftp://mirbase.org/pub/mirbase/19/genomes/*.gff3 -P ./data/")
gff <- import.gff(refgtf)
gff <- split(gff, elementMetadata(gff)$ID)
bams <- names(bampaths); names(bams) <- targets$SampleName
bfl <- BamFileList(bams, yieldSize=50000, index=character())
countDFmiR <- summarizeOverlaps(gff, bfl, mode="Union", ignore.strand=FALSE, inter.feature=FALSE)
# Note: inter.feature=FALSE important since pre and mature miRNA ranges overlap
rpkmDFmiR <- apply(countDFmiR, 2, function(x) returnRPKM(counts=x, gffsub=gff))
write.table(assays(countDFmiR)$counts, "results/countDFmiR.xls", col.names=NA, quote=FALSE, sep="\t")
write.table(rpkmDFmiR, "results/rpkmDFmiR.xls", col.names=NA, quote=FALSE, sep="\t")


library(edgeR)
countDF <- read.delim("countDFeByg.xls", row.names=1, check.names=FALSE)
targets <- read.delim("targets.txt", comment="#")
cmp <- readComp(file="targets.txt", format="matrix", delim="-")
edgeDF <- run_edgeR(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=FALSE, mdsplot="")

# Add custom functional descriptions. Skip this step if desc.xls is not available.

desc <- read.delim("data/desc.xls")
desc <- desc[!duplicated(desc[,1]),]
descv <- as.character(desc[,2]); names(descv) <- as.character(desc[,1])
edgeDF <- data.frame(edgeDF, Desc=descv[rownames(edgeDF)], check.names=FALSE)
write.table(edgeDF, "./results/edgeRglm_allcomp.xls", quote=FALSE, sep="\t", col.names = NA)

#Filter and plot DEG results for up and down regulated genes. 
#The definition of ’up’ and ’down’ is given in the corresponding
# help file. To open it, type ?filterDEGs in the R console.

edgeDF <- read.delim("results/edgeRglm_allcomp.xls", row.names=1, check.names=FALSE)
pdf("results/DEGcounts.pdf")
DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=1))
write.table(DEG_list$Summary, "./results/DEGcounts.xls", quote=FALSE, sep="\t", row.names=FALSE)

vennsetup <- overLapper(DEG_list$Up[6:9], type="vennsets")
vennsetdown <- overLapper(DEG_list$Down[6:9], type="vennsets")
pdf("results/vennplot.pdf")
vennPlot(list(vennsetup, vennsetdown), mymain="", mysub="", colmode=2, ccol=c("blue", "red"))
dev.off()

# GO term enrichment analysis of DEGs
library("biomaRt")
listMarts() # To choose BioMart database
m <- useMart("ENSEMBL_MART_ENSEMBL"); listDatasets(m)
m <- (useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl"))
m <- useMart("ENSEMBL_MART_FUNCGEN"); listDatasets(m)
m <- useMart("ENSEMBL_MART_FUNCGEN", dataset="hsapiens_regulatory_feature" )
go <- getBM(values = genediff[,"gene_short_name"],attribute=c("go_id", "name_1006", "definition_1006"), mart=m)
go <- go[go[,3]!="",]; go[,3] <- as.character(go[,3])
head(go)
go[go[,3]=="molecular_function", 3] <- "F"; go[go[,3]=="biological_process", 3] <- "P"; go[go[,3]=="cellular_component"]
go[1:4,]
dir.create("./GRCh38_gtf_only_data/GO")
write.table(go, "GOannotationsBiomart_mod.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
catdb <- makeCATdb(myfile="GOannotationsBiomart_mod.txt", lib=NULL, org="", colno=c(1,2,3), idconv = T)

# Batch GO term enrichment analysis

BatchResult <- GOCluster_Report(catdb=catdb, method="all",setlist = go, id_type="gene", CLSZ=2)
library("biomaRt"); m <- useMart("ENSEMBL_MART_ENSEMBL", dataset="Hsapiens_gene_ensembl")
goslimvec <- as.character(getBM(attributes=c("goslim_goa_accession"), mart=m)[,1])
BatchResultslim <- GOCluster_Report(catdb=catdb, setlist=go, method="slim", id_type="gene", myslimv=go)

gos <- BatchResultslim[grep("M6-V6_up_down", BatchResultslim$CLID), ]
gos <- BatchResultslim
pdf("GOslimbarplotMF.pdf", height=8, width=10); goBarplot(gos, gocat="MF"); dev.off()
goBarplot(gos, gocat="BP")
goBarplot(gos, gocat="CC")

library(pheatmap)
geneids <- unique(as.character(unlist(genediff[["gene_name"]])))
y <- assay(genediff[geneids,])
pdf("heatmap1.pdf")
pheatmap(y, scale="row", clustering_distance_rows="correlation", clustering_distance_cols="correlation")

## ------------------------------------------------------------------------
hdr <- scanVcfHeader(file)
info(hdr) 
geno(hdr) 
meta(hdr)$META
## ------------------------------------------------------------------------

M<-gvisMotionChart(Fruits, "Fruit", "Year")
plot(M)
cat(M$html$chart, file = "tmp.html")
head(Fruits)

## ------------------------------------------------------------------------
## get entrez ids from gene symbols

genesym <- genediff["gene_short_name"]
geneid <- select(org.Hs.eg.db, keys=genesym, keytype="SYMBOL",
                 columns="ENTREZID")
geneid

## ------------------------------------------------------------------------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb

## ------------------------------------------------------------------------
txdb <- renameSeqlevels(txdb, gsub("chr", "", seqlevels(txdb)))
txdb <- keepSeqlevels(txdb,  )

## ------------------------------------------------------------------------
txbygene = transcriptsBy(txdb, "gene")

## ------------------------------------------------------------------------
gnrng <- unlist(range(txbygene[geneid$ENTREZID]), use.names=FALSE)
names(gnrng) <- geneid$SYMBOL

## ------------------------------------------------------------------------
param <- ScanVcfParam(which = gnrng, info = "DP", geno = c("GT", "cPd"))
param
 
## Extract the TRPV ranges from the VCF file 
vcf <- readVcf(file, "hg19", param)
## Inspect the VCF object with the 'fixed', 'info' and 'geno' accessors
vcf
 
head(fixed(vcf))

geno(vcf)

## ---- eval=FALSE---------------------------------------------------------
#  ## Use the 'region' argument to define the region
#  ## of interest. See ?locateVariants for details.
#  cds <- locateVariants(vcf, txdb, CodingVariants())
#  five <- locateVariants(vcf, txdb, FiveUTRVariants())
#  splice <- locateVariants(vcf, txdb, SpliceSiteVariants())
#  intron <- locateVariants(vcf, txdb, IntronVariants())

## ------------------------------------------------------------------------
all <- locateVariants(vcf, txdb, AllVariants())

## ------------------------------------------------------------------------
## Did any variants match more than one gene?
table(sapply(split(mcols(all)$GENEID, mcols(all)$QUERYID), 
      function(x) length(unique(x)) > 1))
 
## Summarize the number of variants by gene:
idx <- sapply(split(mcols(all)$QUERYID, mcols(all)$GENEID), unique)
sapply(idx, length)
 
## Summarize variant location by gene:
sapply(names(idx), 
    function(nm) {
        d <- all[mcols(all)$GENEID %in% nm, c("QUERYID", "LOCATION")]
        table(mcols(d)$LOCATION[duplicated(d) == FALSE])
    })

## ------------------------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg19)
seqlevelsStyle(vcf) <- "UCSC"
seqlevelsStyle(txdb) <- "UCSC"
aa <- predictCoding(vcf, txdb, Hsapiens)

## ------------------------------------------------------------------------
## Did any variants match more than one gene?
table(sapply(split(mcols(aa)$GENEID, mcols(aa)$QUERYID), 
        function(x) length(unique(x)) > 1))

## Summarize the number of variants by gene:
idx <- sapply(split(mcols(aa)$QUERYID, mcols(aa)$GENEID, drop=TRUE), unique)
sapply(idx, length)

## Summarize variant consequence by gene:
sapply(names(idx), 
       function(nm) {
           d <- aa[mcols(aa)$GENEID %in% nm, c("QUERYID","CONSEQUENCE")]
           table(mcols(d)$CONSEQUENCE[duplicated(d) == FALSE])
       })

## ---- echo=FALSE---------------------------------------------------------
suppressPackageStartupMessages(library(ensemblVEP))

## ------------------------------------------------------------------------
library(ensemblVEP)

## ------------------------------------------------------------------------
dest <- tempfile()
writeVcf(vcf, dest)

## ---- eval=FALSE---------------------------------------------------------
#  gr <- ensemblVEP(file = dest)
## ------------------------------------------------------------------------
VEPParam()

## ------------------------------------------------------------------------
basic(VEPParam())

## ----eval=FALSE----------------------------------------------------------
#  browseVignettes(package="VariantAnnotation")

## ----eval=FALSE----------------------------------------------------------
#  help.start()

## ------------------------------------------------------------------------
sessionInfo()

