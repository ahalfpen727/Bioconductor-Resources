## source("https://bioconductor.org/biocLite.R")
## biocLite("Rsubread")

## library(Rsubread)
## library(limma)
## library(edgeR)

## setwd("/Users/triley/cpct/rna-seq/smeliloti")

                                        # 20 samples

fc = featureCounts(

    files = "./createLinks_results_default/Sample_10211_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_10212_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_10213_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_10214_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_10215_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_1021B1_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_1021B2_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_1021B3_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_1021B4_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_1021B5_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_A1_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_A2_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_A3_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_A4_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_A5_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_AB1_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_AB2_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_AB3_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_AB4_R1.fastq.gz.subread.BAM;./createLinks_results_default/Sample_AB5_R1.fastq.gz.subread.BAM",

    ##                                     # BAM/SAM files
    ## files = "./createLinks_results_default/Sample_10212_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_10214_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_10215_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_1021B2_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_1021B3_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_1021B4_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_1021B5_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_A2_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_A4_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_A5_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_AB1_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_AB2_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_AB3_R1.fastq.gz.subread.BAM;
    ##          ./createLinks_results_default/Sample_AB5_R1.fastq.gz.subread.BAM",

    ## files = "./tophat_results_Smeliloti_1021_unguided/Sample_10212_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_10214_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_10215_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_1021B2_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_1021B3_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_1021B4_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_1021B5_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_A2_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_A4_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_A5_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_AB1_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_AB2_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_AB3_out/accepted_hits.bam;
    ##          ./tophat_results_Smeliloti_1021_unguided/Sample_AB5_out/accepted_hits.bam",

    ## files = c(
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_10212_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_10214_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_10215_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_1021B2_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_1021B3_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_1021B4_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_1021B5_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_A2_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_A4_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_A5_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_AB1_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_AB2_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_AB3_out/accepted_hits.bam,
    ##     ./tophat_results_Smeliloti_1021_unguided/Sample_AB5_out/accepted_hits.bam
    ## ),
                                        # annotation
    annot.inbuilt=NULL,
    annot.ext="Smeliloti_1021_refseq_bowtie_assembly/1021_genome.gtf",
    isGTFAnnotationFile=TRUE,
    GTF.featureType="CDS",
    GTF.attrType="Parent",
    chrAliases=NULL,
                                        # level of summarization
    useMetaFeatures=TRUE,
                                        # overlap between reads and features
    allowMultiOverlap=FALSE,
    minOverlap=1,
    ## fracOverlap=0,
    largestOverlap=FALSE,
    readExtension5=0,
    readExtension3=0,
    read2pos=NULL,
                                        # multi-mapping reads
    countMultiMappingReads=FALSE,
                                        # fractional counting
    fraction=FALSE,
                                        # read filtering
    minMQS=0,
    splitOnly=FALSE,
    nonSplitOnly=FALSE,
    primaryOnly=FALSE,
    ignoreDup=FALSE,
                                        # strandness
    strandSpecific=0,
                                        # exon-exon junctions
    juncCounts=FALSE,
    genome=NULL,
                                        # parameters specific to paired end reads
    isPairedEnd=TRUE,
    requireBothEndsMapped=TRUE,
    checkFragLength=FALSE,
    minFragLength=50,
    maxFragLength=1600,
    countChimericFragments=TRUE,
    autosort=TRUE,
                                        # number of CPU threads
    nthreads=1,
                                        # miscellaneous
    maxMOp=10,
    ## tmpDir=".",
    reportReads=FALSE

)

saveRDS(fc, "featureCounts.list")
## # fc = readRDS("featureCounts.list")

## fc$counts[1:5,]
## ## A_1.bam A_2.bam B_1.bam B_2.bam
## ## 653635        642     522     591     596
## ## 100422834       1       0       0       0
## ## 645520          5       3       0       0
## ## 79501           0       0       0       0
## ## 729737         82      72      30      25

## fc$annotation[1:5,c("GeneID","Length")]
## ## GeneID Length
## ## 1    653635   1769
## ## 2 100422834    138
## ## 3    645520   1130
## ## 4     79501    918
## ## 5    729737   3402

## ## Create a DGEList object.
## x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])

## #Calculate RPKM (reads per kilobases of exon per million reads mapped) values for genes:
## x_rpkm <- rpkm(x,x$genes$Length,prior.count=0)
## x_rpkm[1:5,]
## ## A_1.bam A_2.bam B_1.bam B_2.bam
## ## 653635        939   905.0     709     736
## ## 47
## ## 100422834      19     0.0       0       0
## ## 645520         11     8.1       0       0
## ## 79501           0     0.0       0       0
## ## 729737         62    64.9      19      16

## ## Filtering. Only keep in the analysis those genes which had
## ## 10 reads per million mapped reads in at least two libraries.
## ## isexpr <- rowSums(cpm(x) > 10) >= 2
## isexpr <- rowSums(cpm(x) > 10) >= 5
## x <- x[isexpr,]

## ## Create a design matrix:
## fc$group = c('wt1021','wt1021','wt1021','wt1021','wt1021','wt1021B','wt1021B','wt1021B','wt1021B','wt1021B','A','A','A','A','A','AB','AB','AB','AB','AB')

## groupsFactor <- factor(fc$group)
## design <- model.matrix(~0+groupsFactor)
## colnames(design) <- levels(groupsFactor)
## ## celltype <- factor(targets$CellType)
## ## design <- model.matrix(~0+celltype)
## ## colnames(design) <- levels(celltype)

## ## Normalization. Perform voom normalization:
## ## y <- voom(x,design,plot=FALSE)
## y <- voom(x,design,plot=TRUE)

## ## Multi-dimensional scaling (MDS) plot
## plotMDS(y,xlim=c(-2.5,2.5))

## ## Linear model fitting and differential expression analysis.
## ## Fit linear models to genes and assess differential expression using eBayes moderated t statistic.  Here we compare sample
## ## B vs sample A.
## fit <- lmFit(y,design)
## options(digits=3)








## contr <- makeContrasts(wt1021Bvswt1021=wt1021B - wt1021,levels=design)
## # write.fit(fit.contr, file="wt1021B-over-wt1021.txt")
## fit.contr <- eBayes(contrasts.fit(fit,contr))
## dt <- decideTests(fit.contr)
## summary(dt)

## ## List top 10 differentially expressed genes:
## topTable(fit.contr)

## ## topTable(fit.contr, sort="none",n=Inf)
## topOutput = topTable(fit.contr, n=Inf)
## topOutput$FC = 2^(topOutput$logFC)
## write.table(topOutput, file="wt1021B-over-wt1021.txt", sep="\t", quote=FALSE)








## contr <- makeContrasts(Avswt1021=A-wt1021,levels=design)
## # write.fit(fit.contr, file="A-over-wt1021.txt")
## fit.contr <- eBayes(contrasts.fit(fit,contr))
## dt <- decideTests(fit.contr)
## summary(dt)

## ## List top 10 differentially expressed genes:
## topTable(fit.contr)

## ## topTable(fit.contr, sort="none",n=Inf)
## topOutput = topTable(fit.contr, n=Inf)
## topOutput$FC = 2^(topOutput$logFC)
## write.table(topOutput, file="A-over-wt1021.txt", sep="\t", quote=FALSE)









## contr <- makeContrasts(ABvswt1021=AB-wt1021,levels=design)
## # write.fit(fit.contr, file="AB-over-wt1021.txt")
## fit.contr <- eBayes(contrasts.fit(fit,contr))
## dt <- decideTests(fit.contr)
## summary(dt)

## ## List top 10 differentially expressed genes:
## topTable(fit.contr)

## ## topTable(fit.contr, sort="none",n=Inf)
## topOutput = topTable(fit.contr, n=Inf)
## topOutput$FC = 2^(topOutput$logFC)
## write.table(topOutput, file="AB-over-wt1021.txt", sep="\t", quote=FALSE)









## contr <- makeContrasts(AvsAB=A-AB,levels=design)
## # write.fit(fit.contr, file="A-over-AB.txt")
## fit.contr <- eBayes(contrasts.fit(fit,contr))
## dt <- decideTests(fit.contr)
## summary(dt)

## ## List top 10 differentially expressed genes:
## topTable(fit.contr)

## ## topTable(fit.contr, sort="none",n=Inf)
## topOutput = topTable(fit.contr, n=Inf)
## topOutput$FC = 2^(topOutput$logFC)
## write.table(topOutput, file="A-over-AB.txt", sep="\t", quote=FALSE)








## contr <- makeContrasts(ABvswt1021B=AB-wt1021B,levels=design)
## # write.fit(fit.contr, file="AB-over-wt1021B.txt")
## fit.contr <- eBayes(contrasts.fit(fit,contr))
## dt <- decideTests(fit.contr)
## summary(dt)

## ## List top 10 differentially expressed genes:
## topTable(fit.contr)

## ## topTable(fit.contr, sort="none",n=Inf)
## topOutput = topTable(fit.contr, n=Inf)
## topOutput$FC = 2^(topOutput$logFC)
## write.table(topOutput, file="AB-over-wt1021B.txt", sep="\t", quote=FALSE)









## GeneID Length logFC AveExpr   t P.Value adj.P.Val  B
## 100131754 100131754   1019   1.6      16 113 3.5e-28   6.3e-25 54
## 2023           2023   1812  -2.7      13 -91 2.2e-26   1.9e-23 51
## 2752           2752   4950   2.4      13  82 1.5e-25   9.1e-23 49
## 22883         22883   5192   2.3      12  64 1.8e-23   7.9e-21 44
## 6135           6135    609  -2.2      12 -62 3.1e-23   9.5e-21 44
## 6202           6202    705  -2.4      12 -62 3.2e-23   9.5e-21 44
