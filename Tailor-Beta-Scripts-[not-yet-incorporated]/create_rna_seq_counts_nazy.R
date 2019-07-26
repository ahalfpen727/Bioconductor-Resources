# source("http://bioconductor.org/biocLite.R")

rpackage.dir = "/project/umb_triley/Rpackages/"

# biocLite("Rsamtools", lib.loc=rpackage.dir)
library(Rsamtools, lib.loc= rpackage.dir)

dir = "/project/umb_cpct/data/nazy/hiseq/12.01.2014/141201-A.align/Project_I00005/"

f = read.table("/project/umb_cpct/data/nazy/hiseq/12.01.2014/bam_files_nazy.txt")

filenames = paste(dir, f[,1], ".bam", sep="")
bamfiles = BamFileList(filenames)

# # biocLite("GenomicFeatures", lib.loc= rpackage.dir)
# biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
genes = exonsBy(txdb, by="gene")

# biocLite("GenomicAlignments", lib.loc=rpackage.dir)
library(GenomicAlignments, lib.loc = rpackage.dir)
se = summarizeOverlaps(features=genes, reads=bamfiles,
                       mode="Union",
                       singleEnd=FALSE,
                       ignore.strand=TRUE,
                       fragments=TRUE )

cts = assay(se)
write.table(cts, "analysis/counts_nazy.txt")
print("Done.")