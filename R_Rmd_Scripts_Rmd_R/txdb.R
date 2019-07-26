source("https://bioconductor.org/biocLite.R")
biocLite("TxDb.Hsapiens.UCSC.hg18.knownGene")
library(TxDb.Hsapiens.UCSC.hg18.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene
class(txdb) ## do some digging around!
library(org.Hs.eg.db)
>  cols(org.Hs.eg.db)
select(org.Hs.eg.db, keys="APOE", cols=c("ENTREZID", "SYMBOL", "GENENAME"), keytype="SYMBOL")
tx.by.gene["348"]
# CAREFUL: use levels() to check that you're making new factor names
# that correspond to the old ones!
levels(d$chrom) <- paste("chr", c(1:22, "X", "Y", "M"), sep="")
my.snps <- with(d, GRanges(seqnames=chrom,
                           IRanges(start=position, width=1),
                           rsid=rsid, genotype=genotype)) # this goes into metadata
apoe.i <- findOverlaps(tx.by.gene["348"], my.snps)
hits <- matchMatrix(apoe.i)[, "subject"]
my.snps[hits]
gwrngs.emd <- as.data.frame(elementMetadata(gwrngs))
dm <- merge(d, gwrngs.emd, by.x="rsid", by.y="SNPs")

We can search for the risk allele in the 23andme genotype data with R and attach a vector of i.have.risk to the dm data frame:

  risk.alleles <- gsub("[^\\-]*-([ATCG?])", "\\1", dm$Strongest.SNP.Risk.Allele)
i.have.risk <- mapply(function(risk, mine) {
  risk %in% unlist(strsplit(mine, ""))
}, risk.alleles, dm$genotype)
dm$i.have.risk <- i.have.risk

Now that you have this data frame, you can mine it endlessly. You may want to sort by Risk.Allele.Frequency and whether you have the risk. Because there are quite a few columns in the element metadata, it’s nice to define a quick-summary subset:

  > my.risk <- dm[dm$i.have.risk, ]
> rel.cols <- c(colnames(d), "Disease.Trait", "Risk.Allele.Frequency",
                "p.Value", "i.have.risk", "X95..CI..text.")

> head(my.risk[order(my.risk$Risk.Allele.Frequency), rel.cols], 1)
rsid chrom position genotype Disease.Trait Risk.Allele.Frequency
2553 rs2315504 chr17 36300407       AC        Height                  0.01
p.Value i.have.risk   X95..CI..text.
2553   8e-06        TRUE [NR] cm increase

This is a rare variant, but the most important next question is, rare in who?

  > dm[which(dm$rsid == "rs2315504"), "Initial.Sample.Size"]
[1] 8,842 Korean individuals

So this clearly doesn’t mean much to me. We can use grep to find studies that mention “European”:

  > head(my.risk[grep("European", my.risk$Initial.Sample.Size), rel.cols], 30)


library(ggbio)
> p <- plotOverview(hg19IdeogramCyto, cytoband=FALSE)

Now, let’s take the gwrngs object and subset by my risk alleles. Notice how these assignment function elementMetadata<- is overloaded here:

  (elementMetadata(gwrngs)$my.genotype <-
     d$genotype[(match(elementMetadata(gwrngs)$SNPs, d$rsid))])

elementMetadata(gwrngs)$my.risk <- with(elementMetadata(gwrngs),
                                        mapply(function(risk, mine) {
                                          risk %in% unlist(strsplit(mine, ""))
                                        }, gsub("[^\\-]*-([ATCG?])", "\\1", Strongest.SNP.Risk.Allele), my.genotype))

Now to plot these regions:

  p + geom_hotregion(gwrngs, aes(color=my.risk))
















