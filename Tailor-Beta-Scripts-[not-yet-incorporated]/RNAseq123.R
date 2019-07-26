### R code from vignette source 'goseq.Rnw'
library(RNAseq123)
library(goseq)
filter <- apply(g.cnt.matrix,1,function(x) mean(x)>10)
table(filter)

#gene.vector=as.integer(assayed.genes%in%de.genes)
#names(gene.vector)=assayed.genes
#head(gene.vector)

###################################################
### code chunk number 5: getLengthDataFromUCSC (eval = FALSE)
###################################################
supportedOrganisms()
# txsByGene=transcriptsBy(txdb,"gene")
# lengthData=median(width(txsByGene))


###################################################
### code chunk number 6: edger_1
###################################################
library(edgeR)
table.summary=read.table(system.file("extdata","Li_sum.txt",package='goseq'),
                         sep='\t',header=TRUE,stringsAsFactors=FALSE)
counts=table.summary[,-1]
rownames(counts)=table.summary[,1]


###################################################
### code chunk number 7: edger_2
###################################################
grp <- rep(c("LUTS", "CTRL"), each=9)
summarized=DGEList(G.rep.cnt.matrix,lib.size=colSums(G.rep.cnt.matrix),group=grp)

disp=estimateCommonDisp(summarized)
disp$common.dispersion
tested=exactTest(disp)
topTags(tested)


###################################################
### code chunk number 8: edger_3
###################################################
genes=as.integer(p.adjust(tested$table$PValue[tested$table$logFC!=0],
                          method="BH")<.05)

names(genes)=row.names(tested$table[tested$table$logFC!=0,])
table(genes)

supportedOrganisms()[supportedOrganisms()$Genome=="hg19",]
pwf=nullp(genes,"hg19","geneSymbol")
#pwf=nullp(genes,"hg19","ensGene")
head(pwf)

GO.wall=goseq(pwf,"hg19","geneSymbol")
head(GO.wall)

GO.samp=goseq(pwf,"hg19","geneSymbol",method="Sampling",repcnt=1000)
head(GO.samp)

plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.wall[,1],GO.samp[,1]),2]),
     xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",
     xlim=c(-3,0))
abline(0,1,col=3,lty=2)


###################################################
### code chunk number 16: GO.nobias
###################################################
GO.nobias=goseq(pwf,"hg19","geneSymbol",method="Hypergeometric")
head(GO.nobias)

plot(log10(GO.wall[,2]), log10(GO.nobias[match(GO.wall[,1],GO.nobias[,1]),2]),
     xlab="log10(Wallenius p-values)", ylab="log10(Hypergeometric p-values)",
     xlim=c(-3,0), ylim=c(-3,0))
abline(0,1,col=3,lty=2)


###################################################
### code chunk number 18: GO.limited
###################################################
GO.MF=goseq(pwf,"hg19","geneSymbol",test.cats=c("GO:MF"))
head(GO.MF)

enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue,
                                      method="BH")<.05]
head(enriched.GO)


###################################################
### code chunk number 20: GO_explained
###################################################
library(GO.db)
for(go in enriched.GO[1:10]){
   print(GOTERM[[go]])
   cat("--------------------------------------\n")
}


###################################################
### code chunk number 21: getlength
###################################################
len=getlength(names(genes),"hg19","geneSymbol")
length(len)
length(genes)
head(len)


###################################################
### code chunk number 22: getgo
###################################################
go=getgo(names(genes),"hg19","geneSymbol")
length(go)
length(genes)
head(go)


###################################################
### code chunk number 23: conv_table
###################################################
goseq:::.ID_MAP
goseq:::.ORG_PACKAGES


###################################################
### code chunk number 24: norm_analysis (eval = FALSE)
###################################################
## pwf=nullp(genes,"hg19","ensGene")
## go=goseq(pwf,"hg19","ensGene")


###################################################
### code chunk number 25: verbose_analysis (eval = FALSE)
###################################################
## gene_lengths=getlength(names(genes),"hg19","ensGene")
## pwf=nullp(genes,bias.data=gene_lengths)
## go_map=getgo(names(genes),"hg19","ensGene")
## go=goseq(pwf,"hg19","ensGene",gene2cat=go_map)


###################################################
### code chunk number 26: KEGG_mappings (eval = FALSE)
###################################################
## # Get the mapping from ENSEMBL 2 Entrez
## en2eg=as.list(org.Hs.egENSEMBL2EG)
## # Get the mapping from Entrez 2 KEGG
## eg2kegg=as.list(org.Hs.egPATH)
## # Define a function which gets all unique KEGG IDs
## # associated with a set of Entrez IDs
## grepKEGG=function(id,mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))}
## # Apply this function to every entry in the mapping from
## # ENSEMBL 2 Entrez to combine the two maps
## kegg=lapply(en2eg,grepKEGG,eg2kegg)
## head(kegg)


###################################################
### code chunk number 27: KEGG (eval = FALSE)
###################################################
## pwf=nullp(genes,"hg19","ensGene")
## KEGG=goseq(pwf,gene2cat=kegg)
## head(KEGG)


###################################################
### code chunk number 28: KEGG_goseq
###################################################
pwf=nullp(genes,'hg19','ensGene')
KEGG=goseq(pwf,'hg19','ensGene',test.cats="KEGG")
head(KEGG)


###################################################
### code chunk number 29: KEGG_from_db
###################################################
kegg=as.list(org.Hs.egPATH)
head(kegg)


###################################################
### code chunk number 30: countbias
###################################################
countbias=rowSums(counts)[rowSums(counts)!=0]
length(countbias)
length(genes)


###################################################
### code chunk number 31: GO.counts
###################################################
pwf.counts=nullp(genes,bias.data=countbias)
GO.counts=goseq(pwf.counts,"hg19","ensGene")
head(GO.counts)


###################################################
### code chunk number 32: setup
###################################################
sessionInfo()


