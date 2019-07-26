# load libraries
library(Rsubread);library(limma);library(edgeR);library(DESeq2)
library(DESeq2, quietly=TRUE); library(ape,  warn.conflicts=FALSE)

A.over.AB.genes<-read.delim(file="A-over-AB/genes.read_group_tracking" , as.is = T,stringsAsFactors = T, skip=0)
samplesA_AB<-paste(A.over.AB.genes$condition, A.over.AB.genes$replicate, sep ="_")
condition.A.AB<-as.factor(A.over.AB.genes$condition)
genesA.AB<-A.over.AB.genes$tracking_id
length(genesA.AB)

wt1021B.over.wt1021.genes<-read.delim(file="1021B-over-1021/genes.read_group_tracking",stringsAsFactors = T,  as.is=T, skip=0)
samples1021B_1021<-paste(wt1021B.over.wt1021.genes$condition, wt1021B.over.wt1021.genes$replicate, sep ="_")
condition.1021B.1021<-wt1021B.over.wt1021.genes$condition
genes1021B.1021<-wt1021B.over.wt1021.genes$tracking_id
length(genes1021B.1021)

A.AB.gene.exp<-as.data.frame(cbind(genesA.AB,A.over.AB.genes[,c(2,7)]))
A.AB.Gene.Exp<-cbind(samplesA_AB, A.AB.gene.exp)
head(A.AB.Gene.Exp)

A.AB.gene.exp<-as.data.frame(cbind(genesA.AB,A.over.AB.genes[,c(2,7)]))
A.AB.Gene.Exp<-cbind(samplesA_AB, A.AB.gene.exp)
head(A.AB.Gene.Exp)

A.gene.exp<-A.AB.Gene.Exp[which(A.AB.Gene.Exp$condition ==  'A'),]
head(A.gene.exp)

A0.gene.exp<-A.AB.Gene.Exp[which(A.AB.Gene.Exp$samplesA_AB ==  'A_0'),]
A0.GeneExp<-as.data.frame(A0.gene.exp[,c(-1,-3)])
colnames(A0.GeneExp)<-c('gene','A_0')
head(A0.GeneExp)

A1.gene.exp<-A.AB.Gene.Exp[which(A.AB.Gene.Exp$samplesA_AB ==  'A_1'),]
A1.GeneExp<-as.data.frame(A1.gene.exp[,c(-1,-3)])
colnames(A1.GeneExp)<-c('gene','A_1')
head(A1.GeneExp)

A2.gene.exp<-A.AB.Gene.Exp[which(A.AB.Gene.Exp$samplesA_AB ==  'A_2'),]
A2.GeneExp<-as.data.frame(A2.gene.exp[,c(-1,-3)])
colnames(A2.GeneExp)<-c('gene','A_2')
head(A2.GeneExp)

A3.gene.exp<-A.AB.Gene.Exp[which(A.AB.Gene.Exp$samplesA_AB ==  'A_3'),]
A3.GeneExp<-as.data.frame(A3.gene.exp[,c(-1,-3)])
colnames(A3.GeneExp)<-c('gene','A_3')
head(A3.GeneExp)

A4.gene.exp<-A.AB.Gene.Exp[which(A.AB.Gene.Exp$samplesA_AB ==  'A_4'),]
A4.GeneExp<-as.data.frame(A4.gene.exp[,c(-1,-3)])
colnames(A4.GeneExp)<-c('gene','A_4')
head(A4.GeneExp)

A01<-cbind(A0.GeneExp, A1.GeneExp$A_1) 
A012<-cbind(A01, A2.GeneExp$A_2) 
A0123<-cbind(A012, A3.GeneExp$A_3) 
A01234.gene.exp<-cbind(A0123, A4.GeneExp$A_4) 
head(A01234.gene.exp)

AB.gene.exp<-A.AB.Gene.Exp[which(A.AB.Gene.Exp$condition ==  'AB'),]
head(AB.gene.exp)

AB0.gene.exp<-A.AB.Gene.Exp[which(A.AB.Gene.Exp$samplesA_AB ==  'AB_0'),]
AB0.GeneExp<-as.data.frame(AB0.gene.exp[,c(-1,-3)])
colnames(AB0.GeneExp)<-c('gene','AB_0')
head(AB0.GeneExp)

AB1.gene.exp<-A.AB.Gene.Exp[which(A.AB.Gene.Exp$samplesA_AB ==  'AB_1'),]
AB1.GeneExp<-as.data.frame(AB1.gene.exp[,c(-1,-3)])
colnames(AB1.GeneExp)<-c('gene','AB_1')
head(AB1.GeneExp)

AB2.gene.exp<-A.AB.Gene.Exp[which(A.AB.Gene.Exp$samplesA_AB ==  'AB_2'),]
AB2.GeneExp<-as.data.frame(AB2.gene.exp[,c(-1,-3)])
colnames(AB2.GeneExp)<-c('gene','AB_2')
head(AB2.GeneExp)

AB3.gene.exp<-A.AB.Gene.Exp[which(A.AB.Gene.Exp$samplesA_AB ==  'AB_3'),]
AB3.GeneExp<-as.data.frame(A3.gene.exp[,c(-1,-3)])
colnames(AB3.GeneExp)<-c('gene','AB_3')
head(AB3.GeneExp)

AB4.gene.exp<-A.AB.Gene.Exp[which(A.AB.Gene.Exp$samplesA_AB ==  'AB_4'),]
AB4.GeneExp<-as.data.frame(A4.gene.exp[,c(-1,-3)])
colnames(AB4.GeneExp)<-c('gene','AB_4')
head(AB4.GeneExp)

AB_0<-AB0.GeneExp
AB_1<-AB1.GeneExp$AB_1
AB_2<-AB2.GeneExp$AB_2
AB_3<-AB3.GeneExp$AB_3
AB_4<-AB4.GeneExp$AB_4

AB01<-cbind(AB_0, AB_1) 
AB012<-cbind(AB01, AB_2) 
AB0123<-cbind(AB012, AB_3) 
AB01234<-cbind(AB0123, AB_4)
head(AB01234)

AB0.gene.exp<-AB.gene.exp[which(AB.gene.exp$samplesA_AB ==  'AB_0'),]
AB1.gene.exp<-AB.gene.exp[which(AB.gene.exp$samplesA_AB ==  'AB_1'),]
AB2.gene.exp<-AB.gene.exp[which(AB.gene.exp$samplesA_AB ==  'AB_2'),]
AB3.gene.exp<-AB.gene.exp[which(AB.gene.exp$samplesA_AB ==  'AB_3'),]
AB4.gene.exp<-AB.gene.exp[which(AB.gene.exp$samplesA_AB ==  'AB_4'),]

wt1021B.wt1021.gene.exp<-as.data.frame(cbind(genes1021B.1021, wt1021B.over.wt1021.genes[,c(2,7)]))
wt1021B.wt1021.Gene.Exp<-cbind(samples1021B_1021, wt1021B.wt1021.gene.exp)

wt1021B.Gene.Exp<-wt1021B.wt1021.Gene.Exp[which(wt1021B.wt1021.Gene.Exp$condition == '1021B'),]
head(wt1021B.Gene.Exp)

wt1021B_0.gene.exp<-wt1021B.wt1021.Gene.Exp[which(wt1021B.wt1021.Gene.Exp$genes1021B.1021 ==  '1021B_0'),]
wt1021B_1.gene.exp<-wt1021B.wt1021.Gene.Exp[which(wt1021B.wt1021.Gene.Exp$genes1021B.1021 ==  '1021B_1'),]
wt1021B_2.gene.exp<-wt1021B.wt1021.Gene.Exp[which(wt1021B.wt1021.Gene.Exp$genes1021B.1021 ==  '1021B_2'),]
wt1021B_3.gene.exp<-wt1021B.wt1021.Gene.Exp[which(wt1021B.wt1021.Gene.Exp$genes1021B.1021 ==  '1021B_3'),]
wt1021B_4.gene.exp<-wt1021B.wt1021.Gene.Exp[which(wt1021B.wt1021.Gene.Exp$genes1021B.1021 ==  '1021B_4'),]

wt1021B_0<-wt1021B_0.gene.exp
wt1021B_1<-wt1021B_1.gene.exp["1021B_1",]
wt1021B_2<-wt1021B_2.gene.exp["1021B_2",]
wt1021B_3<-wt1021B_3.gene.exp["1021B_3",]
wt1021B_4<-wt1021B_4.gene.exp["1021B_4",]

wt1021B_01<-cbind(wt1021B_0, wt1021B_1) 
wt1021B_012<-cbind(wt1021B_01, wt1021B_2[(-1),]) 
wt1021B_0123<-cbind(wt1021B_012, wt1021B_3[(-1),]) 
wt1021B_01234<-cbind(wt1021B_0123, wt1021B_4[(-1),])
head(wt1021B_01234)

wt1021.gene.exp<-wt1021B.wt1021.Gene.Exp[which(wt1021B.wt1021.Gene.Exp$condition == '1021'),]
head(wt1021.gene.exp)
wt1021B.wt1021.Gene.Exp<-wt1021.gene.exp[which(wt1021.gene.exp$condition ==  '1021'),]

wt1021_0.gene.exp<-wt1021B.wt1021.Gene.Exp[which(wt1021B.wt1021.Gene.Exp$samples1021B_1021 ==  '1021_0'),]
wt1021_1.gene.exp<-wt1021B.wt1021.Gene.Exp[which(wt1021B.wt1021.Gene.Exp$samples1021B_1021 ==  '1021_1'),]
wt1021_2.gene.exp<-wt1021B.wt1021.Gene.Exp[which(wt1021B.wt1021.Gene.Exp$samples1021B_1021 ==  '1021_2'),]
wt1021_3.gene.exp<-wt1021B.wt1021.Gene.Exp[which(wt1021B.wt1021.Gene.Exp$samples1021B_1021 ==  '1021_3'),]
wt1021_4.gene.exp<-wt1021B.wt1021.Gene.Exp[which(wt1021B.wt1021.Gene.Exp$samples1021B_1021 ==  '1021_4'),]

GeneExp<-merge(wt1021B.wt1021.Gene.Exp, A.AB.Gene.Exp,by.x="genes1021B.1021" ,by.y="genesA.AB")
Gene.Exp<-rbind(wt1021B.wt1021.Gene.Exp, A.AB.Gene.Exp) #,by.x="genes1021B.1021" ,by.y="genesA.AB")
dim(GeneExp);head(GeneExp)
ExpGene<-t(GeneExp)
dim(ExpGene)

GeneNames.1021B.1021<-unique(GeneExp.by.Samp.1021B.1021$tracking_id)
GeneExpr.Sme1021B.Sme1021<-as.data.frame(GeneExp.by.Samp.1021B.1021$FPKM, col.names=samples1021B_1021, rownames = GeneNames.1021B.1021)
dim(GeneExpr.Sme1021B.Sme1021)

GeneExpr.wt1021B.1021<-GeneExp.by.Samp.1021B.1021[,-1]
Gene.Expr.wt1021B.1021<-t(GeneExpr.wt1021B.1021)
rownames(Gene.Expr.wt1021B.1021)

allgene.expr<-merge(AoverAB.gene.expr, wt1021Bover1021.gene.expr, by = intersect(wt1021Bover1021.gene.expr$tracking_id,AoverAB.gene.expr$tracking_id),
      sort = TRUE, suffixes = c(".A.AB",".1021B.1021"))
ncol(AoverAB.genes)

dim(wt1021B.over.wt1021.genes)
head(wt1021B.over.wt1021.genes)
ncol(wt1021B.over.wt1021.genes)

all.samples<-factor(unique(cbind(samplesA_AB, samples1021B_1021)))
length(all.samples)

FPKM.ma<-cbind(condtion.samps,wt1021Bover1021B.genes)
head(FPKM.ma)

sampleXgene<-as.matrix(FPKM.ma, col.names=FPKM.ma$condtion.samps, row.names=FPKM.ma$tracking_id)
dim(sampleXgene)
geneEXPRSma<-sampleXgene[,c(1,8)]
names(sampleXgene)<-unique(condition.rep)
nrow(t(sampleXgene))
gene.by.samp<-(sampleXgene)
dim(gene.by.samp);head(gene.by.samp)
condition.rep<-as.factor(x=condtion.samps)

wt1021B_over_1021genes<-read.delim(file="genes.read_group_tracking")
wt1021b-over1021.genes<-read.table(file="genes.read_group_tracking")
wt1021B_over_1021genes<-read.delim(file="genes.read_group_tracking")

wt1021b-over1021.genes<-read.table(file="genes.read_group_tracking")
wt1021B_over_1021genes<-read.delim(file="genes.read_group_tracking")
head(wt1021B_over_1021genes)
colnames(wt1021B_over_1021genes)

samples<-list()
condition_type<-wt1021B_over_1021genes$condition
replicate_num<-wt1021B_over_1021genes$replicate
stripchart(FPKM ~ replicate_num, wt1021B_over_1021genes)


condition1021_1021B<-factor(wt1021B_over_1021genes$condition,condition<-c("1021","1021B"),ordered = T)
head(condition1021_1021B)
replicate1021_1021B<-factor(wt1021B_over_1021genes$replicate,ordered = T)
max(replicate1021_1021B)
sample.list<-list()
for (i in wt1021B_over_1021genes) {
  sample.condition.num<-merge(wt1021B_over_1021genes[i,"condition"],wt1021B_over_1021genes[i,"replicate"], sep="_")
  print(sample.condition.num)
  append( sample.condition.num, sample.list)}
sample.list
Samp
#####################################################################
# read in target file
options(digits=2)
targets <- readTargets("cuffData.db")
targets

# create a design matrix
celltype <- factor(targets$CellType)
design <- model.matrix(~celltype)

# build an index for reference sequence (Chr1 in hg19)
#buildindex(basename="chr1",reference="hg19_chr1.fa")
# align reads
#align(index="chr1",readfile1=targets$InputFile,input_format="gzFASTQ",output_format="BAM",output_file=targets$OutputFile,unique=TRUE,indels=5)
geneCTS<-count(genes(cuff))
# count numbers of reads mapped to NCBI Refseq genes
fc <- featureCounts(files=targets$geneCTS,annot.inbuilt="hg19")
x <- DGEList(counts=geneCTS[,"count"], genes=geneCTS$annotation[,c("gene_id")])

# filter out low-count genes
isexpr <- rowSums(cpm(x) > 10) >= 2
x <- x[isexpr,]

# perform voom normalization
y <- voom(x, plot=TRUE)

# cluster libraries
plotMDS(x,xlim=c(-2.5,2.5))

# fit linear model and assess differential expression
fit <- eBayes(lmFit(y,design))
topTable(fit,coef=2)

######################################################################################################
## Extracting Gene Expression Matrix from cuffdiff output
######################################################################################################
test<-Sys.geten("PATH")
print(paste("test = ", test, sep=""))


