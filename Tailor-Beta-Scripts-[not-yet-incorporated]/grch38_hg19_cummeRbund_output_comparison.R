library(multtest);library(outliers)
data(golub,package="multtest")

gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
 sh <- apply(golub[,gol.fac=="ALL"], 1, function(x) shapiro.test(x)$p.value)
> sum(sh > 0.05)/nrow(golub) * 100
grch38anov<-anova(lm(rowMeans(grch38.lutsgrp) ~ rowMeans(grch38.ctrlgrp)))
grubbs.test(golub[1042, gol.fac=="ALL"])

grch38.default.cuffdifftable<-file.path("cummeRbund_results_grch38_default/LUTS-over-CTRL/fpkmMatrix.csv")
grch38_default.cufffeatures<-read.table(grch38.default.cuffdifftable,header=T,as.is=T,stringsAsFactors = T,check.names = T)
head(grch38_default.cufffeatures)
dim(grch38_default.cufffeatures)

grch38.expr <- rowSums(grch38_default.cufffeatures > 0) >= 5
head(grch38.expr)
grch38.default_expr <- grch38_default.cufffeatures[grch38.expr,]
head(grch38.default_expr)
dim(grch38.default_expr)

grch38.default_exprs<-grch38.default_expr[,-1]
head(grch38.default_exprs)
dim(grch38.default_exprs)

############################################################################################################
# pairwise cor between samples
############################################################################################################
cor(grch38.default_exprs[,c(1:4)],grch38.default_exprs[,c(5:9)])
cor(grch38.default_exprs[,10:14],grch38.default_exprs[,15:18])
cor(grch38.default_exprs[,1:9],grch38.default_exprs[,10:18])

summary(pc.cr <- princomp(grch38.default_exprs, cor = TRUE))
print(summary(princomp(grch38.default_exprs, cor = TRUE),
              loadings = TRUE, cutoff = 0.2), digits = 2)

pc.cr <- princomp(grch38.default_exprs)  
pc.cr
grch38.pc<-princomp(grch38.default_exprs, cor = TRUE)
grch38.pc
biplot(princomp(grch38.default_exprs, cor = TRUE))
loadings(pc.cr) 
plot(pc.cr) 
biplot(grch38.pc)
princomp(~ ., data = grch38.default_exprs, cor = TRUE)
print(summary(princomp(grch38.default_exprs, cor = TRUE),
              loadings = TRUE, cutoff = 0.2), digits = 2)

# outliers CTRL_2 CTRL_5 LUTS_1 LUTS_5
# remove outliers and rownames (gene_short_name)
grch38.default_no_outs<-grch38.default_expr[,c(-1,-4,-7,-12,-16)]
dim(grch38.default_no_outs)
head(grch38.default_no_outs)
cor(grch38.default_no_outs[,c(1:3)],grch38.default_no_outs[,c(4:7)])
cor(grch38.default_no_outs[,8:11],grch38.default_no_outs[,12:14])
cor(grch38.default_no_outs[,1:7],grch38.default_no_outs[,8:14])

###########################################################################################################
# princomp w/o outliers
############################################################################################################

summary(pc.cr <- princomp(grch38.default_no_outs, cor = TRUE))
print(summary(princomp(grch38.default_no_outs, cor = TRUE),
              loadings = TRUE, cutoff = 0.2), digits = 2)

pc.cr <- princomp(grch38.default_no_outs)  
pc.cr
grch38.pc<-princomp(grch38.default_no_outs, cor = TRUE)
grch38.pc
biplot(princomp(grch38.default_no_outs, cor = TRUE))
loadings(pc.cr) 
plot(pc.cr) 
biplot(grch38.pc)
princomp(~ ., data = grch38.default_no_outs, cor = TRUE)
print(summary(princomp(grch38.default_no_outs, cor = TRUE),
              loadings = TRUE, cutoff = 0.2), digits = 2)

############################################################################################################
# data frames for each group
############################################################################################################

grch38.ctrlgrp<-grch38.expr[,grep(pattern="CTRL",colnames(grch38.expr))]
head(grch38.ctrlgrp)
dim(grch38.ctrlgrp)

grch38.lutsgrp<-grch38.default_exprs[,grep(pattern="LUTS",colnames(grch38.default_exprs))]
head(grch38.lutsgrp)
dim(grch38.lutsgrp)

############################################################################################################
# MDS comparison plots for all groups
############################################################################################################

## Create a design matrix:
fc$group = c('wt1021','wt1021','wt1021','wt1021','wt1021',
             'wt1021B','wt1021B','wt1021B','wt1021B','wt1021B',
             'A','A','A','A','A','AB','AB','AB','AB','AB')

groupsFactor <- factor(fc$group)
design <- model.matrix(~0+groupsFactor)
colnames(design) <- levels(groupsFactor)
cor(grch38.ctrlgrp,grch38.lutsgrp, method = "pearson")
grch38.lutsouts<-apply(grch38.lutsgrp, c(1,2), grubbs.test)


hg19.default.cuffdifftable<-file.path("cummeRbund_results_hg19_default/LUTS-over-CTRL/fpkm.matrix.csv")
hg19_default.cuff.features<-read.table(hg19.default.cuffdifftable,header = T,as.is=T)
hg19_default.cufffeatures<-hg19_default.cuff.features[,-1]
head(hg19_default.cufffeatures)
hg19_default.features<-hg19_default.cufffeatures[which(rowSums(hg19_default.cufffeatures[,-1])!=0),]

hg19.ctrlgrp<-hg19_default.features[,grep("CTRL",colnames(hg19_default.features))]
hg19.lutsgrp<-hg19_default.features[,grep("LUTS",colnames(hg19_default.features))]
cor(hg19.ctrlgrp,method = "spearman")
cor(hg19.lutsgrp,method = "spearman")

data <-read.table("../pnas_expression.txt",header=T,row.names=1)
dim(data)
head(data)
#get rid of the len column
data <- data[,-data$len]
class(data)
#check how many rows are all zeros
table(rowSums(data)==0)

FALSE  TRUE
21877 15558
#create a smaller subset as an example
data_subset <- data[rowSums(data)>500,]
dim(data_subset)
[1] 4010    7

#save row names
my_row <- rownames(data_subset)

#convert into matrix
data_subset_matrix <- as.matrix(data_subset)

#correlation of row 1 and row 2
cor(data_subset_matrix[1,],data_subset_matrix[2,],method="spearman")
[1] 0.5357143

#correlation matrix
correlation_matrix <- cor(t(data_subset_matrix), method="spearman")
correlation_matrix[1,2]
[1] 0.5357143
dim(correlation_matrix)
[1] 4010 4010

#output results
library(MASS)
write.matrix(correlation_matrix, file="pnas_expression_correlation.tsv")