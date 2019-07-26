#drosophila analysis of Rh8_Rh5 [2 samples] against Rh8_Rh6 [2 samples]
#######################################################################
### code chunk number 1: Library Load
#######################################################################
library(igraph);library(GEOquery);library(Hmisc)

#######################################################################
### code chunk number 2: Load the count matrix
#######################################################################

files<-dir()
files
gse116969<-grep(pattern="Rister_",files, ignore.case = T)
gse116969.xlsx<-files[gse116969]
gse116969.xlsx<-readxl::read_xlsx(gse116969.xlsx)
gse116969.xlsx<-as.data.frame(gse116969.xlsx)
colnames(gse116969.xlsx)

transcripts<-unique(gse116969.xlsx[,"transcript_id"])
row.names(gse116969.xlsx)=transcripts
head(gse116969.xlsx)

gse116969<-as.matrix(gse116969.xlsx[,-c(1:3)])
samps<-colnames(gse116969)
genes<-row.names(gse116969)
head(gse116969)
dim(gse116969)

Rh6samps<-grep(pattern="Rh6", colnames(gse116969), ignore.case = T)
Rh5samps<-grep(pattern="Rh5", colnames(gse116969), ignore.case = T)
head(gse116969[,Rh6samps])
head(gse116969[,Rh5samps])

Rh8_Rh6s<-factor(grep(pattern="Rh6", colnames(gse116969), ignore.case = T))
Rh8_Rh5s<-factor(grep(pattern="Rh5", colnames(gse116969), ignore.case = T))
Rh8_Rh6group<-gse116969[,Rh8_Rh6s]
Rh8_Rh5group<-gse116969[,Rh8_Rh5s]

##############################################################
### code chunk number 3: Design Matrix and pair like-samples
##############################################################

gseDesign = data.frame(row.names = colnames(gse116969),
                       condition = c("Rh5", "Rh5", "Rh6","Rh6"))
conditions<-factor(gseDesign$condition)

Rh5group = gseDesign$condition == "Rh5"
Rh6group = gseDesign$condition == "Rh6"

Rh5samples.gene_exp = gse116969[,Rh5group]
Rh6samples.gene_exp = gse116969[,Rh6group]
head(Rh5samples.gene_exp)
head(Rh6samples.gene_exp)

rh5.dist<-dist(Rh5samples.gene_exp)
rh6.dist<-dist(Rh6samples.gene_exp)

r6dist <- as.matrix(dist(Rh6samples.gene_exp))

## Remove the diagonal from consideration
diag(r6dist) <- diag(r6dist) + 100000
 
# Find the index of the points with minimum distance
ind <- which(r6dist == min(r6dist), arr.ind = TRUE)
ind
  
wilcox.test(Rh5samples.gene_exp ~ Rh6samples.gene_exp, paired=TRUE)

###################################################
### code chunk number 4: Descriptive stats
###################################################

gse.sd<-apply(gse116969, 1, sd)
plot(gse.sd)
summary(gse.sd)
hist(gse.sd,breaks = 300,include.lowest = T,freq = T)

gse116969.sd.df<-cbind(gse116969, gse.sd)
head(gse116969.sd.df)

gse.sorted.sd<-order(gse.sd,decreasing = T)
names(gse.sorted.sd)
hist(gse.sorted.sd)

Rh5var<-apply(Rh5samples.gene_exp, 1,var)
Rh6var<-apply(Rh6samples.gene_exp, 1,var)
summary(Rh5var)
summary(Rh6var)
head(Rh5samples.gene_exp)
head(Rh6samples.gene_exp)

Rh5sd<-apply(Rh5samples.gene_exp, 1,sd)
Rh6sd<-apply(Rh6samples.gene_exp, 1,sd)
summary(Rh5sd)
summary(Rh6sd)
head(Rh5sd)
head(Rh6sd)

qqplot(Rh5sd, Rh6sd)
plot(Rh5sd, Rh6sd)

###################################################
### code chunk number 5: pairedSamples
###################################################

gse116969.sd.var.df<-cbind(gse116969.sd.df, Rh5var)
head(gse116969.sd.var.df)
gse116969.sd.var.df<-cbind(gse116969.sd.var.df, Rh6var)
head(gse116969.sd.var.df)


Rh5mean.expr<-apply(Rh5samples.gene_exp, 1,mean)
Rh6mean.expr<-apply(Rh6samples.gene_exp, 1,mean)
head(Rh5mean.expr)
head(Rh6mean.expr)
Rh8mean.expr<-cbind(Rh5mean.expr,Rh6mean.expr)
head(Rh8mean.expr)

logbase2fold<-log(x=c(Rh5mean.expr/Rh6mean.expr), base=2)
plot(logbase2fold)
is(logbase2fold)

RH8expr<-cbind(Rh8mean.expr,logbase2fold)
head(RH8expr)
sorted.logfold<-order(logbase2fold, decreasing = T)
head(RH8expr[sorted.logfold,])

siggenes<-subset(Rh8mean.expr,abs(logbase2fold) >= 100 & rowSums(Rh8mean.expr) >= 10)
head(siggenes)

siggse<-row.names(gse116969)  %in%  row.names(siggenes)
sig.genes.gse116969<-gse116969[siggse,]
head(sig.genes.gse116969)
heatmap(sig.genes.gse116969)


sig.dist<-as.matrix(dist(sig.genes.gse116969))
diag(sig.dist) <- diag(sig.dist) + 100000

ind <- which(sig.dist == min(sig.dist), arr.ind = TRUE)
ind
##################################################
### code chunk number : Step 6: Statistical Tests
##################################################

var.test(RH8expr[,"Rh5mean.expr"],RH8expr[,"Rh6mean.expr"])
t.test(RH8expr[,"Rh5mean.expr"],RH8expr[,"Rh6mean.expr"],var.equal=T)
gsecor<-cor(gse116969)
gse.genescor<-cor(t(gse116969[]))


t.test(x = Rh5mean.expr,
       y = Rh6mean.expr,
       alternative = "two.sided")

t.test(Rh5var.expr ~ Rh6var.expr, data = df)

####################################################################################################
### code chunk number : Step 3: #Hierarchical clustering
####################################################################################################
library(pvclust); library(gplots) # Loads the required packages.
#Hierarchical clustering
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R")
# Import an alternative color scheme for the heatmap function.
mydatascale <- t(scale(t(gse116969))) # Centers and scales data.
hr <- hclust(as.dist(1-cor(t(gse116969), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
hc <- hclust(as.dist(1-cor(gse116969, method="spearman")), method="complete")
# Clusters columns by Spearman correlation.
heatmap(gse116969, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=my.colorFct(), scale="row")
# Plot the data table as heatmap and the cluster results as dendrograms.
mycl <- cutree(hr, h=max(hr$height)/1.5)
mycolhc <- sample(rainbow(256))
mycolhc <- mycolhc[as.vector(mycl)]
heatmap(gse116969, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), 
        col=my.colorFct(), scale="row", RowSideColors=mycolhc)
# Cut the tree at specific height and color the corresponding clusters in the heatmap color bar.

pv <- pvclust(scale(t(gse116969)), method.dist="correlation", method.hclust="complete", nboot=10)
# Perform the hierarchical cluster analysis. 
#Due to time resrictions, we are using here only 10 bootstrap repetitions.
# Usually, one should use at least 1000 repetitions.
plot(pv, hang=-1); pvrect(pv, alpha=0.95)
# Plots result as a dendrogram where the significant clusters are highlighted with red rectangles.
clsig <- unlist(pvpick(pv, alpha=0.95, pv="au", type="geq", max.only=TRUE)$clusters) # Retrieve members of significant clusters.

source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/dendroCol.R") # Import tree coloring function.
dend_colored <- dendrapply(as.dendrogram(pv$hclust), dendroCol, keys=clsig, xPar="edgePar", bgr="black", fgr="red", pch=20)
# Create dendrogram object where the significant clusters are labeled in red.
heatmap(sig.g.rep.matrix, Rowv=dend_colored, Colv=as.dendrogram(hc), col=my.colorFct(), scale="row", RowSideColors=mycolhc)
# Plot the heatmap from above, but with the significant clusters in red and the cluster bins from the tree cutting step in
# the color bar.
x11(height=12); heatmap.2(mydata, Rowv=dend_colored, Colv=as.dendrogram(hc), col=my.colorFct(), scale="row", trace="none", RowSideColors=mycolhc) # Plot heatmap with heatmap.2() function which scales better for many entries.
mydatasort <- mydata[pv$hclust$labels[pv$hclust$order], hc$labels[hc$order]] # Sort rows in data table by 'dend_colored' and its colums by 'hc'.
x11(height=16, width=12); par(mfrow=c(1,2)); plot(dend_colored, horiz=T, yaxt="n"); image(scale(t(mydatasort)), col=my.colorFct(), xaxt="n",yaxt="n") # Plot heatmap with bootstrap tree in larger format using instead of heatmap the image function.
pdf("pvclust.pdf", height=21, width=6); plot(dend_colored, horiz=T, yaxt="n"); dev.off(); pdf("heatmap.pdf", height=20, width=6); image(scale(t(mydatasort)), col=my.colorFct(), xaxt="n",yaxt="n"); dev.off()


####################################################################################################
### code chunk number : Step 2, Read in, normalise, and identify genes with significant effects
####################################################################################################


gse.scaled<-as.matrix(scale(gse116969))
heatmap(gse.scaled, Colv=F, scale='none')

## Box Plot ##
summary(gse116969)
pairs(gse116969)

## Box Plot ##
## ANOVA
summary(aov(R8_Rh5_d1_rep1+R8_Rh5_d1_rep2~R8_Rh6_d1_rep1+R8_Rh6_d1_rep2, data=gse116969) )
R8_Rh5_d1_rep1+R8_Rh5_d1_rep2 R8_Rh6_d1_rep1 R8_Rh6_d1_rep2


gsecor<-cor(gse116969)
round(gsecor,2)
# Do the same, but with colors corresponding to value
plot(hclust(dist(gse116969,method="euclidian"),method="single"))
cl <- kmeans(gse116969, 2)
plot(gse116969, col = cl$cluster)
points(cl$centers, col = 1:2, pch = 8, cex=2)
## Regression ##

model1 = lm(Life.Exp  ~ Population + Income + Illiteracy + Murder + HS.Grad + Frost + Area , data=state.x77)

## make data frame object

st = as.data.frame(state.x77)
dim(st)
## Fit Linear Regression model
model1 = lm(Life.Exp ~ Population + Income + Illiteracy + Murder +
              HS.Grad + Frost + Area , data=st)

summary(model1)
## remove non significant variable
model2 = update(model1, .~.-Population -Illiteracy -Income - Area)
summary(model2)

## Prediction
predict(model2, list(Murder=10, HS.Grad=50, Frost=90))



## Logistic Regression ##

## load data ##

bindata <- read.csv("http://www.ats.ucla.edu/stat/data/binary.csv")

## view the first few rows of the data
head(bindata)

## view the last few rows of the data
tail(bin data)

# obtain basic statistics

summary(bindata)
sapply(bindata, sd)


bindata$rank <- factor(bindata$rank)

## fit logistics model

modlogit <- glm(admit ~ gre + gpa + rank, data = bindata, family = "binomial")

## result
summary(modlogit)
