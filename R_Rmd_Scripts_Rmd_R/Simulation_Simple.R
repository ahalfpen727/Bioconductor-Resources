############################################################################
############################################################################
###									####
### This file provides R code to reproduce a subset of the simulations 	####
### presented in Figures 3, 4 and 5 of					####
###									####
### Law et al. "Voom: precision weights unlock linear model analysis 	####
### tools for RNA-seq read counts", Genome Biology 2014.		####
###									####
############################################################################
############################################################################


############################################################################
### The general settings of the simulation:
### Readers can change the followings for their own purposes

# Number of simulation runs
nsim <- 5

# Use equal or unequal library sizes
equal <- FALSE

# Use inverse chi-square or log-normal dispersion
invChisq <- TRUE

############################################################################

library(limma)
library(edgeR)

# Get distribution function of abundance proportions
# This distribution was generated from a real dataset
load(url("http://bioinf.wehi.edu.au/voom/qAbundanceDist.RData"))

# Generate baseline proportions for desired number of genes
ngenes <- 10000
baselineprop <- qAbundanceDist( (1:ngenes)/(ngenes+1) )
baselineprop <- baselineprop/sum(baselineprop)

# Design
group <- factor(c(1,1,1,2,2,2))
design <- model.matrix(~group)
n1 <- 3
n2 <- 3
nlibs <- n1+n2

# Library size 
if(equal){
	expected.lib.size <- rep(11e6,6)
} else {
	expected.lib.size <- 20e6 * c(1,0.1,1,0.1,1,0.1)
}

# Methods to be compared
methods <- c("voom", "limma trend", "limma notrend", "t-test", "edgeR classic", "edgeR glm")
nmethods <- length(methods)


######################################################################
###### Simulation to produce results in Figure 4 and 5
######################################################################

fd <- ranking <- matrix(0, nrow=ngenes, ncol=nmethods)
colnames(fd) <- colnames(ranking) <- methods
nd <- rep(0,nmethods)
names(nd) <- methods

# Set seed
set.seed(2013)
u <- runif(100)

# BEGIN SIM
for (k in 1:nsim) {
cat("SIM = ",k,"\n")

# Expected counts, group basis
i <- sample(1:ngenes,200)
i1 <- i[1:100]
i2 <- i[101:200]
fc <- 2
baselineprop1 <- baselineprop2 <- baselineprop
baselineprop1[i1] <- baselineprop1[i1]*fc
baselineprop2[i2] <- baselineprop2[i2]*fc
mu0.1 <- matrix(baselineprop1,ngenes,1) %*% matrix(expected.lib.size[1:n1],1,n1)
mu0.2 <- matrix(baselineprop2,ngenes,1) %*% matrix(expected.lib.size[(n1+1):(n1+n2)],1,n2)
mu0 <- cbind(mu0.1,mu0.2)
status <- rep(0,ngenes)
status[i1] <- -1
status[i2] <- 1

# Biological variation
BCV0 <- 0.2+1/sqrt(mu0)
if(invChisq){
	df.BCV <- 40
	BCV <- BCV0*sqrt(df.BCV/rchisq(ngenes,df=df.BCV))
} else {
	BCV <- BCV0*exp( rnorm(ngenes,mean=0,sd=0.25)/2 )
}
if(NCOL(BCV)==1) BCV <- matrix(BCV,ngenes,nlibs)
shape <- 1/BCV^2
scale <- mu0/shape
mu <- matrix(rgamma(ngenes*nlibs,shape=shape,scale=scale),ngenes,nlibs)

# Technical variation
counts <- matrix(rpois(ngenes*nlibs,lambda=mu),ngenes,nlibs)

# Filter
keep <- rowSums(counts)>=10
nkeep <- sum(keep)
counts2 <- counts[keep,]

# voom - ranked by lods
cat("voom\n")
i <- 1
y <- voom(counts2,design,plot=FALSE)
y$genes$Status <- status[keep]
fit <- lmFit(y,design)
fit <- eBayes(fit)
o <- order(fit$lods[,2], decreasing=TRUE)
ranking[1:nkeep,i] <- status[keep][o]
nd[i] <- nd[i] + sum(p.adjust(fit$p.value[,2],method="BH")<0.1)

# limma trend - ranked by lods
cat("limma trend\n")
i <- 2
y <- cpm(counts2,log=TRUE,prior.count=1)
fit <- lmFit(y,design,weights=NULL)
fit <- eBayes(fit,trend=TRUE)
o <- order(fit$lods[,2], decreasing=TRUE)
ranking[1:nkeep,i] <- status[keep][o]
nd[i] <- nd[i] + sum(p.adjust(fit$p.value[,2],method="BH")<0.1)

# limma notrend - ranked by lods
cat("limma notrend\n")
i <- 3
fit <- eBayes(fit,trend=FALSE)
o <- order(fit$lods[,2], decreasing=TRUE)
ranking[1:nkeep,i] <- status[keep][o]
nd[i] <- nd[i] + sum(p.adjust(fit$p.value[,2],method="BH")<0.1)

# t-test
cat("t-test\n")
i <- 4
t.ord <- fit$coef[,2]/fit$stdev.unscaled[,2]/fit$sigma
p.ord <- pt(abs(t.ord),df=4,lower.tail=FALSE)*2
fdr.ord <- p.adjust(p.ord,method="BH")
o <- order(p.ord)
ranking[1:nkeep,i] <- status[keep][o]
nd[i] <- nd[i] + sum(fdr.ord<0.1)

# edgeR classic
cat("edgeR classic\n")
i <- 5
z <- DGEList(counts=counts2,lib.size=colSums(counts2),group=group,genes=status[keep])
z2 <- estimateCommonDisp(z)
z2 <- estimateTagwiseDisp(z2)
et <- exactTest(z2)
o <- order(et$table$PValue)
ranking[1:nkeep,i] <- status[keep][o]
nd[i] <- nd[i] + sum(p.adjust(et$table$PValue,method="BH")<0.1)

# edgeR glm
cat("edgeR glm\n")
i <- 6
z <- estimateGLMTrendedDisp(z,design)
z <- estimateGLMTagwiseDisp(z,design)
fite <- glmFit(z,design)
lrt <- glmLRT(fite)
o <- order(lrt$table$PValue)
ranking[1:nkeep,i] <- status[keep][o]
nd[i] <- nd[i] + sum(p.adjust(lrt$table$PValue,method="BH")<0.1)

# Filtered genes at bottom
ranking[(nkeep+1):ngenes,] <- status[!keep]
fd <- fd + apply(1-abs(ranking),2,cumsum)

} # END SIM

fd <- fd/nsim
nd <- nd/nsim


# Plot Figure 5
i <- 1:120
fd.col <- c("black","blue","lightblue", "red", "blue", "black")
fd.type <- c(1,1,1,3,2,2)
fdmax <- apply(fd,1,max)
plot(i,fdmax[i],type="n", xlab="Genes chosen",ylab="False discoveries", main="Unequal library sizes, chi-square dispersions")
#axis(side=1, at=c(0,40,80,120))
#axis(side=2, at=c(0,20,40,60))
#box()
lines(i,fd[i,1],col="black",lwd=2)
lines(i,fd[i,2],col="blue",lwd=2)
lines(i,fd[i,3],col="lightblue",lwd=2)
lines(i,fd[i,4],col="red",lwd=2, lty=3)
lines(i,fd[i,5],col="blue",lwd=2,lty=2)
lines(i,fd[i,6],col="black",lwd=2, lty=2)
order <- 1:6
legend("topleft",legend=methods[order],lwd=2,col=fd.col[order],lty=fd.type[order], bty="n", cex=0.75)


# Plot Figure 4
nd.fd <- vector(length=nmethods)
for (i in c(1:length(nd.fd))) {
nd.fd[i] <- fd[ceiling(nd+0.00001)[i],i]
}
par(mar=c(7,5,1,1))
barplot(rbind(nd-nd.fd, nd.fd),
beside=FALSE, col=c("lightblue", "salmon"), legend.text=c("true positives", "false positives"), main="Unequal library sizes, chi-square dispersions",
args.legend=list(x="topleft",bty="n"), las=2, border=NA, axes=FALSE, ylab="Number of genes with FDR < 0.1", ylim=c(0,400))
axis(side=2, at=c(0,100,200,300,400), las=2)



######################################################################
###### Simulation to produce results in Figure 3
######################################################################

nd <- rep(0,nmethods)
names(nd) <- methods

# Set seed
set.seed(2013)
u <- runif(100)
cutoff <- 0.01

# BEGIN SIM
for (k in 1:nsim) {
cat("SIM = ",k,"\n")

# Expected counts, group basis
i <- sample(1:ngenes,200)
i1 <- i[1:100]
i2 <- i[101:200]
fc <- 2
baselineprop1 <- baselineprop2 <- baselineprop
mu0.1 <- matrix(baselineprop1,ngenes,1) %*% matrix(expected.lib.size[1:n1],1,n1)
mu0.2 <- matrix(baselineprop2,ngenes,1) %*% matrix(expected.lib.size[(n1+1):(n1+n2)],1,n2)
mu0 <- cbind(mu0.1,mu0.2)

# Biological variation
BCV0 <- 0.2+1/sqrt(mu0)
if(invChisq){
	df.BCV <- 40
	BCV <- BCV0*sqrt(df.BCV/rchisq(ngenes,df=df.BCV))
} else {
	BCV <- BCV0*exp( rnorm(ngenes,mean=0,sd=0.25)/2 )
}
if(NCOL(BCV)==1) BCV <- matrix(BCV,ngenes,nlibs)
shape <- 1/BCV^2
scale <- mu0/shape
mu <- matrix(rgamma(ngenes*nlibs,shape=shape,scale=scale),ngenes,nlibs)

# Technical variation
counts <- matrix(rpois(ngenes*nlibs,lambda=mu),ngenes,nlibs)

# Filter
keep <- rowSums(counts)>=10
nkeep <- sum(keep)
counts2 <- counts[keep,]

# voom - ranked by lods
cat("voom\n")
i <- 1
y <- voom(counts2,design,plot=FALSE)
fit <- lmFit(y,design)
fit <- eBayes(fit)
nd[i] <- nd[i] + sum(fit$p.value[,2]<cutoff)

# limma trend - ranked by lods
cat("limma trend\n")
i <- 2
y <- cpm(counts2,log=TRUE,prior.count=1)
fit <- lmFit(y,design,weights=NULL)
fit <- eBayes(fit,trend=TRUE)
o <- order(fit$lods[,2], decreasing=TRUE)
ranking[1:nkeep,i] <- status[keep][o]
nd[i] <- nd[i] + sum(fit$p.value[,2]<cutoff)

# limma notrend - ranked by lods
cat("limma notrend\n")
i <- 3
fit <- eBayes(fit,trend=FALSE)
o <- order(fit$lods[,2], decreasing=TRUE)
ranking[1:nkeep,i] <- status[keep][o]
nd[i] <- nd[i] + sum(fit$p.value[,2]<cutoff)

# t-test
cat("t-test\n")
i <- 4
t.ord <- fit$coef[,2]/fit$stdev.unscaled[,2]/fit$sigma
p.ord <- pt(abs(t.ord),df=4,lower.tail=FALSE)*2
nd[i] <- nd[i] + sum(p.ord<cutoff)

# edgeR classic
cat("edgeR classic\n")
i <- 5
z <- DGEList(counts=counts2,lib.size=colSums(counts2),group=group,genes=status[keep])
z2 <- estimateCommonDisp(z)
z2 <- estimateTagwiseDisp(z2)
et <- exactTest(z2)
nd[i] <-nd[i] + sum(et$table$PValue<cutoff)

# edgeR glm
cat("edgeR glm\n")
i <- 6
z <- estimateGLMTrendedDisp(z,design)
z <- estimateGLMTagwiseDisp(z,design)
fite <- glmFit(z,design)
lrt <- glmLRT(fite)
nd[i] <- nd[i] + sum(lrt$table$PValue<cutoff)

} # END SIM

nd <- nd/nsim

# Plot Figure 3
nd <- nd/ngenes
par(mar=c(7,5,1,1))
barplot(nd, ylab="Percentage of genes with p-value < 0.01", las=2, border=NA, col="lightblue", axes=FALSE, ylim=c(0,0.05))
axis(side=2, at=c(0,0.01,0.02,0.03, 0.04, 0.05), las=2)
abline(h=0.01, col="red")


