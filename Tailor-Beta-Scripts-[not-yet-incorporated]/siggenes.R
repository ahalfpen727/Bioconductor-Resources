# package for identifying FDR and significantly differentially expressed genes
library(siggenes)
data(golub)
args(sam)
args(d.stat)
args.sam(summary)

###################################################
### code chunk number 7: siggenes.Rnw:264-266
###################################################
n <- 10
rep(1, n)
n1 <- n2 <- 5
rep(c(0, 1), c(n1, n2))
K <- 5
c((-1:-5), 1:5)
K <- 5
rep(1:K, e = 2) * rep(c(-1 ,1), K)
K <- 5
cbind(rep(c(-1, 1), 5), rep(1:5, e = 2))


###################################################
### code chunk number 12: siggenes.Rnw:374-376
###################################################
sam.out <- sam(golub, golub.cl, rand = 123, gene.names = golub.gnames[,3])
sam.out
sam.out2 <- sam(golub, golub.cl, method = wilc.stat, rand = 123)
summary(sam.out)
print(sam.out, seq(1.5, 2.4, 0.1))
sum.sam.out <- summary(sam.out, 3.3)
sum.sam.out
print(sum.sam.out, varNames = "Proteins")
sum.sam.out@row.sig.genes
sum.sam.out@mat.fdr
sum.sam.out@mat.sig

list.siggenes(sam.out, 3.3)
findDelta(sam.out, fdr = 0.05)
findDelta(sam.out, genes = 200)
find.out <- find.a0(golub, golub.cl, rand = 123)
find.out
print(find.out, 0.95)

ebam(find.out)
ebam(find.out, which.a0 = 2)
ebam(golub, golub.cl, a0 = 0, fast = TRUE, rand = 123)
ebam.out <- ebam(golub, golub.cl, a0 = 0, rand = 123)


###################################################
### code chunk number 31: siggenes.Rnw:651-652
###################################################
print(ebam.out, seq(0.91, 0.99, 0.01))
summary(ebam.out, 0.99997)
ebam(golub, golub.cl, a0 = 0, var.equal = TRUE, rand = 123)
ebam(golub, golub.cl, quan.a0 = 0.5, rand = 123)
ebam(golub, golub.cl, method = wilc.ebam, rand =123)
t.stat <- function(data, cl){
    require(genefilter) ||
        stop("genefilter required.")
    cl <- as.factor(cl)
    row.out <- rowttests(data, cl)
    d <- row.out$statistic
    m <- length(na.exclude(d))
    d.bar <- qt(((1:m) - 0.5)/m, length(cl) - 2)
    p.value <- row.out$p.value
    vec.false <- m * p.value/2
    s <- row.out$dm/d
    msg <- paste("SAM Two-Class Analysis",
         "Assuming Normality\n\n")
    list(d = -d, d.bar = d.bar, p.value = p.value,
        vec.false = vec.false, s = s, s0 = 0,
        mat.samp = matrix(numeric(0)),
        msg = msg, fold = numeric(0))}
sam(golub, golub.cl, method = t.stat)


###################################################
### code chunk number 38: siggenes.Rnw:854-871
###################################################
t.find <- function(data, cl, B = 50){
    require(genefilter)
    z.fun <- function(data, cl){
        cl <- as.factor(cl)
        out <- rowttests(data, cl)
        r<- out$dm
        s<- r / out$statistic
        return(list(r = -r, s = s))
    }
    mat.samp <- matrix(0, B, length(cl))
    for(i in 1:B)
        mat.samp[i, ] <- sample(cl)
    z.out <- z.fun(data, cl)
    msg <- paste("EBAM Analysis with a Moderated t-Statistic\n\n")
    list(r = z.out$r, s = z.out$s,
        mat.samp = mat.samp, z.fun = z.fun, msg = msg)
}


###################################################
### code chunk number 39: siggenes.Rnw:881-883
###################################################
t.out <- find.a0(golub, golub.cl, method = t.find, B = 100, rand =123)
t.out
find.a0(golub, golub.cl, var.equal = TRUE, rand =123)
ebam(t.out)


###################################################
### code chunk number 42: siggenes.Rnw:928-941
###################################################
t.ebam<-function(data, cl){
    require(genefilter)
    cl <- as.factor(cl)
    out <- rowttests(data, cl)
    z <- -out$statistic
    z.dens <- denspr(z)$y
    m <- length(z)
    vec.pos <- m * out$p.value / 2
    z.null <- dt(z, length(cl) - 2)
    msg<-paste("EBAM Analysis with t-Statistic Assuming Normality.\n\n")
    list(z = z, ratio = z.null/z.dens, vec.pos = vec.pos,
        vec.neg = vec.pos, msg = msg)}
ebam(golub, golub.cl, method = t.ebam)


