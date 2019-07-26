### R code from vignette source 'sigPathway-vignette.Rnw'
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("sigPathway")
vignette(sigPathway)
data(MuscleExample)
ls()



###################################################
### code chunk number 2: sigPathway-vignette.Rnw:136-140
###################################################
dim(tab)
print(tab[501:504, 1:3])
table(phenotype)



###################################################
### code chunk number 3: PSIDhist
###################################################
statList <- calcTStatFast(tab, phenotype, ngroups = 2)
hist(statList$pval, breaks = seq(0,1,0.025), xlab = "p-value",
     ylab = "Frequency", main = "")


###################################################
### code chunk number 4: sigPathway-vignette.Rnw:164-170
###################################################
set.seed(1234)
res.muscle <-
  runSigPathway(G, 20, 500, tab, phenotype, nsim = 1000,
                weightType = "constant", ngroups = 2, npath = 25, 
                verbose = FALSE, allpathways = FALSE, annotpkg = "hgu133a.db",
                alwaysUseRandomPerm = FALSE)


###################################################
### code chunk number 5: sigPathway-vignette.Rnw:226-227
###################################################
print(res.muscle$df.pathways[1:10, ])


###################################################
### code chunk number 6: sigPathway-vignette.Rnw:243-244
###################################################
print(res.muscle$list.gPS[[7]][1:10, ])


###################################################
### code chunk number 7: sigPathway-vignette.Rnw:280-281
###################################################
print(sessionInfo())


