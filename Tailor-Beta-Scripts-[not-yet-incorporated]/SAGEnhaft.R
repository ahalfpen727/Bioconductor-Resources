### R code from vignette source 'SAGEnhaft.Rnw'
source("https://bioconductor.org/biocLite.R")
biocLite("sagenhaft")
###################################################
### code chunk number 1: SAGEnhaft.Rnw:121-122
###################################################
library(sagenhaft)


###################################################
### code chunk number 2: SAGEnhaft.Rnw:134-142
###################################################
E15post <- read.sage.library(system.file("extdata", "E15postHFI.sage",
                            package="sagenhaft"))
E15post
B6Hypo <- read.sage.library(system.file("extdata", "B6HypothalHFI.sage",
                            package="sagenhaft")) 
libcomp <- compare.lib.pair(B6Hypo, E15post)
plot(libcomp)
libcomp


