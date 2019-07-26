### R code from vignette source 'tutorial.Rnw'
### Encoding: UTF-8
#install.packages("rJava") # if not present already
#install.packages("RCurl") # if not present already
#install.packages("devtools") # if not present already
#install_github("BiGCAT-UM/bridgedb-r", subdir="BridgeDbR")
library(devtools); library(RCurl); library(BridgeDbR); library(rJava)
Hs_code = getOrganismCode("Homo sapiens")
Hs_code
dbLocation <- getDatabase("Homo sapiens",location=getwd())

fullName <- getFullName("Ce")
fullName
code <- getSystemCode("ChEBI")
code

getMatchingSources("HMDB00555")
getMatchingSources("ENSG00000100030")

getBridgeNames()
dbLocation <- getDatabase("Homo sapiens",location=getwd())
mapper <- loadDatabase(dbLocation)
###################################################
### code chunk number 8: tutorial.Rnw:115-118 (eval = FALSE)
###################################################
## location <- getDatabase("Homo sapiens")
## mapper <- loadDatabase(location)
## map(mapper, "L", "196410", "X")


###################################################
### code chunk number 7: tutorial.Rnw:105-106 (eval = FALSE)
###################################################
##


###################################################
### code chunk number 8: tutorial.Rnw:115-118 (eval = FALSE)
###################################################
## location <- getDatabase("Homo sapiens")
## mapper <- loadDatabase(location)
## map(mapper, "L", "196410", "X")

Hs_code = getOrganismCode("Homo sapiens")
Hs_code

fullName <- getFullName("Ce")
fullName
code <- getSystemCode("ChEBI")
code

getMatchingSources("HMDB00555")
getMatchingSources("ENSG00000100030")

getBridgeNames()
dbLocation <- getDatabase("Homo sapiens",location=getwd())
mapper <- loadDatabase(dbLocation)




