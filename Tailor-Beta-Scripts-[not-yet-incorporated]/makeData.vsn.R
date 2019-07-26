
  csv1 = read.csv(gzfile("bacteria.csv"))
  drop = grep(".call", names(csv1))

  ##extract exprs - we find that some are neg - one of the
  ##issues with dChip - so let's run vsn on it and get onto
  ##a potentially better scale
  library(vsn)
  exprs = as.matrix(csv1[,-c(1:4, drop)])
  ex2 = vsn2(exprs)
  nexprs = predict(ex2, newdata=exprs)
  rownames(nexprs) = csv1[,1]

  ##now we need to rename a
  gds680Annot <-read.AnnotatedDataFrame("gds680-pData.txt")
  sampleNames(gds680Annot) = gds680Annot$name

  ##now reorder these experiments to be consistent with the
  ##annotation data
  nameMatch=read.csv(gzfile("nameMatch.csv"))
  m2 = match(gds680Annot$name, nameMatch[,2])
  nexprs = nexprs[,m2]
  colnames(nexprs) = gds680Annot$name

  bES = new("ExpressionSet", exprs=nexprs, phenoData=gds680Annot)

  ## cleaning gene Name
  ## Removing intergenic region, AFFX probe, other probe not related to well established ORFs
  notannotated <- grep("^IG|AFFX|^[A-Z]", featureNames(bES))
  bES <- bES[-notannotated, ]

  ## Simplify features names (affy ID) to ecoli genename
  ## as they are not unique I have add a column in the feature meta-data
  Ids <- read.table("gene-links.dat", skip=16 , sep="\t",
                    header=FALSE, as.is=TRUE)
  colnames(Ids) <- c("geneID", "EG", "B",  "yName",  "CGSC-ID",  "uniprotID",
                   "genename")
  fData(bES)[,1] <- gsub("*\\_[^_]*\\)*", "", featureNames(bES))
  matchingID <- Ids$genename[match(fData(bES)[,1], Ids$B)]
  fData(bES)[!is.na(matchingID), 1] <- matchingID[!is.na(matchingID)]
  fvarLabels(bES) <- c("genename")
  ## restrrict ourself to the single KO
  bES <- bES[, bES$genotype !="arcAfnr"]

  save(bES, file="bES.rda")

  source("http://bioconductor.org/biocLite.R")
  biocLite("org.EcK12.eg.db")
  ##not sure this works for us...
  library("org.EcK12.eg.db")
  accN = csv1[,3]
  EGs = unlist(mget(accN, org.EcK12.egACCNUM2EG, ifnotfound=NA))
  table(is.na(EGs))

  ##so I got the Affymetrix annotation
  ##first we rearrange
  AffyAnnot = read.csv("Ecoli_ASv2.na29.annot.csv")
  mminds = match(csv1[,1], AffyAnnot[,1])
  any(is.na(mminds))
  AA = AffyAnnot[mminds,]


##Next try out the bivariate stuff
##which TFs regulate the same genes and are the effects similar

 sp1 = split(aRes[,2], aRes[,1])
 for(i in 1:5) names(sp1[[i]]) = sp1[[i]]

 spFC = split(aRes[,3], aRes[,1])
 for(i in 1:5) names(spFC[[i]]) = sp1[[i]]

 share = matrix(NA, nr=5, nc=5)
 for(i in 1:5) {
    for(j in 1:5) {
        share[i,j] = sum(sp1[[i]] %in% sp1[[j]]) }}

 share

  ib = sp1[[1]] %in% sp1[[3]]
  nms1 = sp1[[1]][ib]
  plot(spFC[[1]][nms1], spFC[[3]][nms1])

  ib3 = sp1[[3]] %in% sp1[[4]]
  nms2 = sp1[[3]][ib3]
  plot(spFC[[3]][nms2], spFC[[4]][nms2])

  which( (spFC[[3]][nms2] < 0 ) & (spFC[[4]][nms2] > 0 )

  this one: flu_b2000_at  is down for fnr but up for oxyR...
  the only really unusual event here

