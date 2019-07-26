## ---- message=FALSE------------------------------------------------------
library(QuaternaryProd)

# Compute the probability mass function
pmf <- QP_Pmf(q_p = 20, q_m = 20, q_z = 20, q_r = 0, n_p = 20, n_m = 20, n_z = 20)

# Plot the mass function
plot(names(pmf), pmf, col="blue", xlab = "scores", ylab = "probabilities")
lines(names(pmf), pmf, col = "blue")

## ---- message=FALSE------------------------------------------------------
# Get the p-value of score 5
pval <- QP_Pvalue(score = 5, q_p = 20, q_m = 20, q_z = 20, q_r = 0, 
                                                     n_p = 20, n_m = 20, n_z = 20)
pval

# Compue the p-value only if it is statistically significant otherwise
# return -1
pval <- QP_SigPvalue(score = 5, q_p = 20, q_m = 20, q_z = 20, q_r = 0, 
                                                     n_p = 20, n_m = 20, n_z = 20)
pval

## ---- message=FALSE------------------------------------------------------
library(QuaternaryProd)

# Get gene expression data
e2f3 <- system.file("extdata", "e2f3_sig.txt",
                             package = "QuaternaryProd")
e2f3 <- read.table(e2f3, sep = "\t",
                             header = TRUE, stringsAsFactors = FALSE)
myc <- system.file("extdata", "myc_sig.txt",
                             package = "QuaternaryProd")
myc <- read.table(myc, sep = "\t",
                             header = TRUE, stringsAsFactors = FALSE)

# Rename column names appropriately
# and remove duplicated entrez ids in the gene expression data
names(e2f3) <- c("entrez", "pvalue", "fc")
e2f3 <- e2f3[!duplicated(e2f3$entrez),]

names(myc) <- c("entrez", "pvalue", "fc")
myc <- myc[!duplicated(myc$entrez),]

## ---- message=FALSE------------------------------------------------------
# Compute the Quaternary Dot Product Scoring Statistic for only statistically
# significant regulators
quaternary_results <- RunCRE_HSAStringDB(e2f3, method = "Quaternary",
                                   fc.thresh = log2(1.3), pval.thresh = 0.05,
                                   only.significant.pvalues = TRUE, 
                                   significance.level = 0.05,
                                   epsilon = 1e-16)
quaternary_results[1:4, c("uid","symbol","regulation","pvalue")]

## ---- message=FALSE------------------------------------------------------
ternary_results <- RunCRE_HSAStringDB(myc, method = "Ternary",
                                      fc.thresh = log2(1.3), pval.thresh = 0.05,
                                      only.significant.pvalues = TRUE,
                                      significance.level = 0.05,
                                      epsilon = 1e-16)
ternary_results[1:4, c("uid","symbol","regulation","pvalue")]

## ---- message=FALSE------------------------------------------------------
enrichment_results <- RunCRE_HSAStringDB(myc, method = "Enrichment",
                                         fc.thresh = log2(1.3), pval.thresh = 0.05,
                                         only.significant.pvalues = TRUE,
                                         significance.level = 0.05,
                                         epsilon = 1e-16)
enrichment_results[1:10, c("uid","symbol","regulation","pvalue")]

