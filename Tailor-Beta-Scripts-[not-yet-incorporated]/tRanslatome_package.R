### R code from vignette source 'tRanslatome_package.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: tRanslatome_package.Rnw:55-75
###################################################
 ##loading the tRanslatome package
 library(tRanslatome) 
 ##loading the training data set
 data(tRanslatomeSampleData) 
 translatome.analysis <- newTranslatomeDataset(expressionMatrix, 
			c("tot.undiff.a", "tot.undiff.b", "tot.undiff.c"), 
			c("tot.diff.a", "tot.diff.b", "tot.diff.c"), 
			c("pol.undiff.a", "pol.undiff.b", "pol.undiff.c"), 
			c("pol.diff.a", "pol.diff.b", "pol.diff.c"), 
			label.level= c("transcriptome", "translatome"), 
			label.condition=c("undifferentiated", "differentiated"))
 ##identification of DEGs with the use of the limma statistical method
 limma.DEGs <- computeDEGs(translatome.analysis, 
													 method= "limma", mult.cor=TRUE) 
 ##enrichment analysis of the selected DEGs
 CCEnrichment <- GOEnrichment(limma.DEGs,ontology="CC", classOfDEGs="up",
															test.method="elim", test.threshold = 0.05)
 ##performing a comparison of the biological themes enriched 
 ##in the two levels of gene expression
 CCComparison <- GOComparison(CCEnrichment) 


