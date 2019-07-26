### R code from 
vignette(spliceR)
### Encoding: UTF-8

library("spliceR")
cuffDB_spliceR <- prepareCuff(cuff)
xx<-prepareCuff(cuff, fixCufflinksAnnotationProblem=TRUE,removeNonChanonicalChr=TRUE)
xx

Always, load from specified RData directly
SCOP.sf <- dcRDataLoader(RData='SCOP.sf')

'SCOP.sf' (from package 'dcGOR' version 1.0.5) has been loaded into the working environment
Pfam <- dcRDataLoader(RData='Pfam')

'Pfam' (from package 'dcGOR' version 1.0.5) has been loaded into the working environment
InterPro <- dcRDataLoader(RData='InterPro')

'InterPro' (from package 'dcGOR' version 1.0.5) has been loaded into the working environment
Rfam <- dcRDataLoader(RData='Rfam')

'Rfam' (from package 'dcGOR' version 1.0.5) has been loaded into the working environment
onto.GOMF <- dcRDataLoader(RData='onto.GOMF')

'onto.GOMF' (from package 'dcGOR' version 1.0.5) has been loaded into the working environment
# But for annotaion data, there are two ways to do so:
# 1) in a direct way
SCOP.sf2GOMF <- dcRDataLoader(RData='SCOP.sf2GOMF')

'SCOP.sf2GOMF' (from package 'dcGOR' version 1.0.5) has been loaded into the working environment
# 2) in an indirect way: specify both domain and ontology
SCOP.sf2GOMF <- dcRDataLoader(domain='SCOP.sf', ontology='GOMF')

###################################################
### code chunk number 2: Helper_functions
###################################################
myTranscripts <- transcripts(cuffDB_spliceR)
myExons <- exons(cuffDB_spliceR)
conditions(cuffDB_spliceR)
session <- browserSession("UCSC")
genome(session) <- "hg19"
query <- ucscTableQuery(session, "knownGene")
tableName(query) <- "knownGene"
cdsTable <- getTable(query)
tableName(query) <- "kgXref"
kgXref <- getTable(query)


###################################################
### code chunk number 5: GRanges_3
###################################################
knownGeneTranscripts <- GRanges(
	seqnames=cdsTable$"chrom",
	ranges=IRanges(
		start=cdsTable$"txStart",
		end=cdsTable$"txEnd"),
	strand=cdsTable$"strand",
	spliceR.isoform_id = cdsTable$"name",
	spliceR.sample_1="LUTS",   # placeholder1
	spliceR.sample_2="CTRL",      # placeholder2
	spliceR.gene_id=kgXref[match(cdsTable$"name", kgXref$"kgID"), 
		"geneSymbol"],
	spliceR.gene_value_1=1,
	spliceR.gene_value_2=1,
	spliceR.gene_log2_fold_change=log2(1/1),
	spliceR.gene_p_value=1,
	spliceR.gene_q_value=1,
	spliceR.iso_value_1=1,
	spliceR.iso_value_2=1,
	spliceR.iso_log2_fold_change=log2(1/1),
	spliceR.iso_p_value=1,
	spliceR.iso_q_value=1
	)


###################################################
### code chunk number 6: GRanges_4
###################################################
startSplit	<- strsplit(as.character(cdsTable$"exonStarts"), split=",")
endSplit  	<- strsplit(as.character(cdsTable$"exonEnds"), split=",")
startSplit  <- lapply(startSplit, FUN=as.numeric)
endSplit	<- lapply(endSplit, FUN=as.numeric)

knownGeneExons <- GRanges(
	seqnames=rep(cdsTable$"chrom", lapply(startSplit, length)),
	ranges=IRanges(
		start=unlist(startSplit)+1,
		end=unlist(endSplit)),
	strand=rep(cdsTable$"strand", lapply(startSplit, length)),
	spliceR.isoform_id=rep(knownGeneTranscripts$"spliceR.isoform_id",
	lapply(startSplit, length)),
	spliceR.gene_id=rep(knownGeneTranscripts$"spliceR.gene_id",
	lapply(startSplit, length))
	)


###################################################
### code chunk number 7: GRanges_5
###################################################
knownGeneSpliceRList <- SpliceRList(
	transcript_features=knownGeneTranscripts,
	exon_features=knownGeneExons,
	assembly_id="hg19",
	source="granges",
	conditions=c("LUTS", "CTRL")
	)


###################################################
### code chunk number 8: preSpliceRFilter
###################################################
cuffDB_spliceR_filtered <- preSpliceRFilter(
	cuffDB_spliceR,
	filters=c("expressedIso", "isoOK", "expressedGenes", "geneOK")
	)


###################################################
### code chunk number 9: spliceR
###################################################
#Commented out due to problems with vignettes and progress bars.
mySpliceRList <- spliceR(
	cuffDB_spliceR,
	compareTo="preTranscript",
	filters=c("expressedGenes","geneOK", "isoOK", "expressedIso", "isoClass"),
	useProgressBar=F
)


###################################################
### code chunk number 10: annotatePTC
###################################################
ucscCDS <- getCDS(selectedGenome="hg19", repoName="UCSC")

#Annotate with PTCs

PTCSpliceRList <- annotatePTC(cuffDB_spliceR, cds=ucscCDS, Hsapiens, 
                              PTCDistance=50)

###################################################
### code chunk number 11: generateGTF
###################################################
generateGTF(mySpliceRList, filters=
c("geneOK", "isoOK", "expressedGenes", "expressedIso"), 
	scoreMethod="local",
	useProgressBar=F)

###################################################
### code chunk number 12: plot1
###################################################
#Plot the average number of transcripts pr gene
mySpliceRList <- spliceRPlot(mySpliceRList, 
	evaluate="nr_transcript_pr_gene")


###################################################
### code chunk number 13: plot2
###################################################
#Plot the average number of Exon skipping/inclucion evens pr gene
mySpliceRList <- spliceRPlot(mySpliceRList, evaluate="mean_AS_transcript", 
	asType="ESI")


###################################################
### code chunk number 14: plot3
###################################################
#Plot the average number of Exon skipping/inclucion evens pr gene, 
#but only using transcripts that are significantly differntially expressed
mySpliceRList <- spliceRPlot(mySpliceRList, evaluate="mean_AS_transcript", 
	asType="ESI", filters="sigIso",reset=TRUE)

#Total number of alternatively spliced events
txyz<-totalNumberOfAS(mySpliceRList)

#Get top dIF transcripts
topIsoShift(mySpliceRList, n=20)

#Plot number of exon skipping/inclusion events 
mySpliceRList <- spliceRPlot(mySpliceRList, evaluate="nr_AS", asType="ESI")
