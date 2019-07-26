#functions for loading, normalising imagene microarray data and writing files using R
#This script was written for a particular person and was written quickly
#i.e. don't expect good coding practice or general functions in here.
#Microarray data functions used in this script are from Limma. Please reference Limma authors if you
#use their functions.
#
#List of functions used:
#processFiles   -   picks up info from Targest.txt and SpotTypes.txt and reads in data, background corrects,
#                   normalises (within and between if desired) and writes out data for each array into a separate file
#                   in a subdirectory called by default limmaOutput
#e.g. outList<-processFiles(targetsFile="FAL1targets.txt")   saves all limma objects into a list object called outList
#loadImageneData
#bkgCorr
#normWithin
#normBetween
#setContVals
#writeSepFiles
#writeMLinear

#Instructions:
#You can either run the overarching command "processFiles" after sourcing this file in R
#or you can run the "sub" functions - e.g. loadImageneData, by themselves
#Note that if you run the "sub" functions by themselves, you need to put all the options into the command line
#Defaults will not be picked up.

#assumptions: you will read in a targets file and spot types file from the current directory
#see Limma User Guide if you don't know what that means.

library(limma)

processFiles<-function(targetsFile="Targets.txt", spotTypesFile="SpotTypes.txt", bkmethod="subtract", offset=0,
                        normW="printtiploess", between=TRUE, normB="scale", write=TRUE, outDirName="limmaOutput")
{
  RGList<- loadImageneData(targetsFile, spotTypesFile)
  RGList_bkg<-bkgCorr(RGList, bkmethod, offset)
  MA_within<-normWithin(RGList_bkg, normW)
  
  numSlides<-length(colnames(MA_within$M))

  
  if(between)
  {
	MA_between<-normBetween(MA_within, normB)
	if(write)  { 
	
		writeSepFiles(MA_between, numSlides, outDirName )
		writeMLinearFile(MA_between, outDirName)
	
	   }
  }
  else
  {
	MA_between <- NULL
	if(write)  
	{ 
		writeSepFiles(MA_within, numSlides, outDirName) 
		writeMLinearFile(MA_within, outDirName)
	}
  }
  
 return(list(RGList=RGList, RGList_bkg=RGList_bkg, MA_within=MA_within, MA_between=MA_between))
}







######################Functions##################



loadImageneData<- function(targetsFile,spotTypesFile)
{

wt.fun=function(x) as.numeric(x$Flag==0)
targets<-readTargets(file=targetsFile)
files<-as.matrix(targets[,2:3])
spotTypes<-readSpotTypes(file=spotTypesFile)
RGList<-read.maimages(files,source="imagene", wt.fun=wt.fun)

RGList$genes$Status<-controlStatus(spotTypes,RGList)

RGList<-setContVals(RGList)

return(RGList)

}


bkgCorr<-function(RGList, bkmethod, offset)
{
	RGList_bkg <- backgroundCorrect(RGList, method=bkmethod, offset=offset)
	return(RGList_bkg)
}

normWithin<-function(RGList_bkg, normW)
{
 MA_within <- normalizeWithinArrays(RGList_bkg, method=normW)
 return(MA_within)

}

normBetween<-function(MA_within, normB)
{
 MA_between <- normalizeBetweenArrays(MA_within, method=normB)
 return(MA_between)
}



####################################
#set all spots with Status not equal to "gene" to a weight of 0
setContVals<-function(RGList) {

	for (i in 1:ncol(RGList))
	{
		RGList$weights[which(RGList$gene$Status != "gene"),i] <- 0
	}	
	
	return(RGList)
}




##########WRITE OUT DATA############

#assumption - you have loaded up a SpotTypes.txt files

writeSepFiles <- function(MAObj, numSlides, outDirName )
{

if (file.exists(outDirName)!= TRUE) dir.create(outDirName)
oldDir<-getwd()
setwd(outDirName)


indices<-which(MAObj$gene$Status == "gene")
filenames<-paste(colnames(MAObj$M),"_lma",sep="")
 
for (i in 1:numSlides)
 {
   write.table(data.frame(MAObj$gene$"Gene ID"[indices],
   cbind(2^(MAObj$M[indices,i]),
   MAObj$weights[indices,i])),
   file=filenames[i],
#colnames(MAObj$M)[i],
   row.names=FALSE,
   quote=FALSE, col.names=c("Gene ID", "M Value", "Flag"), sep="\t")
 }
 
setwd(oldDir)
}




writeMLinearFile<-function(MAObj, outDirName)
{
if (file.exists(outDirName)!= TRUE) dir.create(outDirName)
oldDir<-getwd()
setwd(outDirName)

indices<-which(MAObj$gene$Status == "gene")
  
write.table(cbind(as.character(MAObj$gene$"Gene ID"[indices]), (2^MAObj$M)[indices,]), quote=FALSE, sep="\t", row.names=FALSE, file="normOutput.txt")
  
setwd(oldDir)
}
