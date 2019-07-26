library("openxlsx")
rm(list=ls(all=TRUE))



path1 = "E:\\HF_AKRJ-over-LF_AKRJ"
diffFiles = list.files(path=path1, pattern="\\.diff$")

#rawd <- RawData

Myfunction <- function(rawd){
  colname <- colnames(rawd, do.NULL = TRUE, prefix = "col")
  if(colname[10]== "log2.fold_change."){
    NewColname <- c(colname[1:10],"Fold Change",colname[11:14])
    NewMatrix<- matrix( nrow=dim(rawd)[1], ncol=(dim(rawd)[2]+1))
    colnames(NewMatrix) <- NewColname
    NewMatrix[ ,1]<- data.matrix(rawd[, 1])
    NewMatrix[ ,2]<- data.matrix(rawd[, 2])
    NewMatrix[ ,3]<- data.matrix(rawd[, 3])
    NewMatrix[ ,4]<- data.matrix(rawd[, 4])
    NewMatrix[ ,5]<- data.matrix(rawd[, 5])
    NewMatrix[ ,6]<- data.matrix(rawd[, 6])
    NewMatrix[ ,7]<- data.matrix(rawd[, 7])
    NewMatrix[ ,8]<- data.matrix(rawd[, 8])
    NewMatrix[ ,9]<- data.matrix(rawd[, 9])
    NewMatrix[ ,10]<- data.matrix(rawd[, 10])
    NewMatrix[ ,11]<- 2^as.numeric(rawd[ , 10])
    NewMatrix[ ,12]<- data.matrix(rawd[, 11])
    NewMatrix[ ,13]<- data.matrix(rawd[, 12])
    NewMatrix[ ,14]<- data.matrix(rawd[, 13])
    NewMatrix[ ,15]<- data.matrix(rawd[, 14])
    data = data.frame(NewMatrix)
  } else {
    data <- RawData
  }
  return(data)
}

wb = createWorkbook()
for (i in 1:length(diffFiles)) {
  RawData = read.table(paste(path1, diffFiles[i], sep="/"), header = TRUE)
  data <- Myfunction(RawData)
  addWorksheet(wb, sheetName=diffFiles[i])
  writeDataTable(wb, i, x = data)
  setColWidths(wb, i, cols=8, widths = 14)
  
}





###################################################
### set working directory to output directory
###################################################
#setwd(outDir)
saveWorkbook(wb, "Task 42S.xlsx",overwrite = TRUE )







