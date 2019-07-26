library(xlsx)
mydata1 = read.csv("C:/Users/Arpa.Samadder001/Desktop/new/CXCL12-over.csv", header = T)
mydata1
mydata2 = read.csv("C:/Users/Arpa.Samadder001/Desktop/new/TGFB-over.CSV", header = T)
mydata2

#This line gives an output file that combines the 2 tables with same number of rows

myfulldata = cbind(mydata1, mydata2[1:1046,]) 

write.table(myfulldata, "C:/Users/Arpa.Samadder001/Desktop/myfulldata.txt", sep="\t")

#Below the code will give the output file that combines the 2 tables with an unequal length:
#CXCL12-over-control and TGFB-over-control have unequal length.
#Here's a simple function that will cbind data frames of uneven length 
#and automatically pad the shorter ones with NAs:

cbindPad <- function(...){
 args <- list(...)
 n <- sapply(args,nrow)
 mx <- max(n)
 pad <- function(x, mx){
   if (nrow(x) < mx){
     nms <- colnames(x)
     padTemp <- matrix(NA, mx - nrow(x), ncol(x))
     colnames(padTemp) <- nms
     if (ncol(x)==0) {
       return(padTemp)
     } else {
       return(rbind(x,padTemp))
     }
   }
   else{
     return(x)
   }
 }
 rs <- lapply(args,pad,mx)
 return(do.call(cbind,rs))
}

myfulldatafinal = cbindPad(mydata1, mydata2)

write.table(myfulldatafinal, "C:/Users/Arpa.Samadder001/Desktop/myfulldatafinal.txt", sep="\t")

#write.xlsx(x = sample.dataframe, file = "test.excelfile.xlsx",
#sheetName = "TestSheet", row.names = FALSE)

write.xlsx(x = myfulldatafinal, file = "myfulldatafinal.xlsx", 
           sheetName="TestSheet", row.names = FALSE)
