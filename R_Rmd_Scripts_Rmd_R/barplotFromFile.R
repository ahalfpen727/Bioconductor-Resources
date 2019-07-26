
setwd("/home/mesude/Desktop")
barcodes<-read.table(file="Barcode.txt",header=T,sep="\t")
barcodes
attach(barcodes)
plot(BarcodeSequence)
barplot(BarcodeSequence)

