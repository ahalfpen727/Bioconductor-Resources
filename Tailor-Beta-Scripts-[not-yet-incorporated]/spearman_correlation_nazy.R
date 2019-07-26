

o = order(nanostring12$ID)
nanostring12 = nanostring12[o,]

fc12 = NULL
fc12 = data.frame(cbind(as.character(nanostring12$ID), nanostring12$FC))
colnames(fc12) = c("Accession", "nanostring.FC")
fc12$Accession = gsub('\\..*$','', as.character(fc12$Accession))

o1  = order(table12$Accession)
table12 = table12[o1,]
table12[as.character(table12$Accession) %in% as.character(fc12$Accession),]$FC
fc12 = fc12[as.character(fc12$Accession) %in% as.character(table12$Accession),]
fc12 = cbind(fc12, table12[as.character(table12$Accession) %in% as.character(fc12$Accession),]$FC)
colnames(fc12) = c("Accession", "nanostring.FC", "hiseq.FC")
test12 = cor.test(as.numeric(fc12$nanostring.FC), as.numeric(fc12$hiseq.FC), method = "spearman")
test12

as.numeric(as.vector(fc12$nanostring.FC))

with(fc12, plot(as.numeric(as.vector(nanostring.FC)), as.numeric(as.vector(hiseq.FC)),
                xlab = "Nanostring FC", ylab = "HiSeq FC", 
                main = "SAMP6 HiSeq and Nanostring Correlation"))
legend(x = 1.8, y = 5.7, legend = "Spearman corr: 0.53", bty = "n")

test12.pearson = cor.test(as.numeric(fc12$nanostring.FC), as.numeric(fc12$hiseq.FC), method = "pearson")
test12.pearson

o = order(nanostring34$ID)
nanostring34 = nanostring34[o,]

fc34 = NULL
fc34 = data.frame(cbind(as.character(nanostring34$ID), nanostring34$FC))
colnames(fc34) = c("Accession", "nanostring.FC")
fc34$Accession = gsub('\\..*$','', as.character(fc34$Accession))
o1  = order(table34$Accession)
table34 = table34[o1,]
table34[as.character(table34$Accession) %in% as.character(fc34$Accession),]$FC
fc34 = fc34[as.character(fc34$Accession) %in% as.character(table34$Accession),]
fc34 = cbind(fc34, table34[as.character(table34$Accession) %in% as.character(fc34$Accession),]$FC)
colnames(fc34) = c("Accession", "nanostring.FC", "hiseq.FC")
test = cor.test(as.numeric(fc34$nanostring.FC), as.numeric(fc34$hiseq.FC), method = "spearman")
with(fc34, plot(as.numeric(as.vector(nanostring.FC)), as.numeric(as.vector(hiseq.FC)),
                xlab = "Nanostring FC", ylab = "HiSeq FC", 
                main = "AKR/J HiSeq and Nanostring Correlation"))
legend(x = 1.25, y = 3.5, legend = "Spearman corr: 0.46", bty = "n")

write.table(fc12, "~/Documents/UMB/Riley_Lab/HiSeq data/nazy/samp6_nanostring_vs_hiseq_fc.txt")
write.table(fc34, "~/Documents/UMB/Riley_Lab/HiSeq data/nazy/akrj_nanostring_vs_hiseq_fc.txt")
