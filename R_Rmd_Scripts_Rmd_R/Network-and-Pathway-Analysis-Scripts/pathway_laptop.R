###########################################################################
# Reconstruction of a signaling regulatory pathway from RNAseq expression data
# Zazil Villaueva-Esperon
# Riley's lab UMASS Boston Biology
############################################################################
#source("http://bioconductor.org/biocLite.R")
biocLite("pathview")
# rpackage.dir="/project/umb_triley/Rpackages/"
# install.packages(rpackage.dir, repos = NULL, type = "source")
library(pathview)

data(paths.hsa) # human pathways IDs from KEGG
data(gse16873.d)
data(demo.paths)
#data(gene.idtype.list) # to map gene IDs to KEGG accesion number 

sessionInfo() # Keep track of versions used, especially packages
dir_in="~"
dir_out="~"

createPathway=function(file, path, species){
  #read file
  table = read.delim(paste(dir, file, sep = "/"), row.names=8) # use gene id from column # in table to map Kegg gene ID
  pvtable = data.frame(table$logFC)
  rownames(pvtable)=rownames(table)
  
  pv.out<- sapply(path, function(x) pathview(gene.data=pvtable, pathway.id = x,species = species, kegg.native = TRUE)) 
  #head(pv.out$plot.data.gene)
  # Graphviz view
  pv.out<- sapply(path, function(x) pathview(gene.data=pvtable, pathway.id = x, species = species, kegg.native = FALSE, sign.pos="bottomleft", same.layer=F)) 
  #head(pv.out$plot.data.gene)
  return(pvtable)
}

# pathway from human genome (hsa) using KEGG ID 04350, common name TGFB1, and diff exp column for logFC in table
# "MAPK"
# Max logFC is 2.57 is 2167 FABP4 NM_001442
# Min logFC is gene 1311 COMP -6.203
# Kegg view with Kegg ID or common name for pathway

path=c("04350", "map04010")
pvtable=createPathway("cxcl12_tgfb_p001_thresh100.txt",path,"hsa")
head(pvtable)

******* Example frm RNAseq ***********
  library(gage)
 data(gse16873)
 cn <- colnames(gse16873)
 hn <- grep('HN',cn, ignore.case =TRUE)
 dcis <- grep('DCIS',cn, ignore.case =TRUE)
 data(kegg.gs)
 #pathway analysis using gage
 gse16873.kegg.p <- gage(gse16873, gsets = kegg.gs, ref = hn, samp = dcis)
 #prepare the differential expression data
 gse16873.d <- gagePrep(gse16873, ref = hn, samp = dcis)
 #equivalently, you can do simple subtraction for paired samples
 gse16873.d <- gse16873[,dcis]-gse16873[,hn]
 #select significant pathways and extract their IDs
 sel <- gse16873.kegg.p$greater[, "q.val"] < 0.1 & !is.na(gse16873.kegg.p$greater[,"q.val"])
 path.ids <- rownames(gse16873.kegg.p$greater)[sel]
 path.ids2 <- substr(path.ids[c(1, 2, 7)], 1, 8)
 #pathview visualization
 pv.out.list <- sapply(path.ids2, function(pid) pathview(gene.data = gse16873.d[,1:2], pathway.id = pid, species = "hsa"))
