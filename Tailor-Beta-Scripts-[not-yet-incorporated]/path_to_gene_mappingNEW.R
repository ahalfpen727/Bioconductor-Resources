 Differential expression analysis
# Note importing BioC pkgs after dplyr requires explicitly using dplyr::select()
library(dplyr)
library(DESeq2)

# Which data do you want to use? Let's use the sailfish counts.
# browseURL("http://dx.doi.org/10.6084/m9.figshare.1601975")
# countDataURL = "http://files.figshare.com/2439061/GSE37704_featurecounts.csv"
countDataURL = "http://files.figshare.com/2600373/GSE37704_sailfish_genecounts.csv"

# Import countdata
countData = read.csv(countDataURL, row.names=1) %>% 
  dplyr::select(-length) %>% 
  as.matrix()

# Filter data where you only have 0 or 1 read count across all samples.
countData = countData[rowSums(countData)>1, ]
head(countData)

# Import metadata
colData = read.csv("http://files.figshare.com/2439060/GSE37704_metadata.csv", row.names=1)
colData
# Set up the DESeqDataSet Object and run the DESeq pipeline
dds = DESeqDataSetFromMatrix(countData=countData,
                              colData=colData,
                              design=~condition)
dds = DESeq(dds)
dds

res = results(dds, contrast=c("condition", "hoxa1_kd", "control_sirna"))
res = res[order(res$pvalue),]
summary(res)

library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)

##  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT" 
##  [5] "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"    
##  [9] "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"       
## [13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"    
## [17] "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
## [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
## [25] "UNIGENE"      "UNIPROT"
res$symbol = mapIds(org.Hs.eg.db,
                     keys=row.names(res), 
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez = mapIds(org.Hs.eg.db,
                     keys=row.names(res), 
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
res$name =   mapIds(org.Hs.eg.db,
                     keys=row.names(res), 
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")

head(res, 10)


library(pathview)
library(gage)
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)

 The gage() function requires a named vector of fold changes, where the names of the values are the Entrez gene IDs.

foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)
# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)

# Look at both up (greater), down (less), and statatistics.
lapply(keggres, head)

# Get the pathways
keggrespathways = data.frame(id=rownames(keggres$greater), 
							 keggres$greater) %>% 
	tbl_df() %>% 
  filter(row_number()<=5) %>% 
  .$id %>% 
  as.character()
keggrespathways

# Get the IDs.
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids

# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges, 
									  pathway.id=pid, species="hsa",
									  new.signature=FALSE)

# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggresids, 
			 function(pid) pathview(gene.data=foldchanges, 
			 					   pathway.id=pid, species="hsa"))


data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]

gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)

lapply(gobpres, head)

library(biomaRt)
library(dplyr)

fix_genes <- . %>% 
  tbl_df %>% 
  distinct %>% 
  rename(ensgene=ensembl_gene_id,
         entrez=entrezgene,
         symbol=external_gene_name,
         chr=chromosome_name,
         start=start_position,
         end=end_position,
         biotype=gene_biotype)

myattributes <- c("ensembl_gene_id",
                  "entrezgene",
                  "external_gene_name",
                  "chromosome_name",
                  "start_position",
                  "end_position",
                  "strand",
                  "gene_biotype",
                  "description")

# Human
grch38 <- useMart("ensembl") %>% 
  useDataset(mart=., dataset="hsapiens_gene_ensembl") %>% 
  getBM(mart=., attributes=myattributes) %>% 
  fix_genes

# Human grch37
grch37 <- useMart("ENSEMBL_MART_ENSEMBL", 
                  host="grch37.ensembl.org") %>% 
  useDataset(mart=., dataset="hsapiens_gene_ensembl") %>% 
  getBM(mart=., attributes=myattributes) %>% 
  fix_genes

# Mouse
grcm38 <- useMart("ensembl") %>% 
  useDataset(mart=., dataset="mmusculus_gene_ensembl") %>% 
  getBM(mart=., attributes=myattributes) %>% 
  fix_genes

# Rat
rnor6 <- useMart("ensembl") %>% 
  useDataset(mart=., dataset="rnorvegicus_gene_ensembl") %>% 
  getBM(mart=., attributes=myattributes) %>% 
  fix_genes

# Chicken
galgal5 <- useMart("ensembl") %>% 
  useDataset(mart=., dataset="ggallus_gene_ensembl") %>% 
  getBM(mart=., attributes=myattributes) %>% 
  fix_genes

# Fly
bdgp6 <- useMart("ensembl") %>% 
  useDataset(mart=., dataset="dmelanogaster_gene_ensembl") %>% 
  getBM(mart=., attributes=myattributes) %>% 
  fix_genes

# Worm
wbcel235 <- useMart("ensembl") %>% 
  useDataset(mart=., dataset="celegans_gene_ensembl") %>% 
  getBM(mart=., attributes=myattributes) %>% 
  fix_genes

fix_txps <- . %>% 
  tbl_df %>% 
  distinct %>% 
  rename(ensgene=ensembl_gene_id,
         enstxp=ensembl_transcript_id)

# Human build 38
grch38_gt <- useMart("ensembl") %>% 
  useDataset(mart=., dataset="hsapiens_gene_ensembl") %>% 
  getBM(mart=., attributes=c("ensembl_gene_id", "ensembl_transcript_id")) %>% 
  fix_txps

# Human build 37
grch37_gt <- useMart("ENSEMBL_MART_ENSEMBL", 
                     host="grch37.ensembl.org") %>% 
  useDataset(mart=., dataset="hsapiens_gene_ensembl") %>% 
  getBM(mart=., attributes=c("ensembl_gene_id", "ensembl_transcript_id")) %>% 
  fix_txps

# Mouse build 38
grcm38_gt <- useMart("ensembl") %>% 
  useDataset(mart=., dataset="mmusculus_gene_ensembl") %>% 
  getBM(mart=., attributes=c("ensembl_gene_id", "ensembl_transcript_id")) 

rm(fix_genes, fix_txps, myattributes)
devtools::use_data(grch38)
devtools::use_data(grch37)
devtools::use_data(grcm38)
devtools::use_data(rnor6)
devtools::use_data(galgal5)
devtools::use_data(bdgp6)
devtools::use_data(wbcel235)
devtools::use_data(grch38_gt)
devtools::use_data(grch37_gt)
devtools::use_data(grcm38_gt)


install.packages("devtools")
devtools::install_github("stephenturner/annotables")

It isnt necessary to load dplyr, but the tables are tbl_df and will print nicely if you have dplyr loaded.

library(dplyr)
library(annotables)

#Look at the human genes table (note the description column gets cut off because the table becomes too wide to print nicely):
grch38

#Look at the human genes-to-transcripts table:
grch38_gt


library(DESeq2)
library(airway)

data(airway)
airway <- DESeqDataSet(airway, design = ~cell + dex)
airway <- DESeq(airway)
res <- results(airway)

# tidy results with biobroom
library(biobroom)
res_tidy <- tidy.DESeqResults(res)
head(res_tidy)

res_tidy %>% 
  arrange(p.adjusted) %>% 
  head(20) %>% 
  inner_join(grch38, by=c("gene"="ensgene")) %>% 
  select(gene, estimate, p.adjusted, symbol, description) %>% 
  pander::pandoc.table(split.table=100, justify="lrrll", style="rmarkdown")

