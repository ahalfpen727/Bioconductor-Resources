### R code from vignette source 'STRINGdb.Rnw'

###################################################
### code chunk number 1: species (eval = FALSE)
###################################################
## get_STRING_species(version="10", species_name=NULL)


###################################################
### code chunk number 2: initialization
###################################################
library(STRINGdb)
string_db <- STRINGdb$new( version="10", species=9606,
                           score_threshold=0, input_directory="" )


###################################################
### code chunk number 3: help
###################################################
STRINGdb$methods()              # To list all the methods available.
STRINGdb$help("get_graph")      # To visualize their documentation.


###################################################
### code chunk number 4: load_data
###################################################
data(diff_exp_example1)
head(diff_exp_example1)


###################################################
### code chunk number 5: map
###################################################
example1_mapped <- string_db$map( diff_exp_example1, "gene", removeUnmappedRows = TRUE )


###################################################
### code chunk number 6: STRINGdb.Rnw:118-121
###################################################
options(SweaveHooks=list(fig=function()
par(mar=c(2.1, 0.1, 4.1, 2.1))))
#par(mar=c(1.1, 0.1, 4.1, 2.1))))


###################################################
### code chunk number 7: get_hits
###################################################
hits <- example1_mapped$STRING_id[1:200]


###################################################
### code chunk number 8: plot_network
###################################################
getOption("SweaveHooks")[["fig"]]()
string_db$plot_network( hits )


###################################################
### code chunk number 9: add_diff_exp_color
###################################################
# filter by p-value and add a color column
# (i.e. green down-regulated gened and red for up-regulated genes)
example1_mapped_pval05 <- string_db$add_diff_exp_color( subset(example1_mapped, pvalue<0.05),
                                                            logFcColStr="logFC" )


###################################################
### code chunk number 10: post_payload
###################################################
# post payload information to the STRING server
payload_id <- string_db$post_payload( example1_mapped_pval05$STRING_id,
                                        colors=example1_mapped_pval05$color )


###################################################
### code chunk number 11: plot_halo_network
###################################################
getOption("SweaveHooks")[["fig"]]()
# display a STRING network png with the "halo"
string_db$plot_network( hits, payload_id=payload_id )


###################################################
### code chunk number 12: STRINGdb.Rnw:183-185
###################################################
options(SweaveHooks=list(fig=function()
par(mar=c(2.1, 2.1, 4.1, 2.1))))


###################################################
### code chunk number 13: plot_ppi_enrichment
###################################################
getOption("SweaveHooks")[["fig"]]()
# plot the enrichment for the best 1000 genes
string_db$plot_ppi_enrichment( example1_mapped$STRING_id[1:1000], quiet=TRUE )


###################################################
### code chunk number 14: enrichment
###################################################
enrichmentGO <- string_db$get_enrichment( hits, category = "Process", methodMT = "fdr", iea = TRUE )
enrichmentKEGG <- string_db$get_enrichment( hits, category = "KEGG", methodMT = "fdr", iea = TRUE )
head(enrichmentGO, n=7)
head(enrichmentKEGG, n=7)


###################################################
### code chunk number 15: background (eval = FALSE)
###################################################
## backgroundV <- example1_mapped$STRING_id[1:2000]   # as an example, we use the first 2000 genes
## string_db$set_background(backgroundV)


###################################################
### code chunk number 16: new_background_inst (eval = FALSE)
###################################################
## string_db <- STRINGdb$new( score_threshold=0, backgroundV = backgroundV )


###################################################
### code chunk number 17: enrichmentHeatmap (eval = FALSE)
###################################################
## eh <- string_db$enrichment_heatmap( list( hits[1:100], hits[101:200]),
##                                     list("list1","list2"), title="My Lists" )


###################################################
### code chunk number 18: clustering1
###################################################
# get clusters
clustersList <- string_db$get_clusters(example1_mapped$STRING_id[1:600])


###################################################
### code chunk number 19: STRINGdb.Rnw:254-256
###################################################
options(SweaveHooks=list(fig=function()
par(mar=c(2.1, 0.1, 4.1, 2.1))))


###################################################
### code chunk number 20: clustering2
###################################################
getOption("SweaveHooks")[["fig"]]()
# plot first 4 clusters
par(mfrow=c(2,2))
for(i in seq(1:4)){
 string_db$plot_network(clustersList[[i]])
}


###################################################
### code chunk number 21: proteins
###################################################
string_proteins <- string_db$get_proteins()


###################################################
### code chunk number 22: atmtp
###################################################
tp53 = string_db$mp( "tp53" )
atm = string_db$mp( "atm" )


###################################################
### code chunk number 23: neighbors (eval = FALSE)
###################################################
## string_db$get_neighbors( c(tp53, atm) )


###################################################
### code chunk number 24: interactions
###################################################
string_db$get_interactions( c(tp53, atm) )


###################################################
### code chunk number 25: pubmedInteractions (eval = FALSE)
###################################################
## string_db$get_pubmed_interaction( tp53, atm )


###################################################
### code chunk number 26: homologs (eval = FALSE)
###################################################
## # get the reciprocal best hits of the following protein in all the STRING species
## string_db$get_homologs_besthits(tp53, symbets = TRUE)


###################################################
### code chunk number 27: homologs2 (eval = FALSE)
###################################################
## # get the homologs of the following two proteins in the mouse (i.e. species_id=10090)
## string_db$get_homologs(c(tp53, atm), target_species_id=10090, bitscore_threshold=60 )


###################################################
### code chunk number 28: benchmark1
###################################################
data(interactions_example)

interactions_benchmark = string_db$benchmark_ppi(interactions_example, pathwayType = "KEGG",
		max_homology_bitscore = 60, precision_window = 400, exclude_pathways = "blacklist")


###################################################
### code chunk number 29: STRINGdb.Rnw:391-393
###################################################
options(SweaveHooks=list(fig=function()
par(mar=c(4.1, 4.1, 4.1, 2.1))))


###################################################
### code chunk number 30: benchmark2
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(interactions_benchmark$precision, ylim=c(0,1), type="l", xlim=c(0,700),
	 xlab="interactions", ylab="precision")


###################################################
### code chunk number 31: benchmark3
###################################################
interactions_pathway_view = string_db$benchmark_ppi_pathway_view(interactions_benchmark, precision_threshold=0.2, pathwayType = "KEGG")
head(interactions_pathway_view)


