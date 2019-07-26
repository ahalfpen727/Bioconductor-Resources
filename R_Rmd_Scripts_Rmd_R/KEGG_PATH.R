###########################################################
###  R Script to Map Arabidpsis AGI Keys KEGG Pathways  ###
###########################################################
# Author: Thomas Girke, Mar 2006
# (A) Usage
# (A.1) Assign your locus IDs (AGI Keys) to vector
# (A.2) Load kegg_stat_fct as shown under demo (A.4)
# (A.3) Run kegg_stat_fct like this:
#	kegg_stat_fct(locus_ids, type=2) 
#		argument 'type': 
#			2 = Pathway
#			4 = PathwayClass
# (A.4) To demo what the script does, run it like this:
#	source("http://bioinfo.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/KEGG_ATH.txt")
# (B) Dependencies
# (B.1) Pathway category file "KEGG_ATH.annot" created from these pages:
#		http://www.genome.jp/kegg/KGML/KGML_v0.6/ath/index.html
# 		http://www.genome.jp/kegg/KGML/KGML_v0.6/ath/index2.html
# (B.2) Script can be used for other organisms after adjusting the regular expressions in part (1),
#     and creating the corresponding pathway category file (KEGG_ATH.annot)

# (C) Script
# (C.1) Import ATH KEGG annotations: the following commands are only necessary to update the info in kegg data frame (see below)
# Import pathway category file "KEGG_ATH.annot"
# kegg <- read.table("http://bioinfo.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/KEGG_ATH.annot", h=T, sep="\t")
# If an object "mylist" is present then remove it, since the following loop will append to it.  
# if(exists("mylist")) { rm(mylist) }
# for loop to import and reformat remote pathway xml files one-by-one by 'PathID' in 'kegg' (103 iterations)
# for(i in 1:length(kegg[,3])) {
# 	zzz <- readLines(paste("http://www.genome.jp/kegg/KGML/KGML_v0.6/ath/", kegg[i,3], ".xml", sep="")) # read lines into vector
# 	zzz <- zzz[grep("at\\dg\\d\\d\\d\\d\\d", as.character(zzz), ignore.case = T, perl=T)] # obtain fields with AGI keys
# 	zzz <- zzz[-grep("link=", as.character(zzz), ignore.case = T, perl=T)] # remove duplications
# 	zzz <- gsub(".*name=\\\"|\".*", "", as.character(zzz), ignore.case = T, perl=T) # remove xml tags
# 	zzz <- gsub("ath:", "", as.character(zzz), ignore.case = T, perl=T) # obtain clean AGI keys
# 	zzz <- unlist(strsplit(zzz, " ")) # split fields with several keys; strsplit generates list of vectors and unlist converts it into one vector
# 	zzz <- unique(zzz) # makes AGIs unique within each PathID
# 	if(!exists("mylist")) # if "mylist" does not exist then create it in next line, otherwise append to "mylist" in second line  
# 		{ mylist <- list(zzz) } else 
# 		{ mylist <- c(mylist, list(zzz)) } 
# 	names(mylist)[i] <- as.vector(kegg[i,3]) # assign names by 'PathID'
# 	cat(i, as.vector(kegg[i,3]), "\n", sep = " ") # print iteration number and 'PathID' to monitor loop status
# }
# keggat <- data.frame(unlist(mylist)) # transform list of vectors into data frame, tracking of one-to-many relationship is resolved by appending counts to PathIDs
# keggat <- data.frame(PathID=gsub("(ath\\d\\d\\d\\d\\d).*", "\\1", row.names(keggat), perl=T), AGI=keggat[,1]) # remove from PathIDs appended counts
# keggat <- data.frame(PathID=keggat$PathID, AGI=gsub("([A-Z])", "\\U\\1", as.character(keggat$AGI), perl=T, ignore.case=T)) # make AGI keys uniform by setting them to upper case
# kegg <- merge(keggat, kegg, by.x="PathID", by.y="PathID", all.x=T) # merge data frames 'keggat' and 'kegg' (KEGG_ATH.annot) by PathIDs

# If no update is required, then read the Arab KEGG mappings from this file:
kegg <- read.table("http://bioinfo.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/kegg_ath", h=T, sep="\t")

# (C.2) Define function: "kegg_stat_fct"
kegg_stat_fct <- function(locus_ids, type=2) { # argument 'type': 2 = Pathway, 4 = PathwayClass (see kegg)
	my_select <- merge(kegg, as.data.frame(locus_ids), by.x="AGI", by.y=1, all.y=T) # obtain kegg data frame that cotains only rows for AGI test keys (locus_ids)
	x <- as.character(my_select[,type]) # create vector with Pathway (2) or PathwayClass (4) strings in my_select
	x[is.na(x)] <- "no_annot" # fill NA fields with "no_annot"
	my_select <- data.frame(count_col=x, my_select[,-type]) # combine x and my_select
	kegg_stat <- data.frame(table(my_select$count_col)) # obtain counts for count_col 
	kegg_stat <- data.frame(kegg_stat, Percent=round(kegg_stat$Freq/sum(kegg_stat$Freq)*100,2)) # append percent calculation for counts
	names(kegg_stat)[1] <- "Cat" # give frist column a descriptive title
	if(type==2) # if type==2 then merge kegg_stat with kegg data frame on column 1, if type==4 then merge kegg_stat with kegg data frame on column 4
		{ kegg_stat <- merge(kegg_stat, kegg[!duplicated(kegg[,type-1]),], by.x="Cat", by.y=type-1, all.x=T) } else
		{ kegg_stat <- merge(kegg_stat, kegg[!duplicated(kegg[,type]),], by.x="Cat", by.y=type, all.x=T) }
	kegg_stat
}

# (C.3) Use function: "kegg_stat_fct"
cat("\n", "You have imported the function 'kegg_stat_fct'.", "\n\n", "Here is a short demo for AGI keys:", "\n") # print usage of function
AGI_vec <- c("AT1G23800", "AT5G63620", "AT1G30120", "AT2G07050", "AT2G22830", "AT2G38700", "AT1G29910")
print(AGI_vec)
kegg_stat <- kegg_stat_fct(locus_ids=AGI_vec, type=4)
print(kegg_stat)
pie(kegg_stat[kegg_stat[,2]>0 ,3], labels=as.vector(kegg_stat[kegg_stat[,2]>0 ,1]), main="KEGG") # plot final result

cat("\n", "Usage:", "\n", "\t kegg_stat <-  kegg_stat_fct(locus_ids=AGI_vec, type=4)", "\n", "\t Arguments: 'locus_ids', vector of locus keys; 'type=2', Pathway; 'type=4', PathwayClass", "\n", "\t pie(kegg_stat[kegg_stat[,2]>0 ,3], labels=as.vector(kegg_stat[kegg_stat[,2]>0 ,1]), main=\"KEGG\")", "\n")
