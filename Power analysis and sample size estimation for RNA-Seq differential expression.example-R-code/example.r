
my_directory =  "C:\\Users\\tching\\Desktop\\sample size paper\\sampleSize_program"  #change this to your working directory
setwd(my_directory)
source("rs_simulations.r")


#Example 1: paired data
load("paired_example.Rdata")
estimate_params(rawdata=rawdata, condition=factor(treatment), pairing = factor(patient), designtype="paired") -> params
save(params, file="paired_params.Rdata")
RS_simulation(sims=5, params=params, budget=3000, designtype = "paired", nmax = 5, nmin = 2, program="DESeq2") -> results
plot(rownames(results),rowMeans(results, na.rm=T), main="DESeq2 simulations on paired example", xlab = "number of replicates", ylab = "Power")




#Example 2: one factor design
read.table("bottomly_count_table.txt", sep = "\t", header=T) -> rawdata
rownames(rawdata) <- rawdata$gene
rawdata <- as.matrix(rawdata[,-1])
condition = c(rep("A", 10),rep("B", 11))
estimate_params(rawdata=rawdata, condition=factor(condition), designtype="one factor") -> params
save(params, file="bottomly_params.Rdata")
RS_simulation(sims=5, params=params, budget=3000, designtype = "one factor", nmax = 5, nmin = 2, program="DESeq") -> results
plot(rownames(results),rowMeans(results, na.rm=T), main="DESeq simulations on Bottomly Dataset", xlab = "number of replicates", ylab = "Power")
