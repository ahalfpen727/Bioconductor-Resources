## ---- fig.keep='last', fig.width=10, fig.height=6.5, message=FALSE, warning=FALSE, results='hide', tidy=TRUE, error=TRUE----
# Load the GenVisR package
library("GenVisR")
set.seed(426)

# Plot with the MAF file type specified (default)
# The mainRecurCutoff parameter is described in the next section
waterfall(brcaMAF, fileType="MAF", mainRecurCutoff=.05)

## ---- eval=FALSE, error=TRUE---------------------------------------------
#  # read in a file from the genome modeling system
#  file <- read.delim("file.anno.tsv")
#  
#  # Plot the variant information via waterfall
#  waterfall(file, fileType="MGI")

## ---- fig.keep='all', fig.width=10, fig.height=7, message=FALSE, warning=FALSE, error=TRUE, results='hide', tidy=TRUE, fig.show='hold', fig.cap="In cell e/e (second row/first column) two variants are present \'z\' and \'x\'. In the first plot variant \'z\' is considered more deleterious (top panel), In the second plot variant \'x\' is considered more deleterious (bottom panel).", out.width="50%", out.height="50%"----
# make sure seed is set to 426 to reproduce!
set.seed(426)

# Create a data frame of random elements to plot
inputData <- data.frame("sample"=sample(letters[1:5], 20, replace=TRUE), "gene"=sample(letters[1:5], 20, replace=TRUE), variant_class=sample(c("x", "y", "z"), 20, replace=TRUE))

# choose the most deleterious to plot with y being defined as the most deleterious
most_deleterious <- c("y", "z", "x")

# plot the data with waterfall using the "Custom" parameter
waterfall(inputData, fileType="Custom", variant_class_order = most_deleterious, mainXlabel = TRUE)

# change the most deleterious order
waterfall(inputData, fileType="Custom", variant_class_order = rev(most_deleterious), mainXlabel = TRUE)

## ----kable 1, echo=FALSE, error=TRUE, fig.cap="Hierarchy of variant sub types from most to least deleterious."----
library(knitr)
MGI <- c("nonsense", "frame_shift_del",
         "frame_shift_ins", "splice_site_del",
         "splice_site_ins", "splice_site",
         "nonstop", "in_frame_del", "in_frame_ins",
         "missense", "splice_region_del",
         "splice_region_ins", "splice_region",
         "5_prime_flanking_region",
         "3_prime_flanking_region",
         "3_prime_untranslated_region",
         "5_prime_untranslated_region", "rna",
         "intronic", "silent")
MAF <- c("Nonsense_Mutation", "Frame_Shift_Ins",
         "Frame_Shift_Del", "Translation_Start_Site",
         "Splice_Site", "Nonstop_Mutation",
         "In_Frame_Ins", "In_Frame_Del",
         "Missense_Mutation", "5\'Flank",
         "3\'Flank", "5\'UTR", "3\'UTR", "RNA", "Intron",
         "IGR", "Silent", "Targeted_Region", "", "")

kable(as.data.frame(cbind(MAF, MGI)))

## ---- fig.keep='last', fig.width=10, fig.height=6.5, message=FALSE, warning=FALSE, results='hide', tidy=TRUE, error=TRUE----
# Plot the genes with mutatations in >= 20% of samples
waterfall(brcaMAF, fileType="MAF", mainRecurCutoff=.2)

## ---- fig.keep='last', fig.width=10, fig.height=6.5, message=FALSE, warning=FALSE, results='hide', tidy=TRUE, error=TRUE----
# Define specific genes to plot
genes_to_plot <- c("ERBB2", "MAPK1", "CDKN1B", "PIK3CA")

# Plot the genes defined above
waterfall(brcaMAF, plotGenes=genes_to_plot)

## ---- fig.keep='last', fig.width=10, fig.height=6.5, message=FALSE, warning=FALSE, results='hide', tidy=TRUE, error=TRUE----
# Define specific genes to plot
samples_to_plot <- c("TCGA-A1-A0SO-01A-22D-A099-09", "TCGA-A2-A0EU-01A-22W-A071-09", "TCGA-A2-A0ER-01A-21W-A050-09", "TCGA-A1-A0SI-01A-11D-A142-09", "TCGA-A2-A0D0-01A-11W-A019-09")

# Plot the samples defined above
waterfall(brcaMAF, plotSamples=samples_to_plot, mainRecurCutoff=.25)

## ---- fig.keep='last', fig.width=10, fig.height=6.5, message=FALSE, warning=FALSE, results='hide', tidy=TRUE, error=TRUE----

# plotting all genes with a mutation recurrence above 5%, limit to plot only the top 25 and remove silent mutations
waterfall(brcaMAF, mainRecurCutoff=.05, maxGenes=25, rmvSilent=TRUE)

## ---- fig.keep='last', fig.width=10, fig.height=6.5, message=FALSE, warning=FALSE, error=TRUE, results='hide', tidy=TRUE, fig.cap="Altering the coverage space dramatically affects the mutation burden calculation"----

# Alter the coverage space to whole genome space
waterfall(brcaMAF, mainRecurCutoff=.05, maxGenes=25, coverageSpace=3200000000)

## ---- fig.keep='last', fig.width=10, fig.height=6.5, message=FALSE, warning=FALSE, error=TRUE, results='hide', tidy=TRUE----
# Create a data frame specifying the mutation burden for each sample
tumor_sample <- unique(brcaMAF$Tumor_Sample_Barcode)
mutation_burden <- sample(1:10, length(tumor_sample), replace=TRUE)
mutation_rate <- data.frame(sample=tumor_sample, mut_burden=mutation_burden)

# Alter the coverage space to whole genome space
waterfall(brcaMAF, mutBurden=mutation_rate, mainRecurCutoff=.05, maxGenes=25)

## ---- fig.keep='last', fig.width=10, fig.height=6.5, message=FALSE, warning=FALSE, error=TRUE, results='hide', tidy=TRUE----
# Turn off plotting of the mutation burden subplot
waterfall(brcaMAF, plotMutBurden=FALSE, mainRecurCutoff=.05, maxGenes=25)

## ---- fig.keep='last', fig.width=14, fig.height=10, message=FALSE, warning=FALSE, error=TRUE, results='hide', tidy=TRUE----
# Create clinical data
subtype <- c('lumA', 'lumB', 'her2', 'basal', 'normal')
subtype <- sample(subtype, 50, replace=TRUE)
age <- c('20-30', '31-50', '51-60', '61+')
age <- sample(age, 50, replace=TRUE)
sample <- as.character(unique(brcaMAF$Tumor_Sample_Barcode))
clinical <- as.data.frame(cbind(sample, subtype, age))

# Melt the clinical data into "long" format.
library(reshape2)
clinical <- melt(clinical, id.vars=c('sample'))

# create the waterfall plot with the corresponding clinical data
waterfall(brcaMAF, clinDat=clinical, mainRecurCutoff=.05, maxGenes=25)

## ---- fig.keep='last', fig.width=12, fig.height=7.5, message=FALSE, warning=FALSE, error=TRUE, results='hide', tidy=TRUE----

# Create clinical data
subtype <- c('lumA', 'lumB', 'her2', 'basal', 'normal')
subtype <- sample(subtype, 50, replace=TRUE)
age <- c('20-30', '31-50', '51-60', '61+')
age <- sample(age, 50, replace=TRUE)
sample <- as.character(unique(brcaMAF$Tumor_Sample_Barcode))
clinical <- as.data.frame(cbind(sample, subtype, age))

# Melt the clinical data into "long" format.
library(reshape2)
clinical <- melt(clinical, id.vars=c('sample'))

# create the waterfall plot altering various aesthetics in the clinical data
waterfall(brcaMAF, clinDat=clinical,
          clinVarCol=c('lumA'='blue4', 'lumB'='deepskyblue', 
                            'her2'='hotpink2', 'basal'='firebrick2',
                            'normal'='green4', '20-30'='#ddd1e7',
                            '31-50'='#bba3d0', '51-60'='#9975b9',
                            '61+'='#7647a2'), 
          mainRecurCutoff=.05, maxGenes=25,
          clinLegCol=2,
          clinVarOrder=c('lumA', 'lumB', 'her2', 'basal', 'normal',
                         '20-30', '31-50', '51-60', '61+'))


## ---- fig.keep='last', fig.width=14, fig.height=6.5, message=FALSE, warning=FALSE, error=TRUE, results='hide', tidy=TRUE----
# Use the chromosome column in brcaMAF to label cells
waterfall(brcaMAF, mainRecurCutoff=.05, maxGenes=10, mainLabelCol="Chromosome")

## ---- fig.keep='last', fig.width=14, fig.height=8.5, message=FALSE, warning=FALSE, results='hide', error=TRUE, tidy=TRUE----
# Use the amino_acid change column in brcaMAF to label cells 
waterfall(brcaMAF, mainRecurCutoff=.05, maxGenes=10, mainLabelCol="amino_acid_change_WU", mainLabelAngle=90, mainLabelSize=3)

## ---- fig.keep='last', fig.width=14, fig.height=8.5, message=FALSE, warning=FALSE, results='hide', error=TRUE, tidy=TRUE----
# Label the x-axis 
waterfall(brcaMAF, mainRecurCutoff=.05, maxGenes=10, mainXlabel=TRUE)

## ---- fig.keep='last', fig.width=14, fig.height=8.5, message=FALSE, warning=FALSE, results='hide', tidy=TRUE, error=TRUE----
# Drop unused mutation types from the legend 
waterfall(brcaMAF, mainRecurCutoff=.05, maxGenes=10, mainDropMut=TRUE)

## ---- fig.keep='last', fig.width=14, fig.height=8.5, message=FALSE, warning=FALSE, results='hide', tidy=TRUE, error=TRUE----
# Increase the gene label size
waterfall(brcaMAF, mainRecurCutoff=.05, maxGenes=10, mainDropMut=TRUE, main_geneLabSize=14)

## ---- fig.keep='last', fig.width=14, fig.height=8.5, message=FALSE, warning=FALSE, results='hide', tidy=TRUE, error=TRUE----
# make a custom colour pallete
custom_pallete <- c("#A069C7", "#9CD05B", "#C46839", "#97BDBD", "#513C4D", "#6B7644", "#C6587F")

# provide a custom colour pallete
waterfall(brcaMAF, mainRecurCutoff=.05, maxGenes=10, mainDropMut=TRUE, mainPalette=custom_pallete)

## ---- fig.keep='last', fig.width=14, fig.height=8.5, message=FALSE, warning=FALSE, results='hide', tidy=TRUE, fig.cap="As can be seen we have changed the theme via ggplot2 however this has overwritten the previously defined theme which had suppressed the x-axis labels and removed the rotation.", error=TRUE----
# load ggplot2
library(ggplot2)

# suppress the y axis labels in the mutation burden plot
mut_burden_layer <- theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())

# change the ggplot theme back to default in the main plot
main_layer <- theme_grey()

# Run waterfall with the new layer
waterfall(brcaMAF, mainRecurCutoff=.05, maxGenes=10, mainDropMut=TRUE, mainLayer=main_layer, mutBurdenLayer=mut_burden_layer)

## ---- fig.keep='last', fig.width=14, fig.height=10, message=FALSE, warning=FALSE, results='hide', tidy=TRUE, error=TRUE----
# Create clinical data
subtype <- c('lumA', 'lumB', 'her2', 'basal', 'normal')
subtype <- sample(subtype, 50, replace=TRUE)
age <- c('20-30', '31-50', '51-60', '61+')
age <- sample(age, 50, replace=TRUE)
sample <- as.character(unique(brcaMAF$Tumor_Sample_Barcode))
clinical <- as.data.frame(cbind(sample, subtype, age))

# Melt the clinical data into "long" format.
library(reshape2)
clinical <- melt(clinical, id.vars=c('sample'))

# Obtain a sample order corresponding to the clinical data
new_samp_order <- as.character(unique(clinical[order(clinical$variable, clinical$value),]$sample))

# create the waterfall plot with the corresponding clinical data
waterfall(brcaMAF, clinDat=clinical, mainRecurCutoff=.05, maxGenes=25, sampOrder=new_samp_order)

## ---- fig.keep='last', fig.width=14, fig.height=8.5, message=FALSE, warning=FALSE, results='hide', tidy=TRUE, error=TRUE----
# Define a custom gene order
new_gene_order <- c("MUC16", "MUC17", "MUC12", "RYR2", "PIK3CA", "TP53", "USH2A", "MLL3", "TTN", "LRP2")

# Increase the gene label size
waterfall(brcaMAF, mainRecurCutoff=.05, maxGenes=10, geneOrder=new_gene_order)

## ---- eval=FALSE, error=TRUE---------------------------------------------
#  # Save a GenVisR plot
#  pdf(file="myplot.pdf", height=10, width=15)
#  waterfall(brcaMAF, mainRecurCutoff=.05, maxGenes=10)
#  dev.off()

## ---- fig.keep='last', fig.width=5, fig.height=5, message=FALSE, warning=FALSE, error=TRUE, results='hide', tidy=TRUE, fig.cap="A small graphics device may cause grob collisions, here the device size is 5 by 5"----
# A GenVisR plot on a small graphics device
waterfall(brcaMAF, mainRecurCutoff=.05, maxGenes=10)

## ---- fig.keep='last', fig.width=10, fig.height=7, message=FALSE, warning=FALSE, results='hide', error=TRUE, tidy=TRUE, fig.cap="Increasing the size of the graphic device will alleviate this issue, here the device size is 10 by 7"----
# A GenVisR plot on a small graphics device
waterfall(brcaMAF, mainRecurCutoff=.05, maxGenes=10)

