install.packages("gplots", dependencies=TRUE)
library("gplots")
install.packages("RColorBrewer")
library("RColorBrewer")

data <- read.csv("C:\\users\\travi_000\\Desktop\\cluster.csv", sep=",", comment.char="#")
data
rnames = data[,1]
mat_data = data.matrix(data[,2:ncol(data)])
rownames(mat_data)=rnames

my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

col_breaks = c(seq(-1,0,length=100),
               seq(0,0.8,length=100), 
               seq(0.8,1,length=100))  

png("C:\\Users\\travi_000\\Desktop\\cluster.csv", 
    width = 5*300,
    height = 5*300,
    res = 300,  
    pointsize = 8) 
    
heatmap.2(mat_data,
          cellnote = mat_data,  # same data set for cell labels
          main = "Correlation", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA")            # turn off column clustering

dev.off()               # close the PNG device
    
    
    
    
