###########################################################################
# Subnetwork built using bum model of optimal and largest component of connected modules from the human PPI network (interactome)
# Zazil Villaueva-E. and Cory C. 
# Riley's lab -UMASS Boston Biology
############################################################################

# installRpkgs.sh installs all packages needed. This script needs to be run once prior to running this Rscript
rpackage.dir="/project/umb_triley/Rpackages/"

"Starting to load libraries..."
rpackage.dir

library(BioNet)
library(DLBCL,lib.loc= rpackage.dir)
library(RCytoscape)
#library(xlsx)

#We need to keep track of versions used, especially packages
#sessionInfo() 
"Loading interactome..."
data(interactome) 

"Finished loading all libraries..."
args <- commandArgs(TRUE)
args
dir=args[1]
file=args[2]
dir_out=args[3]

table = read.delim(paste(dir, file, sep = "/"))

"Finished reading input file..."

# reformat Symbol column to match node labels in the graph
reformatSym = function(table)
{
  table$Symbol = as.character(table$Symbol)
  table$ID = as.character(table$ID)
  for (node1 in interactome@nodes)
  {
    sym = interactome@nodeData@data[node1][[1]]$geneSymbol
    if (!is.null(table["Symbol"][table["Symbol"] == sym]))
    {
      table["Symbol"][table["Symbol"] == sym] = node1
    }
  }
  return(table)
}

# reformat the 'Symbol' column so that genes in the data can be matched with nodes in the network
table = reformatSym(table)

"dimensions of table..."
dim(table)

# perform Bonferroni correction on p values, which are too low for BUM algorithm
"defining table$PValue..."
#table$PValue = as.numeric(levels(table$PValue))[table$PValue]
"Bonferroni correction..."
table$PValue = p.adjust(table$PValue, method = "bonferroni", n = nrow(table))

# get the subnetwork
"Get subnetwork from interactome..."
subnet = subNetwork(table$Symbol, interactome)
"largestComponent..."
subnet = largestComp(subnet) #Function added to identify modules as connected subgraphs only
"remove self loops..."
subnet = rmSelfLoops(subnet)

"Subnet..."
subnet

# get p values from differential expression data
# p-values in the diff exp dataset were way to small for model to be created.
# Therefore, a correction factor was added to pump these values up
"Applying correction factor..."
pval = table$PValue
pval = pval + (1 * 10 ^ -300)
names(pval) = table$Symbol

# fit Beta-uniform mixture model
"BUM..."
bum = fitBumModel(pval, plot = F)

# score nodes based on model
"Scores..."
scores = scoreNodes(subnet, bum, fdr = 0.001)

# create module using runFastHeinz to calculate the max. scoring subnetwork (or optimal solution)
#logFC = as.numeric(levels(table$logFC))[table$logFC]
logFC=as.numeric(table$logFC)
names(logFC) = table$Symbol
"Module..."
module = runFastHeinz(subnet, scores) 

# display subnetwork
options(bitmapType='cairo')
filename = paste(dir_out, file,sep= "/")
outfile = paste(filename, ".png", sep = "")
png(file=outfile)
plotModule(module, scores = scores, diff.expr = logFC)
dev.off()

################### RCytoscape #######################################################
# RCytoscape requires a port to be listening using rpc protocol. Default is port 9000.
# Need to see what happens running the R script in the cluster and Cytoscape in my laptop.
# Most likely host will need to be changed to my IP address.
#######################################################################################
# 3,104 red (positive) nodes    # 6 green (negative) nodes    #0 nodes are white
#sum(scores[nodes(module)]>0)   #sum(scores[nodes(module)]<0) #sum(scores[nodes(module)]==0)

#g<-new('graphNEL', edgemode='directed')
#cw<-new.CytoscapeWindow(title="Jose's data", graph=g, overwriteWindow=TRUE, host = "localhost", rpcPort = 9000, )

#g<-initNodeAttribute (module, attribute.name='geneID', attribute.type='numeric', default.value=0)
#g<-initNodeAttribute (g, 'geneSymbol', 'char', 'x')
#g<-initNodeAttribute (g, 'score', 'numeric', scores)
#g<-initNodeAttribute(g,'lfc','numeric',logFC)

##g<-initEdgeAttribute(g,"weight",'integer',0)
##g<-initEdgeAttribute (g, "similarity", 'numeric', 0)
#nodeData(g,table$Symbol,"geneID")<-table$ID
#nodeData(g,table$Symbol,"geneSymbol")<-table$Symbol
#nodeData(g,table$Symbol,"score")<-scores
#nodeData(g,table$Symbol,"lfc")<-logFC
#cw=setGraph(cw,g)

#set node colors based on logFC
#green when logFC is negative, red is positive, white is zero and shades are in between
#setNodeColorRule(cw,"lfc", c (-3.0, 0.0, 3.0), c('#00AA00', '#00FF00', '#FFFFFF', '#FF0000', '#AA0000'), mode='interpolate')
#displayGraph(cw)
#redraw(cw)
#layoutNetwork(cw, hlp[18]) #default fruchterman.reingold layout is 18

# Apply values to some of the properties and plot the layout
#cy <- CytoscapeConnection()
#hlp <-getLayoutNames(cy) 

#print(noa.names(getGraph(cw))) #what data attributes are defined? noa=node attributes
##print(noa(getGraph(cw),'geneID'))
##print(noa(getGraph(cw),'geneSymbol'))
##print(noa(getGraph(cw),'score'))
##print(noa(getGraph(cw),'lfc'))

#saveNetwork(g, "Jose_RCytoscape", 'XGMML')  #doesn't work, so I just saved it directly from Cytoscape
############################################################

