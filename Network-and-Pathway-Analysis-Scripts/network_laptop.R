###########################################################################
# Subnetwork built using bum model of optimal and largest component of connected modules from the human PPI network (interactome)
# Zazil Villaueva-E. and Cory C. 
# Riley's lab -UMASS Boston Biology
############################################################################
library(BioNet)
library(DLBCL)
library(RCytoscape)

#We need to keep track of versions used, especially packages
sessionInfo() 
data(interactome) 

############# Functions ###############
reformatSym = function(table)
{
  # reformat Symbol column to match node labels in the graph
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

createModule = function(table)
{
  table$PValue = p.adjust(table$PValue, method = "bonferroni", n = nrow(table))
  pval = table$PValue # get p values from differential expression data
  pval = pval + (1 * 10 ^ -300) #applying correction factor
  names(pval) = table$Symbol #name pvalues with gene symbols
  
  subnet = subNetwork(table$Symbol, interactome) #Get subnetwork from interactome
  subnet = largestComp(subnet) #include connected nodes only
  subnet = rmSelfLoops(subnet) #remove self loops
  bum = fitBumModel(pval, plot = F) # fit Beta-uniform mixture model
  scores = scoreNodes(subnet, bum, fdr = 0.0001) # score nodes FDR .001 or which?
  logFC = as.numeric(table$logFC)
  names(logFC) = table$Symbol
  module = runFastHeinz(subnet, scores) #create optimal module based on scores
  FDR = scanFDR(bum, fdr=c(.0001, .00001), labels=names(pval))
  
  myList = list(module, pval, scores, logFC, FDR)
  return(myList)
}

saveModule = function(module, pval, scores, logFC, dir_out, filename) 
{
  #options(bitmapType ='cairo')
  dev.set()
  pathfile = paste(dir_out, filename, sep= "/")
  outfile = paste(pathfile, ".png", sep = "")
  png(file=outfile)
  plotModule(module, scores = scores, diff.expr = logFC)
  plot3dModule(module, diff.or.scores=logFC, red=c("positive"))
  saveNetwork(module, name=filename, file=outfile, type=c("table"))
  saveNetwork(module, name=filename, file=outfile, type=c("XGMML"))
  saveNetwork(module, name=filename, file=outfile, type=c("sif"))
  dev.off()  
}
########################################

dir="~"
inputFile="edger_table_.001_thresh100.txt"
dir_out="~"

table = read.delim(paste(dir, inputFile, sep = "/"))
# reformat the 'Symbol' column so that genes in the data can be matched with nodes in the network
table = reformatSym(table)
dim(table) #table with all pvalues from differential expression
max(table$logFC)
min(table$logFC)

Wholetable=table[table$PValue<1e-200,]
Whole = createModule(Wholetable)
Whole[[5]]
saveModule(module=Whole[[1]], pval=Whole[[2]], scores=Whole[[3]], logFC=Whole[[4]], dir_out=dir_out, filename="Whole")

UpregTable = table[table$logFC>0,]
#UpregTable = UpregTable[UpregTable$PValue<1e-200,]
max(UpregTable$logFC)
min(UpregTable$logFC)
dim(UpregTable)
Upreg = createModule(UpregTable)
saveModule(module=Upreg[[1]], pval=Upreg[[2]], scores=Upreg[[3]], logFC=Upreg[[4]], dir_out=dir_out, filename="Upreg")

DownregTable = table[table$logFC<0, ]
DownregTable = DownregTable[DownregTable$PValue<1e-200,]
max(DownregTable$logFC)
min(DownregTable$logFC)
dim(DownregTable)
Downreg = createModule(DownregTable)
saveModule(module=Downreg[[1]], pval=Downreg[[2]], scores=Downreg[[3]], logFC=Downreg[[4]], dir_out=dir_out, filename="Downreg")

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

