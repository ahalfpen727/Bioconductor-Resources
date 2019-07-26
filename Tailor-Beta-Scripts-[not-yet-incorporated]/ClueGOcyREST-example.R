#
# ClueGO cyREST Example Workflow
#

library(RJSONIO)
library(httr)


text.to.data.frame <- function(table.text) {
    table <- NULL
    rows <- unlist(strsplit(result.table.text, split="\n"))
    header <- t(unlist(strsplit(rows[1], split="\t")))
    for(i in 2:length(rows)) {
        if(is.null(table)) {
            table <- t(unlist(strsplit(rows[i], split="\t")))
        } else {
            table <- rbind(table,t(unlist(strsplit(rows[i], split="\t"))))
        }
    }
    table <- as.data.frame(table)
    names(table) <- header
    return(table)
}

#Basic settings for cyREST

home.folder = "/home/gabi" # Linux / MacOSX
cluego.home.folder = paste(home.folder,"ClueGOConfiguration/v2.5.0",sep="/") # Linux / MacOSX

# home.folder = "C:\Users\UserName" # Windows
# cluego.home.folder = paste(home.folder,"ClueGOConfiguration\v2.5.0",sep="/") # Windows

port.number = 1234
host.address <- "localhost" # "192.168.0.20" #

cytoscape.base.url = paste("http://",host.address,":", toString(port.number), "/v1", sep="")
cluego.base.url = paste(cytoscape.base.url,"apps/cluego/cluego-manager", sep="/")

# list all available organisms
response <- GET(paste(cluego.base.url,"organisms","get-all-installed-organisms",sep="/"))
print(content(response))

# set organism to 'Homo Sapiens'
organism.name = "Homo Sapiens"
response <- PUT(url=paste(cluego.base.url,"organisms","set-organism",URLencode(organism.name),sep="/"), encode = "json")
#print(response)

# set gene list for cluster 1
cluster1 = "1"
file.location = paste(cluego.home.folder,"ClueGOExampleFiles/GSE6887_Bcell_Healthy_top200UpRegulated.txt",sep="/") # Linux / MacOSX
# file.location = paste(cluego.home.folder,"ClueGOExampleFiles\GSE6887_Bcell_Healthy_top200UpRegulated.txt",sep="\\") # Windows
gene.list <- toJSON(read.table(file.location,as.is=TRUE,header=FALSE)[[1]])
response <- PUT(url=paste(cluego.base.url,"cluster","upload-ids-list",URLencode(cluster1),sep="/"), body=gene.list, encode = "json", content_type_json())
#print(response)

# list all available ontologies
response <- GET(paste(cluego.base.url,"ontologies","get-ontology-info",sep="/"))
print(content(response))

# set ontologies
selected.ontologies <- toJSON(c("2;Ellipse","7;Triangle","8;Rectangle"))
response <- PUT(url=paste(cluego.base.url,"ontologies","set-ontologies",sep="/"), body=selected.ontologies, encode = "json", content_type_json())
#print(response)

# run the analysis
analysis.name <- "ClueGO Example analysis"
response <- GET(paste(cluego.base.url,URLencode(analysis.name),sep="/"))
#print(response)
print(content(response, encode = "text"))

# get network id
response <- GET(paste(cytoscape.base.url,"networks","currentNetwork",sep="/"))
current.network.suid <- content(response, encode = "json")$data$networkSUID
print(current.network.suid)

# get network graphics
image.type = "svg" # png, pdf
response <- GET(paste(cytoscape.base.url,"networks",current.network.suid,"views",paste("first.",image.type,sep=""),sep="/"))
image.file.name = paste(home.folder,paste("ClueGOExampleNetwork.",image.type,sep=""),sep="/") # Linux / MacOSX
# image.file.name = paste(home.folder,paste("ClueGOExampleNetwork.",image.type,sep=""),sep="\\") # Windows
writeBin(content(response, encode = "raw"),image.file.name)

# get analysis results
response <- GET(paste(cluego.base.url,"analysis-results","get-cluego-table",current.network.suid,sep="/"))
result.table.text <- content(response, encode = "text", encoding = "UTF-8")
result.table.data.frame <- text.to.data.frame(result.table.text)
table.file.name = paste(home.folder,"ClueGOExampleResultTable.txt",sep="/") # Linux / MacOSX
# table.file.name = paste(home.folder,"ClueGOExampleResultTable.txt",sep="\\") # Windows
write.table(result.table.data.frame,file=table.file.name,row.names=FALSE, na="",col.names=TRUE, sep="\t")
# print(result.table.data.frame)


