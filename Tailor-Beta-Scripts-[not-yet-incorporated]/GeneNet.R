library(GeneNet)
??GeneNet
??corpcor
install.packages( c("corpcor",
                    "longitudinal", "fdrtool",
                    "locfdr", "GeneNet") )
library(corpcor);library(longitudinal);library(fdrtool)
library(locfdr)
data(ecoli)
dim(ecoli)
str(ecoli)
inferred.pcor<-ggm.estimate.pcor(ecoli)
dim(inferred.pcor)

test.results<-ggm.test.edges(inferred.pcor)
dim(test.results)
signif<-test.results$prob > 0.80
sum(signif)

test.results[signif,]

node.labels <- colnames(ecoli)

gr <- ggm.make.graph(
   test.results[signif,],
   node.labels)
show.edge.weights(gr)

ggm.plot.graph(gr,
               show.edge.labels=FALSE,
               layoutType="fdp")

layoutType="fdp"
layoutType="neato"
show.edge.labels=TRUE



