## Bioinformatics for biologists
#This Readme is just to describe the purpose; the site itself can be
#viewed [here](https://rwhetten.github.io/BIT815).

library(igraph)

m <- expand.grid(LETTERS[1:10], LETTERS[1:10])
m$Weight <- runif(nrow(m), 0.01, max = 5)
m <- m[m$Var1!=m$Var2, ] ##remove loop edges

g <- graph_from_data_frame(m)
plot(g, edge.curved=TRUE)

g1 <- induced_subgraph(g, c("A","B","C"))
plot(g1, edge.curved=TRUE)
#E(g1)$Weight
for(xcrit in sort(E(g1)$Weight)){
   if(any(degree(subgraph.edges(g1, E(g1)[Weight > xcrit], delete.vertices = F))==0)){
      print(paste0("Xcrit = ", xcrit))
      break
   }
}
#[1] "Xcrit = 4.5026016204874"

-#Often times a data sample is needed when seeking help, but it is also good practice to work with your own random data matrix.
#A little background of our control samples
CON = data.frame(Type=rep("RA", length=7), Time=c("0", "2", "4", "6", "8", "10", "12"))
CON = data.frame(Sample=paste(CON$Type, CON$Time, sep = "_"))
CON
#The Y matrix which describes the experimental design aka the metadata.
CON_Y <- as.matrix(data.frame(Duration = c(0, 2, 4, 6, 8, 10, 12)))

#The X matrix of their "raw" or non-normalized count values
set.seed(0408)
CON_X <- matrix(rnorm(17500, mean=5000, sd=500), ncol=7, nrow=2500)
min(CON_X)
max(CON_X)

#Add the samples to Y and X
#Let's tighten it up with the N and P descriptors for each matrix, respectively.
rownames(CON_Y) = CON$Sample
colnames(CON_X) = CON$Sample
rownames(CON_X) = paste0("Gene", 1:2500)
head(CON_X)
tail(CON_X)

#Logbase10 scale the raw count data
logCON <- log(CON_X, 10)
#DO PLS
library(mixOmics)
plsCON <- spls(t(logCON), CON_Y, ncomp = 6)
#notice how the logtransformed CON_X matrix was transformed 90 degrees, and I set the number of components to two.
attributes(plsCON)
plsCON$ncomp

set.seed(2017)
RA.valid <- perf(plsCON, validation = "loo")
#?perf to acknowledge and read about the available validation arguments
attributes(RA.valid)
RA.valid$Q2.total
RA.valid$R2
#Most of the variation was explained in the first components
Let's see what happens when I do this:
plsCON <- spls(t(logCON), CON_Y, ncomp = 2)
set.seed(2017)
RA.valid <- perf(plsCON, validation = "loo")
RA.valid$Q2.total
RA.valid$R2

#The individuals plot display
plotLoadings(plsCON, comp = 1:2)
#The variable plot display
plotVar(plsCON, comp = 1:2, plot=TRUE, Y.label = "D1", X.label = "D2", cex = c(0.5, 0.8))
#The biological network
network(plsCON, comp = 1:2, cutoff = 0.95, shape.node = c("rectangle", "rectangle"),
color.node = c("white", "pink"))

#Since we are exploring here, let's reduce the cutoff to a lower value, circle the meta data
network(plsCON, comp = 1:2, cutoff = 0.75, shape.node = c("rectangle", "circle"),
color.node = c("white", "pink"))

#Let's add the treatment samples
#A little background
TMT = data.frame(Type=rep("TA", length=7), Time=c("0", "2", "4", "6", "8", "10", "12"))
TMT = data.frame(Sample=paste(TMT$Type, TMT$Time, sep = "_"))

#The Y matrix which describes the experimental design aka the metadata.
TMT_Y <- as.matrix(data.frame(Duration = c(0, 2, 4, 6, 8, 10, 12)))

#The X matrix of their "raw" or non-normalized count values
set.seed(0408)
TMT_X <- matrix(rnorm(17500, mean=10000, sd=1000), ncol=7, nrow=2500)
min(TMT_X)
max(TMT_X)

#Add the samples to Y and X
#Let's tighten it up with the N and P descriptors for each matrix, respectively.
rownames(TMT_Y) = TMT$Sample
colnames(TMT_X) = TMT$Sample
rownames(TMT_X) = paste0("Gene", 1:2500)
head(TMT_X)
tail(TMT_X)

#Logbase10 scale the raw count data
logTMT <- log(TMT_X, 10)
#DO PLS
plsTMT <- spls(t(logTMT), TMT_Y, ncomp = 2)
set.seed(2017)
TA.valid <- perf(plsTMT, validation = "loo")
TA.valid$Q2.total
TA.valid$R2

#The individuals plot display
plotIndiv(plsTMT, comp = 1:2)
#The variable plot display
plotVar(plsTMT, comp = 1:2, plot=TRUE, Y.label = "D1", X.label = "D2", cex = c(0.5, 0.8))
#The biological network
network(plsTMT, comp = 1:2, cutoff = 0.95, shape.node = c("rectangle", "rectangle"),
color.node = c("white", "pink"))

#Since we are exploring here, let's reduce the cutoff to a lower value, circle the meta data
network(plsTMT, comp = 1:2, cutoff = 0.75, shape.node = c("rectangle", "circle"),
color.node = c("white", "pink"))


Well, since we have both CON and TMT log10 normalized data let's build a bigger network.
X <- cbind(logCON, logTMT)
class(X)
Y <- as.matrix(data.frame(Duration = rep(c(0, 2, 4, 6, 8, 10, 12), times=2),
Treatment=rep(c(0,1), each=7),
#and let's add in a third variable to account for the different freezer space used to store each sample
Storage=rep(c(1,2), each=7)))

#uh oh, I need to name the rows...!!??
rnms <- rbind(CON, TMT)
rownames(Y) = rnms$Sample

plsX <- spls(t(X), Y, ncomp = 6)

#The individuals plot display
plotIndiv(plsX, comp = 1:2)
#The variable plot display
plotVar(plsX, comp = 1:2, plot=TRUE, Y.label = "D1", X.label = "D2", cex = c(0.5, 0.8))
plotVar(plsX, comp = 2:3, plot=TRUE, Y.label = "D2", X.label = "D3", cex = c(0.5, 0.8))
#Oh, Ok, so most of the variation is explained in the first and second dimensions
#The biological network
network(plsX, comp = 1:2, cutoff = 0.65, shape.node = c("rectangle", "rectangle"),
color.node = c("white", "pink"))

#Since we are exploring here, let's increase the cutoff to a higher value
network(plsX, comp = 1:2, cutoff = 0.989999, shape.node = c("rectangle", "rectangle"),
color.node = c("white", "pink"))
