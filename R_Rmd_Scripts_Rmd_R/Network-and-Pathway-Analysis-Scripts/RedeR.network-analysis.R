### R code from vignette source 'RedeR.Rnw'
library (RedeR)
rdp <- RedPort()
calld(rdp)

??RedeR
###################################################
### code chunk number 3: Add graph
###################################################
g1 <- graph.lattice(c(5,5,5))
addGraph( rdp, g1, layout.kamada.kawai(g1) )


###################################################
### code chunk number 4: Get graph
###################################################
g2 <- getGraph(rdp)
resetd(rdp)

g3 <- barabasi.game(10)
g4 <- barabasi.game(10)
V(g3)$name<-paste("sn",1:10,sep="")
V(g4)$name<-paste("sm",1:10,sep="")
addGraph(rdp, g3, isNest =TRUE, gcoord=c(25,25), gscale=50)
addGraph(rdp, g4, isNest =TRUE, gcoord=c(75,75), gscale=50)


###################################################
### code chunk number 6: Get subgraph
###################################################
selectNodes(rdp,"N0")
g5 <- getGraph(rdp, status= "selected")
resetd(rdp)


g6 <- barabasi.game(500)
addGraph(rdp, g6, zoom=20)


###################################################
### code chunk number 8: Start relax
###################################################
relax(rdp,p2=400,p5=30,ps=T)


###################################################
### code chunk number 9: Map clic communities
###################################################
g <- getGraph(rdp, status= "selected")
if(vcount(g)>0)plot(degree.distribution(g), xlab = "k", ylab = "P(k)", pch=19)
resetd(rdp)


###################################################
### code chunk number 10: Workflow 1: get a dataframe and an interactome
###################################################
data(ER.limma)
data(hs.inter)
dt <- ER.limma
gi <- hs.inter


###################################################
### code chunk number 11: Workflow 1: extract a subgraph and set attributes to RedeR
###################################################
gt3  <- subg(g=gi, dat=dt[dt$degenes.t3!=0,], refcol=1)
gt3  <- att.setv(g=gt3, from="Symbol", to="nodeAlias")
gt3  <- att.setv(g=gt3, from="logFC.t3", to="nodeColor", breaks=seq(-2,2,0.4), pal=2)


###################################################
### code chunk number 12: Workflow 1: extract another subgraph and set attributes to RedeR
###################################################
gt6  <- subg(g=gi, dat=dt[dt$degenes.t6!=0,], refcol=1)
gt6  <- att.setv(g=gt6, from="Symbol", to="nodeAlias")
gt6  <- att.setv(g=gt6, from="logFC.t6", to="nodeColor", breaks=seq(-2,2,0.4), pal=2)


###################################################
### code chunk number 13: Workflow 1: extract another subgraph and set attributes to RedeR
###################################################
gt12 <- subg(g=gi, dat=dt[dt$degenes.t12!=0,], refcol=1)
gt12 <- att.setv(g=gt12, from="Symbol", to="nodeAlias")
gt12 <- att.setv(g=gt12, from="logFC.t12", to="nodeColor", breaks=seq(-2,2,0.4), pal=2)


###################################################
### code chunk number 14: Workflow 1: add subgraphs to the app
###################################################
addGraph(rdp, gt3, gcoord=c(10,25), gscale=20, isNest=TRUE, theme='tm1', zoom=30)
addGraph(rdp, gt6, gcoord=c(20,70), gscale=50, isNest=TRUE, theme='tm1', zoom=30)
addGraph(rdp, gt12, gcoord=c(70,55), gscale=80, isNest=TRUE, theme='tm1', zoom=30)


###################################################
### code chunk number 15: Workflow 1: nest subgraphs
###################################################
nestNodes(rdp, nodes=V(gt3)$name, parent="N1", theme='tm2')
nestNodes(rdp, nodes=V(gt6)$name, parent="N2", theme='tm2')
nestNodes(rdp, nodes=V(gt3)$name, parent="N4", theme='tm3')


###################################################
### code chunk number 16: Workflow 1: assign edges to containers
###################################################
mergeOutEdges(rdp)


###################################################
### code chunk number 17: Workflow 1: relax the network
###################################################


















































































































































































































































relax(rdp,50,400)


###################################################
### code chunk number 18: Workflow 1: add a color legend (other types are available)
###################################################
scl <- gt3$legNodeColor$scale
leg <- gt3$legNodeColor$legend
addLegend.color(rdp, colvec=scl, labvec=leg, title="node color (logFC)")


###################################################
### code chunk number 19: Workflow 1: select a gene
###################################################
selectNodes(rdp,"RET")


###################################################
### code chunk number 20: Workflow 1: reset graph
###################################################
resetd(rdp)


data(ER.deg)
dt <- ER.deg$dat
sg <- ER.deg$ceg


sg <- att.mapv(sg, dat=dt, refcol=1)


###################################################
### code chunk number 23: Workflow 2: set attributes to RedeR
###################################################
sg <- att.setv(sg, from="Symbol", to="nodeAlias")
sg <- att.setv(sg, from="logFC.t3", to="nodeColor", breaks=seq(-1,1,0.2), pal=2)
sg <- att.setv(sg, from="ERbdist", to="nodeSize", nquant=10, isrev=TRUE, xlim=c(5,40,1))


###################################################
### code chunk number 24: Workflow 2: add graph to the app
###################################################
addGraph(rdp,sg)


###################################################
### code chunk number 25: Workflow 2: compute a hierarchical clustering using standard R functions
###################################################
hc <- hclust(dist(get.adjacency(sg, attr="weight")))


nesthc(rdp,hc, cutlevel=3, nmemb=5, cex=0.3, labels=V(sg)$nodeAlias)


mergeOutEdges(rdp,nlev=2)


###################################################
### code chunk number 28: Workflow 2: relax the network
###################################################
relax(rdp)

scl <- sg$legNodeColor$scale
leg <- sg$legNodeColor$legend
addLegend.color(rdp, colvec=scl, labvec=leg, title="diff. gene expression (logFC)")

scl <- sg$legNodeSize$scale
leg <- sg$legNodeSize$legend
addLegend.size(rdp, sizevec=scl, labvec=leg, title="bd site distance (kb)")

resetd(rdp)

print(sessionInfo(), locale=FALSE)


