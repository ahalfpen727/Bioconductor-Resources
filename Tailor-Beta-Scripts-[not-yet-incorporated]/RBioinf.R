### R code from vignette source 'RBioinf.Rnw'

###################################################
### code chunk number 1: ex1
###################################################
library(RBioinf)
setClass("object")
setClass("grid-layout", contains="object")
setClass("horizontal-grid", contains="grid-layout")
setClass("vertical-grid", contains="grid-layout")
setClass("hv-grid", contains=c("horizontal-grid", "vertical-grid"))

LPO("hv-grid")



###################################################
### code chunk number 2: cpotest
###################################################
computeClassLinearization("object")
computeClassLinearization("grid-layout")
computeClassLinearization("vertical-grid")


###################################################
### code chunk number 3: s1
###################################################
setClass("vh-grid", contains=c("vertical-grid", "horizontal-grid"))

setClass("confused", contains=c("hv-grid", "vh-grid"))

LPO("vh-grid")
tryCatch(LPO("confused"), error=function(x) "this one failed")



###################################################
### code chunk number 4: plotconfG
###################################################
 library(Rgraphviz)
 confG = class2Graph("confused")

 cGa = makeNodeAttrs(confG, shape="ellipse", fill="grey", width=4)

 plot(confG, nodeAttrs=cGa)



###################################################
### code chunk number 5: showExtends
###################################################

 setClass("a")
 setClass("b")
 setClass("c", contains = c("a", "b"))
 setClass("d", contains = c("b", "a"))

 extends("c")
 extends("d")

 setClass("e", contains=c("c", "d"))



###################################################
### code chunk number 6: sCdemo
###################################################

getAllSuperClasses(getClass("e"))

cD = superClassDepth(getClass("e"))
cD$label
cD$depth

superClasses(getClass("e"))



###################################################
### code chunk number 7: c2G
###################################################

cH = class2Graph("e")



###################################################
### code chunk number 8: editWin
###################################################

 setClass("pane", contains="object")
 setClass("editing-mixin", contains="object")
 setClass("scrolling-mixin", contains="object")

 setClass("scrollable-pane", contains=c("pane", "scrolling-mixin"))
 setClass("editable-pane", contains=c("pane", "editing-mixin"))

 setClass("editable-scrollable-pane",
         contains=c("scrollable-pane", "editable-pane"))



###################################################
### code chunk number 9: LPOseW
###################################################

LPO("editable-scrollable-pane")
LPO("editable-scrollable-pane", C3=TRUE)



###################################################
### code chunk number 10: eWgraph
###################################################
 eWG = class2Graph("editable-scrollable-pane")

 eWGattrs = makeNodeAttrs(eWG, shape="ellipse", fill="grey", width=4)

 plot(eWG, nodeAttrs=eWGattrs)



###################################################
### code chunk number 11: RBioinf.Rnw:324-325
###################################################
sessionInfo()


