library(RColorBrewer)
library("gplots")


#############
# Explanation for the new keep_order parameter:
# keep_order can accept the numbers 1 or 2, in which case it will make sure the resulting merged data.frame will be ordered according to the original order of rows of the data.frame entered to x (if keep_order=1) or to y (if keep_order=2).  If keep_order is missing, merge will continue working as usual.  If keep_order gets some input other then 1 or 2, it will issue a warning that it doesn't accept these values, but will continue working as merge normally would.  Notice that the parameter "sort" is practically overridden when using keep_order (with the value 1 or 2).

# an example is offered at the end of this code chunk

merge.data.frame <- function (x, y, by = intersect(names(x), names(y)), by.x = by, 
                              by.y = by, all = FALSE, all.x = all, all.y = all, sort = TRUE, 
                              suffixes = c(".x", ".y"), incomparables = NULL, keep_order, ...) 
{
                                        # if we use the "keep_order" parameter, we might need to modify either the x or y data.frame objects: (either by placing 1 or "x", or by placing 2 or "y")
    if(!missing(keep_order))
	{
            
                                        # some functions we will use soon:
            add.id.column.to.data <- function(DATA)
		{
                    data.frame(DATA, id... = seq_len(nrow(DATA)))
		}
                                        # example:
                                        # add.id.column.to.data(data.frame(x = rnorm(5), x2 = rnorm(5)))
            order.by.id...and.remove.it <- function(DATA)
		{
                                        # gets in a data.frame with the "id..." column.  Orders by it and returns it
                    if(!any(colnames(DATA)=="id...")) 
                        {
                            warning("The function order.by.id...and.remove.it is useful only for data.frame objects that includes the 'id...' order column")
                            return(DATA)
                        }
                    
                    
                    ss_r <- order(DATA$id...)
                    ss_c <- colnames(DATA) != "id..."
                    DATA[ss_r, ss_c]		
		}
                                        # example: 
                                        # set.seed(3424)
                                        # x  <- data.frame(x = rnorm(5), x2 = rnorm(5))
                                        # x2 <- add.id.column.to.data(x)[c(1,4,2,5,3),]
                                        # x2
                                        # order.by.id...and.remove.it(x2)
            if(keep_order == "x") keep_order <- 1
            if(keep_order == "y") keep_order <- 2

            if(keep_order == 1) x<-add.id.column.to.data(x)
            if(keep_order == 2) y<-add.id.column.to.data(y)
                                        # if you didn't get 1 or 2 - issue a warning:
            if(!(any(keep_order == c(1,2)) ))  warning("The parameter 'keep_order' in the function merge.data.frame only accepts the values 1 (for the x data.frame) or 2 (for the y data.frame)")

                                        # sort <- FALSE
                                        # notice that if sort was TRUE, using the keep_order parameter will eventually override it...
	}

    fix.by <- function(by, df) {
        if (is.null(by)) 
            by <- numeric()
        by <- as.vector(by)
        nc <- ncol(df)
        if (is.character(by)) 
            by <- match(by, c("row.names", names(df))) - 1L
        else if (is.numeric(by)) {
            if (any(by < 0L) || any(by > nc)) 
                stop("'by' must match numbers of columns")
        }
        else if (is.logical(by)) {
            if (length(by) != nc) 
                stop("'by' must match number of columns")
            by <- seq_along(by)[by]
        }
        else stop("'by' must specify column(s) as numbers, names or logical")
        if (any(is.na(by))) 
            stop("'by' must specify valid column(s)")
        unique(by)
    }
    nx <- nrow(x <- as.data.frame(x))
    ny <- nrow(y <- as.data.frame(y))
    by.x <- fix.by(by.x, x)
    by.y <- fix.by(by.y, y)
    if ((l.b <- length(by.x)) != length(by.y)) 
        stop("'by.x' and 'by.y' specify different numbers of columns")
    if (l.b == 0L) {
        nm <- nm.x <- names(x)
        nm.y <- names(y)
        has.common.nms <- any(cnm <- nm.x %in% nm.y)
        if (has.common.nms) {
            names(x)[cnm] <- paste(nm.x[cnm], suffixes[1L], sep = "")
            cnm <- nm.y %in% nm
            names(y)[cnm] <- paste(nm.y[cnm], suffixes[2L], sep = "")
        }
        if (nx == 0L || ny == 0L) {
            res <- cbind(x[FALSE, ], y[FALSE, ])
        }
        else {
            ij <- expand.grid(seq_len(nx), seq_len(ny))
            res <- cbind(x[ij[, 1L], , drop = FALSE], y[ij[, 
                                           2L], , drop = FALSE])
        }
    }
    else {
        if (any(by.x == 0L)) {
            x <- cbind(Row.names = I(row.names(x)), x)
            by.x <- by.x + 1L
        }
        if (any(by.y == 0L)) {
            y <- cbind(Row.names = I(row.names(y)), y)
            by.y <- by.y + 1L
        }
        row.names(x) <- NULL
        row.names(y) <- NULL
        if (l.b == 1L) {
            bx <- x[, by.x]
            if (is.factor(bx)) 
                bx <- as.character(bx)
            by <- y[, by.y]
            if (is.factor(by)) 
                by <- as.character(by)
        }
        else {
            bx <- x[, by.x, drop = FALSE]
            by <- y[, by.y, drop = FALSE]
            names(bx) <- names(by) <- paste("V", seq_len(ncol(bx)), 
                                            sep = "")
            bz <- do.call("paste", c(rbind(bx, by), sep = "\r"))
            bx <- bz[seq_len(nx)]
            by <- bz[nx + seq_len(ny)]
        }
        comm <- match(bx, by, 0L)
        bxy <- bx[comm > 0L]
        xinds <- match(bx, bxy, 0L, incomparables)
        yinds <- match(by, bxy, 0L, incomparables)
        if (nx > 0L && ny > 0L) 
            m <- .Internal(merge(xinds, yinds, all.x, all.y))
        else m <- list(xi = integer(), yi = integer(), x.alone = seq_len(nx), 
                       y.alone = seq_len(ny))
        nm <- nm.x <- names(x)[-by.x]
        nm.by <- names(x)[by.x]
        nm.y <- names(y)[-by.y]
        ncx <- ncol(x)
        if (all.x) 
            all.x <- (nxx <- length(m$x.alone)) > 0L
        if (all.y) 
            all.y <- (nyy <- length(m$y.alone)) > 0L
        lxy <- length(m$xi)
        has.common.nms <- any(cnm <- nm.x %in% nm.y)
        if (has.common.nms) 
            nm.x[cnm] <- paste(nm.x[cnm], suffixes[1L], sep = "")
        x <- x[c(m$xi, if (all.x) m$x.alone), c(by.x, seq_len(ncx)[-by.x]), 
               drop = FALSE]
        names(x) <- c(nm.by, nm.x)
        if (all.y) {
            ya <- y[m$y.alone, by.y, drop = FALSE]
            names(ya) <- nm.by
            ya <- cbind(ya, x[rep.int(NA_integer_, nyy), nm.x, 
                              drop = FALSE])
            x <- rbind(x, ya)
        }
        if (has.common.nms) {
            cnm <- nm.y %in% nm
            nm.y[cnm] <- paste(nm.y[cnm], suffixes[2L], sep = "")
        }
        y <- y[c(m$yi, if (all.x) rep.int(1L, nxx), if (all.y) m$y.alone), 
               -by.y, drop = FALSE]
        if (all.x) {
            for (i in seq_along(y)) is.na(y[[i]]) <- (lxy + 1L):(lxy + 
                                                                     nxx)
        }
        if (has.common.nms) 
            names(y) <- nm.y
        res <- cbind(x, y)
        if (sort) 
            res <- res[if (all.x || all.y) 
                do.call("order", x[, seq_len(l.b), drop = FALSE])
            else sort.list(bx[m$xi]), , drop = FALSE]
    }
    attr(res, "row.names") <- .set_row_names(nrow(res))

	if(!missing(keep_order) && any(keep_order == c(1,2))) return(order.by.id...and.remove.it(res))	
	# notice how it is essential to use && here, since if the first argument is false, it will not be possible to evaluate the second argument
	
    res
}

library(readxl)
A.over.wt1021 <- read_excel("A-over-1021.xlsx")
AB.over.wt1021 <- read_excel("AB-over-1021.xlsx")
head(A.over.wt1021)
head(AB.over.wt1021)

colnames(A.over.wt1021)
head(AB.over.wt1021)


## Delta = "Delta"

## A       = eval(paste0(Delta, "cbrA pTH1227"))
## A       = expression(paste0(Delta, "cbrA pTH1227"))
A       = expression(Delta * "cbrA pTH1227")
## A       = paste0(Delta, "cbrA pTH1227")
## A       = "Delta~cbrA pTH1227"

## AB      = eval(paste0(Delta, "cbrA pTH1227::morA"))
## AB      = expression(paste0(Delta, "cbrA pTH1227::morA"))
AB      = expression(Delta * "cbrA pTH1227::morA")
## AB      = "Delta~cbrA pTH1227::morA"

wt1021  = "WT pTH1227"
wt1021B = "WT pTH1227::morA"


## AB.over.wt1021 = read.table(file="AB-over-wt1021.names.csv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
## A.over.wt1021 = read.table(file="A-over-wt1021.names.csv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
AB.over.wt1021 = read.table(file="AB-over-1021_gene_exp.diff", sep="\t", header=TRUE, stringsAsFactors=FALSE)
A.over.wt1021 = read.table(file="A-over-1021_gene_exp.diff", sep="\t", header=TRUE, stringsAsFactors=FALSE)

AB.over.wt1021 = read.table(file="AB-over-1021", sep="\t", header=TRUE, stringsAsFactors=FALSE)
A.over.wt1021 = read.table(file="A-over-1021_gene_exp.diff", sep="\t", header=TRUE, stringsAsFactors=FALSE)


AB.over.wt1021$gene_id = gsub("gene:", "", AB.over.wt1021$gene_id)
A.over.wt1021$gene_id = gsub("gene:", "", A.over.wt1021$gene_id)

AB.over.wt1021$FC = 2^AB.over.wt1021$log2.fold_change.
A.over.wt1021$FC = 2^A.over.wt1021$log2.fold_change.


## AB.over.wt1021$gene_lc = tolower(AB.over.wt1021$gene)
## A.over.wt1021$gene_lc = tolower(A.over.wt1021$gene)

head(AB.over.wt1021)
head(A.over.wt1021)


#####################################################################################################################
#####################################################################################################################
##  Gene Table 1
#####################################################################################################################
#####################################################################################################################
data1 = read.table(file="heatmap.data.new.1.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
## colnames(data1) = c("GeneName","locus","GOterm", "rescued")
data1=data2
colnames(data1) = c("GeneName","gene_id","GOterm", "rescued")
head(data1)

## data1$gene_lc = tolower(data1$gene)
## data1$locus_lc = tolower(data1$locus)

## AB.over.wt1021.data1 = AB.over.wt1021[(AB.over.wt1021$GeneName %in% data1$gene) | (AB.over.wt1021$GeneName %in% data1$locus),]
## A.over.wt1021.data1 = A.over.wt1021[(A.over.wt1021$GeneName %in% data1$gene) | (A.over.wt1021$GeneName %in% data1$locus),]

## AB.over.wt1021.data1 = AB.over.wt1021[(AB.over.wt1021$gene %in% data1$gene) | (AB.over.wt1021$gene %in% data1$locus),]
## A.over.wt1021.data1 = A.over.wt1021[(A.over.wt1021$gene %in% data1$gene) | (A.over.wt1021$gene %in% data1$locus),]

## AB.over.wt1021.data1 = AB.over.wt1021[(AB.over.wt1021$gene_id %in% data1$gene) | (AB.over.wt1021$gene_id %in% data1$locus),]
## A.over.wt1021.data1 = A.over.wt1021[(A.over.wt1021$gene_id %in% data1$gene) | (A.over.wt1021$gene_id %in% data1$locus),]

## AB.over.wt1021.data1 = AB.over.wt1021[(AB.over.wt1021$gene_lc %in% data1$gene_lc) | (AB.over.wt1021$gene_lc %in% data1$locus_lc),]
## A.over.wt1021.data1 = A.over.wt1021[(A.over.wt1021$gene_lc %in% data1$gene_lc) | (A.over.wt1021$gene_lc %in% data1$locus_lc),]


##############################################################################################
# inner join with keeping order in data1
##############################################################################################
AB.over.wt1021.data1 = merge(AB.over.wt1021, data1, by=c("gene_id"), keep_order=2) # inner join with keeping order in data1
A.over.wt1021.data1 = merge(A.over.wt1021, data1, by=c("gene_id"), keep_order=2) # inner join with keeping order in data1


head(AB.over.wt1021.data1)
head(A.over.wt1021.data1)

dim(data1)
dim(AB.over.wt1021)
dim(AB.over.wt1021.data1)

dim(A.over.wt1021)
dim(A.over.wt1021.data1)

## heatmap.data = merge(A.over.wt1021.data1, AB.over.wt1021.data1, by=c("GeneID","GeneName","Length"), suffixes=c(".A",".AB"))
heatmap.data = merge(A.over.wt1021.data1, AB.over.wt1021.data1, by=c("gene_id","gene","GeneName","test_id","locus","GOterm","rescued"), suffixes=c(".A",".AB"))
dim(heatmap.data)
head(heatmap.data)

notSig = (heatmap.data$q_value.A > 0.05) | ((heatmap.data$FC.A > 0.75) & (heatmap.data$FC.A < 1.5))
heatmap.data$newFC.A  = heatmap.data$FC.A
heatmap.data[notSig,]$newFC.A  = NA

notSig = (heatmap.data$q_value.AB > 0.05) | ((heatmap.data$FC.AB > 0.75) & (heatmap.data$FC.AB < 1.5))
heatmap.data$newFC.AB  = heatmap.data$FC.AB
heatmap.data[notSig,]$newFC.AB  = NA

## heatmap.data2 = t(cbind(heatmap.data$FC.A, heatmap.data$FC.AB))
heatmap.data2 = t(cbind(heatmap.data$newFC.A, heatmap.data$newFC.AB))
colnames(heatmap.data2) = heatmap.data$GeneName
## colnames(heatmap.data2) = heatmap.data$gene
rownames(heatmap.data2) = c(A, AB)
## rownames(heatmap.data2) = c(bquote(.(parse(text=A))), bquote(.(parse(text=AB))))
## rownames(heatmap.data2) = c(eval(parse(text=A)), eval(parse(text=AB)))
# rownames(heatmap.data2) = c(expression(paste(Delta, "cbrA pTH1227")), expression(paste(Delta, "cbrA pTH1227")))

dim(heatmap.data2)
head(heatmap.data2)

quartz() # new window
heatmap.2(heatmap.data2, # values
          Rowv = NA, # no row clustering
          Colv = NA, # no column clustering
          scale = "row", # row-wise scaling
          na.color="black",
          ## col=greenred(75), # green to red heatmap
          ## col=blueyellow(75), # green to red heatmap
          col=colorRampPalette(c('blue', 'yellow')),
          labRow=c(A, AB), # MUST USE LABROW TO PROPERLY SHOW AN EXPRESSION!!!
          dendrogram="none", # no dendrogram
          key=TRUE, # add color key
          symkey=FALSE,
          density.info="none",
          trace="none",
          cexRow=0.6 # make row labels small
          )

##############################################################################################
# sort by GOterm and then by rescued
##############################################################################################
## heatmap.data = heatmap.data[order(heatmap.data$GOterm, -heatmap.data$rescued),]
heatmap.data = heatmap.data[order(heatmap.data$GOterm, -as.numeric(as.factor(heatmap.data$rescued))),]

## heatmap.data2 = t(cbind(heatmap.data$FC.A, heatmap.data$FC.AB))
heatmap.data2 = t(cbind(heatmap.data$newFC.A, heatmap.data$newFC.AB))
colnames(heatmap.data2) = heatmap.data$GeneName
## colnames(heatmap.data2) = heatmap.data$gene
rownames(heatmap.data2) = c(A, AB)
## rownames(heatmap.data2) = c(bquote(.(parse(text=A))), bquote(.(parse(text=AB))))
## rownames(heatmap.data2) = c(eval(parse(text=A)), eval(parse(text=AB)))
# rownames(heatmap.data2) = c(expression(paste(Delta, "cbrA pTH1227")), expression(paste(Delta, "cbrA pTH1227")))

dim(heatmap.data2)
head(heatmap.data2)

quartz() # new window
heatmap.2(heatmap.data2, # values
          
          Rowv = NA, # no row clustering
          Colv = NA, # no column clustering
          scale = "row", # row-wise scaling
          na.color="black",
          ## col=greenred(75), # green to red heatmap
          ## col=blueyellow(75), # green to red heatmap
          col=colorRampPalette(c('blue', 'yellow')),
          labRow=c(A, AB), # MUST USE LABROW TO PROPERLY SHOW AN EXPRESSION!!!
          dendrogram="none", # no dendrogram
          key=TRUE, # add color key
          symkey=FALSE,
          density.info="none",
          trace="none",
          cexRow=0.6 # make row labels small
          )


#####################################################################################################################
#####################################################################################################################
##  Gene Table 2
#####################################################################################################################
#####################################################################################################################
data2 = read.table(file="heatmap.data.new.2.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
## colnames(data2) = c("gene","locus","GOterm")
colnames(data2) = c("GeneName","gene_id","Rescued", "Both")
head(data2)

data2 <- read_excel("heatmap data final manuscript.xlsx", 
                                            sheet = "CtrA regulated")
colnames(data2) = c("GeneName","gene_id","Rescued", "Both")
head(data2)

## AB.over.wt1021.data2 = AB.over.wt1021[(AB.over.wt1021$GeneName %in% data2$gene),]
## AB.over.wt1021.data2 = AB.over.wt1021[(AB.over.wt1021$GeneName %in% data2$locus),]
## AB.over.wt1021.data2 = AB.over.wt1021[(AB.over.wt1021$GeneName %in% data2$gene) | (AB.over.wt1021$GeneName %in% data2$locus),]
## A.over.wt1021.data2 = A.over.wt1021[(A.over.wt1021$GeneName %in% data2$gene) | (A.over.wt1021$GeneName %in% data2$locus),]

##############################################################################################
# inner join with keeping order in data2
##############################################################################################
AB.over.wt1021.data2 = merge(AB.over.wt1021, data2, by=c("gene_id"), keep_order=2) # inner join with keeping order in data2
A.over.wt1021.data2 = merge(A.over.wt1021, data2, by=c("gene_id"), keep_order=2) # inner join with keeping order in data2

head(AB.over.wt1021.data2)
head(A.over.wt1021.data2)

dim(data2)
dim(AB.over.wt1021)
dim(AB.over.wt1021.data2)

dim(A.over.wt1021)
dim(A.over.wt1021.data2)

## heatmap.data = merge(A.over.wt1021.data2, AB.over.wt1021.data2, by=c("GeneID","GeneName","Length"), suffixes=c(".A",".AB"))
heatmap.data = merge(A.over.wt1021.data2, AB.over.wt1021.data2, 
                     by=c("gene_id","gene","GeneName","test_id","locus","Both","Rescued"), suffixes=c(".A",".AB"))
heatmap.data = merge(A.over.wt1021.data2, AB.over.wt1021.data2, 
                     by=c("gene_id","gene","GeneName","test_id","locus","Both","Rescued"), suffixes=c(".A",".AB"), keep_order = 1)

dim(heatmap.data)
heatmap.data[1:50,]

notSig = (heatmap.data$q_value.A > 0.05) | ((heatmap.data$FC.A > 0.75) & (heatmap.data$FC.A < 1.5))
heatmap.data$newFC.A  = heatmap.data$FC.A
heatmap.data[notSig,]$newFC.A  = NA

notSig = (heatmap.data$q_value.AB > 0.05) | ((heatmap.data$FC.AB > 0.75) & (heatmap.data$FC.AB < 1.5))
heatmap.data$newFC.AB  = heatmap.data$FC.AB
heatmap.data[notSig,]$newFC.AB  = NA

## heatmap.data2 = t(cbind(heatmap.data$FC.A, heatmap.data$FC.AB))
heatmap.data2 = t(cbind(heatmap.data$newFC.A, heatmap.data$newFC.AB))
colnames(heatmap.data2) = heatmap.data$GeneName
rownames(heatmap.data2) = c(A, AB)
## rownames(heatmap.data2) = c(bquote(.(parse(text=A))), bquote(.(parse(text=AB))))
## rownames(heatmap.data2) = c(eval(parse(text=A)), eval(parse(text=AB)))
# rownames(heatmap.data2) = c(expression(paste(Delta, "cbrA pTH1227")), expression(paste(Delta, "cbrA pTH1227")))

dim(heatmap.data2)
head(heatmap.data2)

quartz() # new window
pdf(file = "heatmap2-part2.pdf")
heatmap.2(heatmap.data2, # values
          Rowv = NA, # no row clustering
          Colv = NA, # no column clustering
          scale = "row", # row-wise scaling
          na.color="black",
          ## col=greenred(75), # green to red heatmap
          ## col=blueyellow(75), # green to red heatmap
          col=colorRampPalette(c('blue', 'yellow')),
          labRow=c(A, AB), # MUST USE LABROW TO PROPERLY SHOW AN EXPRESSION!!!
          dendrogram="none", # no dendrogram
          key=TRUE, # add color key
          symkey=FALSE,
          density.info="none",
          trace="none",
          cexRow=0.6 # make row labels small
          )
dev.off()

##############################################################################################
# sort by GOterm and then by rescued
##############################################################################################
## heatmap.data = heatmap.data[order(heatmap.data$Rescued, -heatmap.data$Both),]
heatmap.data = heatmap.data[order(heatmap.data$Rescued, -as.numeric(as.factor(heatmap.data$Both))),]

## heatmap.data2 = t(cbind(heatmap.data$FC.A, heatmap.data$FC.AB))
heatmap.data2 = t(cbind(heatmap.data$newFC.A, heatmap.data$newFC.AB))
colnames(heatmap.data2) = heatmap.data$GeneName
rownames(heatmap.data2) = c(A, AB)
## rownames(heatmap.data2) = c(bquote(.(parse(text=A))), bquote(.(parse(text=AB))))
## rownames(heatmap.data2) = c(eval(parse(text=A)), eval(parse(text=AB)))
# rownames(heatmap.data2) = c(expression(paste(Delta, "cbrA pTH1227")), expression(paste(Delta, "cbrA pTH1227")))

dim(heatmap.data2)
head(heatmap.data2)

quartz() # new window
heatmap.2(heatmap.data2, # values
          Rowv = NA, # no row clustering
          Colv = NA, # no column clustering
          scale = "row", # row-wise scaling
          na.color="black",
          ## col=greenred(75), # green to red heatmap
          ## col=blueyellow(75), # green to red heatmap
          col=colorRampPalette(c('blue', 'yellow')),
          labRow=c(A, AB), # MUST USE LABROW TO PROPERLY SHOW AN EXPRESSION!!!
          dendrogram="none", # no dendrogram
          key=TRUE, # add color key
          symkey=FALSE,
          density.info="none",
          trace="none",
          cexRow=0.6 # make row labels small
          )


#####################################################################################################################
#####################################################################################################################
##  Gene Table 3
#####################################################################################################################
#####################################################################################################################
data3 = read.table(file="heatmap.data.new.3.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
## colnames(data3) = c("gene","locus","ClusterNumber")
colnames(data3) = c("GeneName","gene_id","ClusterNumber", "Rescued")
head(data3)

## AB.over.wt1021.data3 = AB.over.wt1021[(AB.over.wt1021$GeneName %in% data3$gene),]
## AB.over.wt1021.data3 = AB.over.wt1021[(AB.over.wt1021$GeneName %in% data3$locus),]
## AB.over.wt1021.data3 = AB.over.wt1021[(AB.over.wt1021$GeneName %in% data3$gene) | (AB.over.wt1021$GeneName %in% data3$locus),]
## A.over.wt1021.data3 = A.over.wt1021[(A.over.wt1021$GeneName %in% data3$gene) | (A.over.wt1021$GeneName %in% data3$locus),]

##############################################################################################
# inner join with keeping order in data3
##############################################################################################
AB.over.wt1021.data3 = merge(AB.over.wt1021, data3, by=c("gene_id"), keep_order=2) # inner join with keeping order in data3
A.over.wt1021.data3 = merge(A.over.wt1021, data3, by=c("gene_id"), keep_order=2) # inner join with keeping order in data3

head(AB.over.wt1021.data3)
head(A.over.wt1021.data3)

dim(data3)
dim(AB.over.wt1021)
dim(AB.over.wt1021.data3)

dim(A.over.wt1021)
dim(A.over.wt1021.data3)

## heatmap.data = merge(A.over.wt1021.data3, AB.over.wt1021.data3, by=c("GeneID","GeneName","Length"), suffixes=c(".A",".AB"))
heatmap.data = merge(A.over.wt1021.data3, AB.over.wt1021.data3, by=c("gene_id","gene","GeneName","test_id","locus","ClusterNumber","Rescued"), suffixes=c(".A",".AB"))
head(heatmap.data)
dim(heatmap.data)

notSig = (heatmap.data$q_value.A > 0.05) | ((heatmap.data$FC.A > 0.75) & (heatmap.data$FC.A < 1.5))
heatmap.data$newFC.A  = heatmap.data$FC.A
heatmap.data[notSig,]$newFC.A  = NA

notSig = (heatmap.data$q_value.AB > 0.05) | ((heatmap.data$FC.AB > 0.75) & (heatmap.data$FC.AB < 1.5))
heatmap.data$newFC.AB  = heatmap.data$FC.AB
heatmap.data[notSig,]$newFC.AB  = NA

## heatmap.data3 = t(cbind(heatmap.data$FC.A, heatmap.data$FC.AB))
heatmap.data3 = t(cbind(heatmap.data$newFC.A, heatmap.data$newFC.AB))
colnames(heatmap.data3) = heatmap.data$GeneName
rownames(heatmap.data3) = c(A, AB)
## rownames(heatmap.data3) = c(bquote(.(parse(text=A))), bquote(.(parse(text=AB))))
## rownames(heatmap.data3) = c(eval(parse(text=A)), eval(parse(text=AB)))
# rownames(heatmap.data3) = c(expression(paste(Delta, "cbrA pTH1227")), expression(paste(Delta, "cbrA pTH1227")))

dim(heatmap.data3)
head(heatmap.data3)

quartz() # new window
pdf(file="heatmap3.pdf")
heatmap.2(heatmap.data3, # values
          Rowv = NA, # no row clustering
          Colv = NA, # no column clustering
          scale = "row", # row-wise scaling
          na.color="black",
          ## col=greenred(75), # green to red heatmap
          ## col=blueyellow(75), # green to red heatmap
          col=colorRampPalette(c('blue', 'yellow')),
          labRow=c(A, AB), # MUST USE LABROW TO PROPERLY SHOW AN EXPRESSION!!!
          dendrogram="none", # no dendrogram
          key=TRUE, # add color key
          symkey=FALSE,
          density.info="none",
          trace="none",
          cexRow=0.6, # make row labels small
          cexCol=0.5 # make column labels small
          )

##############################################################################################
# sort by ClusterNumber and then by Rescued
##############################################################################################
## heatmap.data = heatmap.data[order(heatmap.data$ClusterNumber, -heatmap.data$Rescued),]
heatmap.data = heatmap.data[order(heatmap.data$ClusterNumber, -as.numeric(as.factor(heatmap.data$Rescued))),]

## heatmap.data3 = t(cbind(heatmap.data$FC.A, heatmap.data$FC.AB))
heatmap.data3 = t(cbind(heatmap.data$newFC.A, heatmap.data$newFC.AB))
colnames(heatmap.data3) = heatmap.data$GeneName
rownames(heatmap.data3) = c(A, AB)
## rownames(heatmap.data3) = c(bquote(.(parse(text=A))), bquote(.(parse(text=AB))))
## rownames(heatmap.data3) = c(eval(parse(text=A)), eval(parse(text=AB)))
# rownames(heatmap.data3) = c(expression(paste(Delta, "cbrA pTH1227")), expression(paste(Delta, "cbrA pTH1227")))

dim(heatmap.data3)
head(heatmap.data3)

quartz() # new window
heatmap.2(heatmap.data3, # values
          Rowv = NA, # no row clustering
          Colv = NA, # no column clustering
          scale = "row", # row-wise scaling
          na.color="black",
          ## col=greenred(75), # green to red heatmap
          ## col=blueyellow(75), # green to red heatmap
          col=colorRampPalette(c('blue', 'yellow')),
          labRow=c(A, AB), # MUST USE LABROW TO PROPERLY SHOW AN EXPRESSION!!!
          dendrogram="none", # no dendrogram
          key=TRUE, # add color key
          symkey=FALSE,
          density.info="none",
          trace="none",
          cexRow=0.6, # make row labels small
          cexCol=0.5 # make column labels small
          )

