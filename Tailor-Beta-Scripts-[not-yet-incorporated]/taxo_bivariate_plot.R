#!/bin/Rscript
# ***************************************************************
# Name:      generateBivariatePlot.R
# Purpose:   Generates bivariate plots with histograms on the diagonals,
#            scatter plots with smooth curves below the diagonals and
#            correlations with significance levels above diagonals.
#       	   Datafile has the following organization:
#     			        Var_1	Var_2	Var_3 .. Var_R
#		         Sample_1
#		         Sample_2
#		         Sample_3
#		         ...
#		         Sample_N
#
#            Currently, it has support for three correlation measures:
#                 -Pearson
#                 -Spearman
#                 -Kendall
#	           The significance levels are as follows:
#	   	            if pvalue< 0.05 display: *
# 		            if pvalue <= 0.01 display: **
# 		            if pvalue <= 0.001 display: ***
#	           Moreover, the variables are reordered in the plots with any two consecutive
#            variables being most similar
# Version:   0.4
# History:   0.1 to 0.2 fixed the bug that produces exta spaces in the files by including sep="" in paste command
#		         0.2 to 0.3 included Q Mode and R Mode
#		         0.3 to 0.4 bug in the panel.cor. Method was not assigned properly
#		         0.4 to 0.5 transparent png bug fixed
# Authors:   Umer Zeeshan Ijaz (Umer.Ijaz@glasgow.ac.uk)
#                 http://userweb.eng.gla.ac.uk/umer.ijaz
#            Christopher Quince (Christopher.Quince@glasgow.ac.uk)
#                 http://userweb.eng.gla.ac.uk/christopher.quince
# Created:   2012-09-10
# License:   Copyright (c) 2012 Computational Microbial Genomics Group, University of Glasgow, UK
#
#            This program is free software: you can redistribute it and/or modify
#            it under the terms of the GNU General Public License as published by
#            the Free Software Foundation, either version 3 of the License, or
#            (at your option) any later version.
#
#            This program is distributed in the hope that it will be useful,
#            but WITHOUT ANY WARRANTY; without even the implied warranty of
#            MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#            GNU General Public License for more details.
#
#            You should have received a copy of the GNU General Public License
#            along with this program.  If not, see <http://www.gnu.org/licenses/>.
# **************************************************************/

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("gclus"))

#specify desired options in a list
option_list <- list(
   make_option("--ifile", action="store",default=NULL, help="CSV file"),
   make_option("--opath", action="store",default=NULL, help="Output path"),
   make_option("--fsize", action="store", default="1.2", help="Font size [default %default]"),
   make_option("--width", type="integer",default=800, help="Width of jpeg files [default %default]"),
   make_option("--height", type="integer",default=800, help="Height of jpeg files [default %default]"),
   make_option("--correlation", type="integer",default=1, help="Correlation to use: 1=pearson, 2=spearman, 3=kendall [default %default]"),
   make_option("--rmode", action="store_true",default=FALSE, help="Mode: TRUE=R mode, FALSE=Q mode [default %default]")
)

#get command line options
opt<-parse_args(OptionParser(usage="%prog [options] file", option_list=option_list))

if(is.null(opt$ifile))
   quit()

#Import data now
data <- read.csv(opt$ifile,header=TRUE,row.names=1)

#use complete.cases function to filter out those rows that have NA entry in Env
data<-data[complete.cases(data),]

#Transpose the data if R Mode
if(opt$rmode){
   data=t(data)
}
method<-switch(opt$correlation,"pearson","spearman","kendall")
#Generate panels for the graph
panel.cor <- function(x, y, method=method, digits=3, cex.cor=as.numeric(opt$fsize))
{
   usr <- par("usr"); on.exit(par(usr))
   par(usr = c(0, 1, 0, 1))
   r <- cor(x, y, method=method)
   ra <- cor.test(x, y, method=method)$p.value
   txt <- round(r, digits)
   sig <- 1
   prefix <- ""
   if(ra <= 0.1) prefix <- "."
   if(ra <= 0.05) prefix <- "*"
   if(ra <= 0.01) prefix <- "**"
   if(ra <= 0.001) prefix <- "***"
   if(ra <= 0.001) sig <- 2
   color <- 2
   if(r < 0) color <- 4
   txt <- paste(txt, prefix, sep="\n")
   text(0.5, 0.5, txt, cex = cex.cor, font=sig, col=color)
}

# Put histograms on the diagonal
panel.hist <- function(x, ...)
{
   usr <- par("usr"); on.exit(par(usr))
   par(usr = c(usr[1:2], 0, 1.5) )
   h <- hist(x, plot = FALSE)
   breaks <- h$breaks; nB <- length(breaks)
   y <- h$counts; y <- y/max(y)
   rect(breaks[-nB], 0, breaks[-1], y, col="yellow", ...)
}

#Find the correlations
data.correlation <- cor(data,method=method)

# Reorder the data prior to plotting
data.o <- order.single(data.correlation)

# Plot the data
jpeg(filename = paste(opt$opath,substr(basename(opt$ifile),1,nchar(basename(opt$ifile))-4),"_BP.jpg",sep=""),width = as.numeric(opt$width), height = as.numeric(opt$height), quality=100)

op <- par(mfrow=c(1,1), pty="s",cex.axis=(as.numeric(opt$fsize)*0.6))
pairs(data[,data.o], lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist,cex.labels=as.numeric(opt$fsize),cex=as.numeric(opt$fsize),method=method)
par(op)
dev.off()


