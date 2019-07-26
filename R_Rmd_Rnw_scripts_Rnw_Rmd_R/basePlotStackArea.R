# -*- Mode:R; Coding:utf-8; fill-column:160 -*-

################################################################################################################################################################
# @file      stackedAreaChart.R
# @author    Mitch Richling <https://www.mitchr.me>
# @Copyright Copyright 2015 by Mitch Richling.  All rights reserved.
# @brief     Draw a filled in area graph with base R.@EOL
# @Keywords  R base area plot graphics
#
# In my opinion, this kind of graph is a almost always a poor way to communicate information -- ranking with a pie chart.  Individual time series (in one
# facet with colors or in multiple facets) along with a single cumulative sum is almost always more clear.  But, you ask for it... ;)
#

################################################################################################################################################################
# First we need some data

x <- 1:30
y <- cbind((30:1)*(1:30)*3,
           (1:30)*(1:30),
           (30:1)*(30:1)
           )

################################################################################################################################################################

par(mar=c(5,5,5,5))

plot(c(x,x), c(rowSums(y), y[,1]),                 # Setup the plot window and ranges
     col=NA,                                       # Do not actually draw anything.
     main='Area Plot', xlab='x', ylab='y')         # Need to do something reasonable for the y label

for(i in (dim(y)[2]):1) {
  yd <- rowSums(matrix(y[,1:i], ncol=i))           # Compute sum
  polygon(c(x, rev(x)),                            # Draw filled area
          c(yd, rep(par('usr')[3], dim(y)[1])),
          col=i,  border=NA)
  points(x, yd)                                    # Draw the points
}

box()                                              # Redraw the box because the polygons may have drawn over the axis box edges


