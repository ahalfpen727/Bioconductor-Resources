# -*- Mode:R; Coding:utf-8; fill-column:160 -*-

################################################################################################################################################################
# @file      plotAnno.R
# @author    Mitch Richling <https://www.mitchr.me>
# @Copyright Copyright 2015 by Mitch Richling.  All rights reserved.
# @brief     How I annotate the lower margin of plots.@EOL
# @Keywords  base R annotations plot graphics
#
# I like to annotate the bottom of my plots with authorship information, date created, assumptions, and references/links to the original data sources used.
#

################################################################################################################################################################
# Set a large enough bottom margins for our annotations (I like to have 4 usable lines)

par(mar=c(8,4,4,2)+0.1)                          

################################################################################################################################################################
# Plot some random numbers using a semi-transparent 2D scatter-plot technique to show point density

plot(rnorm(2000),rnorm(2000, sd=2),              
     xlab='x', ylab='y', main='Random Points',   # Provide a title and axis labels
     pch=20,                                     # pch 20 (a filled circle) has broad support across PDF readers
     cex=3,                                      # Bigger plotting points work better with transparency
     col=rgb(0,0,0,.1))                          # Set the transparency to 1/10 because we have lots of points

################################################################################################################################################################
# Annotations in the bottom margin (side=1) justified on the left side of the page (adj=0).

mtext('Some left side text on line  4',
      col='blue',                                # I like blue.  Contrasts with titles and labels, but isn't garish
      cex=0.75,                                  # I like to shrink the text a bit to make it less intrusive
      side=1,                                    # side=1 => annotation in the BOTTOM margin
      adj=0,                                     # adj=0  => left justify the text
      line=4)                                    # line=4 => put the text on the fourth margin line
mtext('Some left side text on line  5', col='blue', cex=0.75, side=1, adj=0, line=5)
mtext('Some left side text on line  6', col='blue', cex=0.75, side=1, adj=0, line=6)
mtext('Some left side text on line  7', col='blue', cex=0.75, side=1, adj=0, line=7)

################################################################################################################################################################
# Annotations in the bottom margin (side=1) justified on the right side of the page (adj=1).

mtext('Data Source: Unpublished made up around 1998 by Mitch',        col='blue', cex=0.75, side=1, adj=1, line=5)
mtext('Graph By: Mitch Richling <https://www.mitchr.me> (2003-02-17)', col='blue', cex=0.75, side=1, adj=1, line=6)
