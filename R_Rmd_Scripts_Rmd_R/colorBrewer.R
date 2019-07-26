# -*- Mode:R; Coding:us-ascii-unix; fill-column:160 -*-

################################################################################################################################################################
##
# @file      colorBrewer.R
# @author    Mitch Richling <https://www.mitchr.me>
# @Copyright Copyright 2015 by Mitch Richling.  All rights reserved.
# @brief     Using color brewer colors in R.@EOL
# @Keywords  colorbrewer package cran brewer colors R
#
# See ggplot2.R for how to use colorbrewer with ggplot2.
#            

################################################################################################################################################################
# Load the library

library(RColorBrewer)

################################################################################################################################################################
# The "Sequential palettes"

display.brewer.all(type="div")

################################################################################################################################################################
# The "Diverging palettes"

display.brewer.all(type="seq")

################################################################################################################################################################
# The "Qualitative palettes"

display.brewer.all(type="qual")

################################################################################################################################################################
# Using color brewer with base graphics barplot

barplot(1:10, col=mypalette<-brewer.pal(10,"Set3"))

################################################################################################################################################################
# Using color brewer with base graphics image

x <- seq(-pi, pi, len = 300)
y <- x
r <- sqrt(outer(x^2, y^2, "+"))
z <- sin(r^2)
image(z, col=brewer.pal(11,"Spectral"))

################################################################################################################################################################
# Sometimes you need more colors.  Here is one way to expand the number of colors

image(z, col=colorRampPalette(brewer.pal(11,"Spectral"))(100))
