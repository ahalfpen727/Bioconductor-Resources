# -*- Mode:R; Coding:us-ascii-unix; fill-column:160 -*-

################################################################################################################################################################
##
# @file      basePlotPDF.R
# @author    Mitch Richling <https://www.mitchr.me>
# @Copyright Copyright 2015 by Mitch Richling.  All rights reserved.
# @brief     Get a base graphics plot into a PDF.@EOL
# @Keywords  base R pdf plot graphics
#

################################################################################################################################################################

pdf('basePlotPDF_out.pdf', width=10, height=8)
plot(1:100, rnorm(100), type='l')
dev.off()
system('xpdf -fullscreen basePlotPDF_out.pdf &')
