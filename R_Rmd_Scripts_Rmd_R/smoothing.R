# -*- Mode:R; Coding:us-ascii-unix; fill-column:160 -*-

################################################################################################################################################################
# @file      dataSmooth.R
# @author    Mitch Richling <https://www.mitchr.me>
# @Copyright Copyright 2015 by Mitch Richling.  All rights reserved.
# @brief     Smoothing data.@EOL
# @Keywords  r base smooth filter lowess running median mean cleveland tukey legend
#
# We demonstrate several smoothers built into R:
#  - Cleveland's lowess procedure (lowess)
#  - Running averages (stats::filter)
#  - Running median (runmed)
#  - The Tukey 3RS3R smoother (smooth)
#  - Friedman's Super Smoother (supsmu)
# We also demonstrate several lowess smoothings with various f parameter values
#

################################################################################################################################################################
# Create some data
numPts <- 150
daDat <- data.frame(x=1:numPts,
                    y=rnorm(numPts, mean=0, sd=10)+(1:numPts)*.2)

################################################################################################################################################################
# Running Median
smDatRunMed   <- runmed(daDat$y, 15)

################################################################################################################################################################
# Two sided, running average -- Should have an odd length convolution
smDatRunAvg2  <- stats::filter(daDat$y, rep(1/15,15), method='convolution', sides=2)

################################################################################################################################################################
# One sided, running average -- must have an odd length convolution
smDatRunAvg1  <- stats::filter(daDat$y, rep(1/7,7), method='convolution', sides=1)

################################################################################################################################################################
# Tukey Smoother -- other 'kind' values may be used
smDatTuk3RS3R <- smooth(daDat$y, kind = "3RS3R")

################################################################################################################################################################
# Cleveland lowess smoother
smDatLowess0d10   <- lowess(daDat$x, daDat$y, f=.1)

################################################################################################################################################################
# Friedman's SuperSmoother
smDatSupSmu   <- supsmu(daDat$x, daDat$y, span=.1)

################################################################################################################################################################
# Put it all in a data.frame for ggplot
allDat <- rbind(data.frame(smoother=rep('actual',   numPts), x=daDat$x,           y=daDat$y),
                data.frame(smoother=rep('runMed',   numPts), x=daDat$x,           y=smDatRunMed),
                data.frame(smoother=rep('runAvg2',  numPts), x=daDat$x,           y=smDatRunAvg2),
                data.frame(smoother=rep('runAvg1',  numPts), x=daDat$x,           y=smDatRunAvg1),
                data.frame(smoother=rep('tuk3RS3R', numPts), x=daDat$x,           y=as.vector(smDatTuk3RS3R)),
                data.frame(smoother=rep('lowess',   numPts), x=smDatLowess0d10$x, y=smDatLowess0d10$y),
                data.frame(smoother=rep('supsmu',   numPts), x=smDatSupSmu$x,     y=smDatSupSmu$y)
                )

################################################################################################################################################################
# Plot them all
ggplot(allDat, aes(x=x, y=y, col=smoother)) +
  geom_line(data=subset(allDat,  smoother!='actual')) +
  geom_point(data=subset(allDat, smoother=='actual'))

################################################################################################################################################################
# Now we compute several lowess smoothings with different f values
smDatLowess0d05   <- lowess(daDat$x, daDat$y, f=.05)
smDatLowess0d20   <- lowess(daDat$x, daDat$y, f=.20)
smDatLowess0d50   <- lowess(daDat$x, daDat$y, f=.50)

################################################################################################################################################################
# Put all of our lowess curves into a data.frame with the original data so we can plot it all with ggplot
allDat <- rbind(data.frame(smoother=rep('actual',       numPts), x=daDat$x,           y=daDat$y),
                data.frame(smoother=rep('lowess0d05',   numPts), x=smDatLowess0d05$x, y=smDatLowess0d05$y),
                data.frame(smoother=rep('lowess0d10',   numPts), x=smDatLowess0d10$x, y=smDatLowess0d10$y),
                data.frame(smoother=rep('lowess0d20',   numPts), x=smDatLowess0d20$x, y=smDatLowess0d20$y),
                data.frame(smoother=rep('lowess0d50',   numPts), x=smDatLowess0d50$x, y=smDatLowess0d50$y)
                )

################################################################################################################################################################
# Plot them all
ggplot(allDat, aes(x=x, y=y, col=smoother)) +
  geom_line(data=subset(allDat,  smoother!='actual')) +
  geom_point(data=subset(allDat, smoother=='actual')) +
  labs(title='Lowess At Various f Values')

