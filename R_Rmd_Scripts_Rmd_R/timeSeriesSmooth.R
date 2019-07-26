# -*- Mode:R; Coding:us-ascii-unix; fill-column:160 -*-

################################################################################################################################################################
# @file      timeSeriesSmooth.R
# @author    Mitch Richling <https://www.mitchr.me>
# @Copyright Copyright 2015 by Mitch Richling.  All rights reserved.
# @brief     Smoothing time series plots.@EOL
# @Keywords  base r smooth filter lowess
#
# Smoothers (especially ones found in CRAN packages) frequently have one or more of the following characteristics:
#   - They may require both x & y data
#   - They might only work when the x data is numeric -- they simply fail when given dates for example.
#   - They might coerce x data into numeric values
#   - They return a composite data type or list containing $x and $y values in the returned object
#   - They may change the number of observations
#
# Here we demonstrate how to use such smoothers with time series data where the x data are dates.
#

################################################################################################################################################################
# We need some data, so let's fake some. Note that we use POSIXct date objects -- I find them easy to work with as each date is an atomic object.
numPts <- 150
daDat <- data.frame(xpdate=as.POSIXct('2012-01-01')+(1:numPts)*(60*60*24),
                    y=rnorm(numPts, mean=0, sd=10)+(1:numPts)*.2)

################################################################################################################################################################
# Now we transform our POSIXct time stamps into integers so that they will work with any smoother.  The smoother we demonstrate, lowess, would do this step for
# us, but not all smoothers are so nice.
daDat$xpint <- as.integer(daDat$xpdate)

head(daDat)

#################################################################################################################################################################
# Now we smooth our data (we give it the x integers we constructed from our time data as well as the y data)
smoothedData <- lowess(daDat$xpint, daDat$y, f=.3)
head(smoothedData)

################################################################################################################################################################
# Now we put everything in a data.frame so we can graph it with ggplot.  Notice that we convert the integers we got from lowess at the same time.
allDat <- rbind(data.frame(smoother=rep('actual', numPts),
                           x=daDat$xpdate,
                           y=daDat$y),
                data.frame(smoother=rep('lowess', numPts),
                           x=as.POSIXct(smoothedData$x, origin='1970-01-01 00:00:00 UTC'),
                           y=smoothedData$y))

################################################################################################################################################################
# Plot them all
ggplot(allDat, aes(x=x, y=y, col=smoother)) +
  geom_line() +
  labs(title='Smoothing Time Series With A Generic Smoother')

