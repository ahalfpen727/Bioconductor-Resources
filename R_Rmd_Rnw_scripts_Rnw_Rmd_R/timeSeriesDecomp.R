# -*- Mode:R; Coding:utf-8; fill-column:160 -*-

################################################################################################################################################################
# @file      timeSeriesDecomp.R
# @author    Mitch Richling <https://www.mitchr.me>
# @Copyright Copyright 2015 by Mitch Richling.  All rights reserved.
# @brief     How to decompose a time series with built in R commands.@EOL
# @Keywords  r base timeseries decompose forecast stats::ts stats::decompose stats::ts.plot stats::arima
#
# This example 1) constructs a time series with a linear trend, a 7 day seasonal component, and a normally distributed random component, 2) decomposes that time
# series, and 3) plots the decomposition along with the original components.  Allows one to explore the implications of a positive seasonal component.
#

################################################################################################################################################################
# First we create some data

daData          <- data.frame(date=as.POSIXct('2012-01-01')+(1:365)*(60*60*24))
daData$x        <- as.numeric(daData$date)
daData$x        <- (daData$x-min(daData$x))/(60*60*24)
daData$trend    <- daData$x/50
daData$seasonal <- sin(pi*daData$x/3.5)                         ######## TRY THIS: equal positive and negative components
#daData$seasonal <- abs(1+sin(pi*daData$x/3.5))                 ######## TRY THIS: positive seasonal component
daData$random   <- rnorm(daData$x, sd=.25)
daData$val      <- daData$trend+daData$seasonal+daData$random

################################################################################################################################################################
# Now we construct a TS from daData$val with a frequency of 7 days, and then we decompose into an additive seasonal model

daDataSeries <- ts(daData$val, frequency=7)
daDataDecomp <- decompose(daDataSeries, type='add')

daDataDecompDF <- data.frame(date=daData$date, val=daData$val,
                             trend=daDataDecomp$trend, seasonal=daDataDecomp$seasonal, random=daDataDecomp$random)
daDataDecompDF <- melt(daDataDecompDF, id="date")


################################################################################################################################################################
# Plot things with base

plot(daDataDecomp)

################################################################################################################################################################
# Plot things with lattice

xyplot(daData$val + daDataDecomp$trend + daDataDecomp$seasonal + daDataDecomp$random ~ daData$date, type='l', outer=TRUE, horizontal=FALSE, layout=c(1,4))

################################################################################################################################################################
# Plot things with ggplot2

ggplot(data=daDataDecompDF, aes(x=date)) +
    geom_line(aes(y=value))  +
    facet_grid(variable ~ ., scales = "free")

################################################################################################################################################################
# Plot things with base the KRAZY way

par(mfcol=c(4,1))
par(mar=c(.5,2.5,.5,.5))
plot(daData$date, daData$val, type='l', ylab='', xaxt='n')
text(mean(par('usr')[1:2]), par('usr')[4], 'Value', pos=1, cex=3, col='blue')
par(mar=c(.5,2.5,0,.5))
plot(as.POSIXct('2012-01-01'), 0,
     xlim=range(daData$date), ylim=range(c(daDataDecomp$trend, daData$trend), na.rm=TRUE),
     col=NA, ylab='', xaxt='n')
points(daData$date, daDataDecomp$trend, type='l', xaxt='n')
points(daData$date, daData$trend,       type='l', col='red')
text(mean(par('usr')[1:2]), par('usr')[4], 'Trend', pos=1, cex=3, col='blue')
plot(as.POSIXct('2012-01-01'), 0,
     xlim=range(daData$date), ylim=2*range(c(daDataDecomp$seasonal, daData$seasonal), na.rm=TRUE),
     col=NA, ylab='', xaxt='n')
points(daData$date, daDataDecomp$seasonal, type='l', xaxt='n')
points(daData$date, daData$seasonal,       type='l', col='red')
text(mean(par('usr')[1:2]), par('usr')[4], 'Seasonal', pos=1, cex=3, col='blue')
par(mar=c(2.5,2.5,0,.5))
plot(as.POSIXct('2012-01-01'), 0,
     xlim=range(daData$date), ylim=range(c(daDataDecomp$random, daData$random), na.rm=TRUE),
     col=NA, xlab='', ylab='')
points(daData$date, daData$random,       type='p', col='red', pch=20)
points(daData$date, daDataDecomp$random, type='l', xaxt='n')
text(mean(par('usr')[1:2]), par('usr')[4], 'Random', pos=1, cex=3, col='blue')
    
################################################################################################################################################################
# Fit an arima model, and predict the future

fit <- arima(daDataSeries, order=c(5,0,0), seasonal=list(order=c(2,1,0), period=7))
fore <- predict(fit, n.ahead=7*5)

################################################################################################################################################################
# # error bounds at 95% confidence level

U <- fore$pred + 2*fore$se
L <- fore$pred - 2*fore$se

################################################################################################################################################################
# Plot our prediction

par(mfcol=c(1,1))
par(mar=c(5,5,5,5))
ts.plot(daDataSeries, fore$pred, U, L, col=c(1,2,4,4), lty = c(1,1,2,2))
legend("topleft", c("Actual", "Forecast", "Error Bounds (95% Confidence)"), col=c(1,2,4), lty=c(1,1,2))
