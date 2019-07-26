# -*- Mode:R; Coding:utf-8; fill-column:160 -*-

################################################################################################################################################################
# @file      polynomialRegression.R
# @author    Mitch Richling <https://www.mitchr.me>
# @Copyright Copyright 2015 by Mitch Richling.  All rights reserved.
# @brief     Polynomial regression and graphics.@EOL
# @Keywords  polynomial regression multiple base r cran package ggplot2
#
# This example highlights some of the particulars of doing "simple" polynomial regression R when compared to simple linear regression:
#  1) Automating multiple fits
#  2) Specifying the model formula for the fit,
#  3) Selecting the model degree using ANOVA, and graphically with simultaneous model plots
#  4) Visually exploring confidence and prediction intervals vs model degree with ggplot2.
#

################################################################################################################################################################
# Our data.  

# TRUE => fixed data, FALSE => randomly generated data.
if(TRUE) {
  daDat   <- data.frame(x=c(  0.0000000,  0.4444444,  0.8888889, 1.3333333, 1.7777778,  2.2222222,   2.6666667,  3.1111111,  3.5555556, 4.0000000),
                        f=c( -6.0000000, -2.2085048, -0.2606310, 0.3703704, 0.2112483, -0.2112483,  -0.3703704,  0.2606310,  2.2085048, 6.0000000),
                        e=c( -0.5747532, -0.1930761,  1.7809493, 1.2608157, 0.8505483,  1.23321618, -1.5896486, -0.0399569, -1.7746238, 0.0406907),
                        y=c( -6.5747532, -2.4015809,  1.5203184, 1.6311861, 1.0617966,  1.0219679,  -1.9600190,  0.2206741,  0.4338810, 6.0406907))
} else {
  daDat   <- data.frame(x=seq(0, 4, length.out=20))  # x-data
  daDat$f <- with(daDat, x^3-6*x^2+11*x-6)           # Function
  daDat$e <- rnorm(daDat$x)                          # Error
  daDat$y <- daDat$f + daDat$e                       # y-data
}

head(daDat, 20)

################################################################################################################################################################
# We fit four models using progressively higher degree polynomials (degree 1 up to maxFdeg)

maxFdeg <- 4                                            ######## TRY THIS: Value of 6 vs. 4
if(maxFdeg == 4) {
  ## For illustrative purposes we demonstrate hand coded formulas for the maxFdeg==4 case
  daFits <- list(lm(y ~ x,                    data=daDat),  # The degree 1 case is simple linear regression
                 lm(y ~ x + I(x^2),           data=daDat),  # Note the use of the "I" function here
                 lm(y ~ x + I(x^2) + I(x^3),  data=daDat),
                 lm(y ~ poly(x, 4, raw=TRUE), data=daDat))  # The formula gets a bit long by deg>4, so we use the "poly" function
} else {
  ## This is how one generic fitting formulas can be constructed
  daFits <- lapply(1:maxFdeg, function (i) lm(y ~ poly(x, i, raw=TRUE), data=daDat))
}

daFits

################################################################################################################################################################
# Compute model values at 100 points between the min and max x values

newx  <- data.frame(x=seq(min(daDat$x), max(daDat$x), length.out=100))
newy  <- NULL
for(daFitDeg in 1:maxFdeg)
  newy  <- rbind(newy, data.frame(x=newx,
                                  degree=rep(daFitDeg, length(newx)),
                                  y=predict(daFits[[daFitDeg]], newdata=newx)))
newy$degree <- factor(newy$degree)
  
ggplot() +
  geom_line(data=newy, aes(x=x, y=y, col=degree)) +
  geom_line(data=daDat, aes(x=x, y=y), lwd=2)

################################################################################################################################################################
# Use ANOVA to determine which fit seems best

do.call(anova, daFits)

################################################################################################################################################################
# Compute prediction intervals and confidence intervals.  Draw the result.

newx  <- data.frame(x=seq(from=min(daDat$x),                                    
                          to=max(daDat$x)+diff(range(daDat$x))/4, # go out only 1/4 as far as we did in linearRegression.R
                          length.out=100))
newy  <- NULL
for(daFitDeg in 1:maxFdeg) {
  tmpp <- predict(daFits[[daFitDeg]], newdata=newx, interval="prediction")
  tmpc <- predict(daFits[[daFitDeg]], newdata=newx, interval="confidence")
  newy  <- rbind(newy, data.frame(x=newx, degree=rep(daFitDeg, length(newx)),
                                  fit=tmpp[,'fit'], pLow=tmpp[,'lwr'], pUp=tmpp[,'upr'], cLow=tmpc[,'lwr'], cUp=tmpc[,'upr']))
}
newy$degree <- factor(paste('degree', newy$degree))

ggplot(newy, aes(x=x, y=fit, group=degree)) +
  facet_wrap(~degree, ncol=2) +
  geom_ribbon(aes(ymin=pLow, ymax=pUp), alpha=.5, fill='pink', col='red') +
  geom_ribbon(aes(ymin=cLow, ymax=cUp), alpha=.5, fill='red', col='pink') +
  geom_line()

