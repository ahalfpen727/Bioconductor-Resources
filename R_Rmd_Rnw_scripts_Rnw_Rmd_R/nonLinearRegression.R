# -*- Mode:R; Coding:utf-8; fill-column:160 -*-

################################################################################################################################################################
# @file      nonLinearRegression.R
# @author    Mitch Richling <https://www.mitchr.me>
# @Copyright Copyright 2015 by Mitch Richling.  All rights reserved.
# @brief     @EOL
# @Keywords  r base nonlinear regression multiple
#
# The topic of this example is non-linear regression; however, as most of what one might do with a non-linear model is identical to what one might do with a
# linear model, the focus is really on how to fit the data.  For the other stuff see linearRegression.R and polynomialRegression.R.  In addition to fitting the
# model, we demonstrate how to generate random data that follows a predetermined model and then how to compare our fitted model to that original model -- just
# for fun.
#


################################################################################################################################################################
## Make up some data: y=sin(x)+e where e is random (normal or uniform)

numPts   <- 100
daDat    <- data.frame(x=1:numPts/(numPts/20*pi))    # x-data
daDat$f  <- sin(daDat$x)                  # Function
daDat$e  <- rnorm(daDat$x, mean=0, sd=.5) # Identically distributed Normal Error
#daDat$e <- runif(daDat$x, -1, 1)         # Identically distributed Uniform Error
#daDat$e <- runif(daDat$x, 0, 1)          # Identically distributed, but asymmetric, Uniform Error
#daDat$e <- daDat$x*rnorm(daDat$x, sd=.5) # Non-Identically distributed Normal Error
daDat$y  <- daDat$f + daDat$e             # y-data

################################################################################################################################################################
# Compute the non-linear model with independent variable of x, dependent variable of y, and function of a*sin(b*x+c)+d

daFit   <- nls(y~a*sin(b*x+c)+d,           # Specify our model formula
               data=daDat,                 # Always use data= to simplify the formula argument, and future predict calls
               start=list(a=1,b=1,c=0,d=0) # Initial conditions (we set them to the true model values)
               )

summary(daFit)

################################################################################################################################################################
# This plot is not terribly useful for the practicing data modeler as one wouldn't be modeling the data in the first place if the true model from which the data
# was generated available!  That said, it is an interesting way to explore how pushing the envelope of the various theoretical requirements impacts the accuracy
# of the fit (for example, try adding one sided, positive errors or non-normal ones).


daDat$fit <- fitted(daFit)                                # Put the fitted valeus in the data.frame so ggplot is simple to use

ggplot(daDat) +
  geom_ribbon(aes(x=x,
                  ymin=pmin(daDat$f, daDat$fit),
                  ymax=pmax(daDat$f, daDat$fit)),
              fill='pink') +
  geom_line(aes(x=x, y=f,   col='function')) +           # Note col is an aes -- so we get a legend
  geom_point(aes(x=x, y=y,  col='data')) +
  geom_line(aes(x=x, y=fit, col='fit')) +
  labs(title='Fit vs Actual',
       x='x',
       y='y') +
  scale_color_manual(values=c("black",     "red", "blue"))

################################################################################################################################################################
# How we might do it with base R

par(mar=c(8,4,4,2)+0.1)
plot(0, 0, col='NA',                                                                       # Setup the plot window
     xlim=range(daDat$x),                                                                    # Set the x-range
     ylim=range(daDat$y, daDat$f, fitted(daFit)),                                            # Set the y-range
     main='Fit vs Actual', xlab='x', ylab='y')                                               # Set the plot titles
polygon(c(daDat$x, rev(daDat$x)), c(daDat$f, rev(fitted(daFit))), col='pink', border=NA)   # Fill in area between fit and Actual
points(daDat$x, daDat$y, type='p')                                                         # Draw data (with random noise)
points(daDat$x, daDat$f,       type='l', col='green')                                      # Draw Actual
points(daDat$x, fitted(daFit), type='l', col='red')                                        # Draw Fit





