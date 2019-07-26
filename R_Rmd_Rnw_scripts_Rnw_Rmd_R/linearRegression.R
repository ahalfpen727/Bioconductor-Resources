# -*- Mode:R; Coding:utf-8; fill-column:160 -*-

################################################################################################################################################################
# @file      linearRegression.R
# @author    Mitch Richling <https://www.mitchr.me>
# @Copyright Copyright 2015 by Mitch Richling.  All rights reserved.
# @brief     Simple linear models in R.@EOL
# @Keywords  Simple linear models regression goodness fit normality confidence intervals prediction forecasting Shpiro Wilk r base
#
#            Topics:
#                1) Generating linear model data with random errors of different types (normal, uniform, positive uniform)
#                2) Fitting a linear model to some data
#                3) exploring the model (predictions, goodness of fit tests, confidence intervals),
#                4) graphical display (model, data, goodness fit, confidence intervals, & predictions).
#

################################################################################################################################################################
# Our data

daDat <- data.frame(x=c(1, 2, 3, 4, 7, 8, 9, 10, 12, 13, 17, 18, 20, 21, 22, 23, 24, 27, 28, 29, 30, 31, 32, 34, 35, 36, 37, 40, 41,
                        42, 43, 46, 47, 48, 49, 50, 51, 52, 53, 54, 57, 59, 60, 61, 62, 63, 67, 68, 70, 71, 72, 73, 74, 77, 78, 79,
                        80, 81, 82, 84, 85, 86, 87, 90, 91, 92, 93, 96, 97, 98, 99, 100)
                    )                        # x-data
daDat$f <- 10+daDat$x*1.5                    # Function
daDat$e <- rnorm(daDat$x, mean=0, sd=15)     # error                     ######## TRY THIS: Normally distributed errors
#daDat$e <- runif(daDat$x, min=-40, max=40)  # error                     ######## TRY THIS: Uniformly distributed errors
#daDat$e <- runif(daDat$x, min=0, max=80)    # error                     ######## TRY THIS: Uniformly distributed POSITIVE errors
daDat$y <- daDat$f + daDat$e                 # y-data

################################################################################################################################################################
# Compute the linear model with x being the independent variable and y being the dependent one

daFit <- lm(y ~ x, data=daDat)     # data= simplifies the formula and make future use of newdata= in 'predict' calls less painful.
summary(daFit)

################################################################################################################################################################
# plot the data with lines connecting points, draw the fitted line, mark the data points and fitted points with red dots and connect
# them with red lines.

ggplot(data=daDat, aes(x=x)) +
  geom_abline(intercept=coef(daFit)[1], slope=coef(daFit)[2], size=3, col='lightblue') +  # Fitted LINE
  geom_point(aes(y=y), col='red') +                                                       # Data points
  geom_line(aes(y=fitted(daFit)), col='blue', size=1) +                                   # Fitted SEGMENT
  geom_segment(aes(x=x, y=fitted(daFit), xend=x, yend=y), col='red') +                    # Error Bars
  labs(title='Data & Fitted Model', x='x', y='y')+                                        # Labels
  geom_rug(aes(y=y))                                                                      # Rugs for x and y data

################################################################################################################################################################
# Predict model values beyond our x-data, compute intervals, and make a nice graph showing it all

# Compute prediction intervals and confidence over the x-data and a an extended interval beyond the data...
newx <- data.frame(x=                                                      # Compute new x points on which intervals will be computed
                   seq(from=min(daDat$x),                                  # Start at first data point:   x_min
                       to=max(daDat$x)+diff(range(daDat$x)),               # End one "range" beyond last: x_max+(x_max-x_min)
                       length.out=100))                                    # Compute 100 points
tmpp <- predict(daFit, newdata=newx, interval="prediction")                # Compute prediction intervals
tmpc <- predict(daFit, newdata=newx, interval="confidence")                # Compute confidence intervals
newx <- data.frame(x=newx$x, fit=tmpc[,'fit'],                             # Put everything in a data.frame
                   plwr=tmpp[,'lwr'], pupr=tmpp[,'upr'],
                   clwr=tmpc[,'lwr'], cupr=tmpc[,'upr'])

# Plot everything
ggplot() +
  geom_ribbon(data=newx, aes(x=x, ymin=plwr, ymax=pupr), fill='yellow') +                      # prediction intervals
  geom_ribbon(data=newx, aes(x=x, ymin=clwr, ymax=cupr), fill='gold') +                        # confidence intervals
  geom_abline(intercept=coef(daFit)[1], slope=coef(daFit)[2], size=3, col='lightblue') +       # Fitted LINE
  geom_line(data=newx, aes(x=x, y=fit), col='red', size=1) +                                   # Prediction SEGMENT
  geom_line(data=daDat, aes(x=x, y=fitted(daFit)), col='blue', size=1) +                       # Fitted SEGMENT
  geom_point(data=daDat, aes(x=x, y=y), col='red') +                                           # Data points
  geom_segment(data=daDat, aes(x=x, y=fitted(daFit), xend=x, yend=y), col='red') +             # Error Bars
  labs(title='Data & Fitted Model', x='x', y='y')

# Clean up temporary variables.
rm(newx, tmpp, tmpc)

################################################################################################################################################################
# Check the residuals for normality with the Shpiro-Wilk test
# If p-value<0.1, then we reject the hypothesis that the residuals are normally distributed.

shapiro.test(residuals(daFit))    # NOTE: residuals(daFit) works on more "fit-like" objects than daFit$residuals

################################################################################################################################################################
# Check the residuals for normality with a Q-Q plot

slope=(quantile(residuals(daFit),p=.75)-quantile(residuals(daFit),.25))/(qnorm(.75)-qnorm(.25))
intercept = quantile(residuals(daFit),.25) - slope*qnorm(.25)
ggplot() +
  geom_point(aes(sample=residuals(daFit)), stat="qq", distribution=qnorm) +   # Draw a Q-Q plot against a normal
  geom_abline(intercept=intercept,slope=slope)                                # Draw the Q-Q plot *LINE* against a normal

# Clean up temporary variables.
rm(slope, intercept)

################################################################################################################################################################
# Check the residuals for normality with a scatter plot (lines at mean, +sd, -sd, 2*sd, * -2*sd) and a histogram

daScat <- ggplot(data=daDat, aes(x=x)) +
  geom_point(aes(y=residuals(daFit)), col='red') +                               # Residuals
  geom_segment(aes(x=x, y=0, xend=x, yend=residuals(daFit)), col='red') +        # Residuals Bars
  labs(x='x', y='Residuals') +                                                   # Labels
  coord_cartesian(ylim = range(1.1*residuals(daFit)))

daHist <- ggplot() +
  geom_histogram(aes(x=residuals(daFit)), col='black', fill='red')+
  theme(axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),           
        axis.title.y=element_blank()
        ) +
  coord_flip(xlim = range(1.1*residuals(daFit)))

grid.arrange(daScat, daHist, ncol=2, nrow=1, widths=c(5, 2))

# Clean up temporary variables.
rm(daScat, daHist)

################################################################################################################################################################
# Clean up

rm(daDat, daFit)
