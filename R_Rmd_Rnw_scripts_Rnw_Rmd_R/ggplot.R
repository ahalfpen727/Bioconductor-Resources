# -*- Mode:R; Coding:us-ascii-unix; fill-column:160 -*-

################################################################################################################################################################
##
# @file      ggplot.R
# @author    Mitch Richling <https://www.mitchr.me>
# @Copyright Copyright 2015 by Mitch Richling.  All rights reserved.
# @brief     ggplot2 examples.@EOL
# @Keywords  ggplot2 package r cran examples
#

################################################################################################################################################################
# Load libraries.  ggplot2 is the graphics library, and gridExtra helps us lay things out.

library(ggplot2)
library(gridExtra)
library(reshape2)
library(scales)

################################################################################################################################################################
# Create some data sets we can graph...

numPnt <- 1000
numCat <- 4
someData <- data.frame(xi=(1:numPnt)/numPnt,
                       xi2=(1:numPnt)/numPnt+2,
                       xd=seq(as.POSIXct('2012-01-01'), by='day', length.out=numPnt),
                       rn=runif(numPnt, -.6, .6),
                       m=rep(1:numCat,(numPnt/numCat)))
someData$rm <- someData$rn*someData$m
someData$y  <- someData$m*(someData$xi-0.5)+0.5+someData$rm
someData$y2 <- 2.5*(someData$xi-0.5)+0.5
someData$k1 <- paste('K1F', someData$m, sep='_')

someMoreData <- read.table(header=T, text='
 quarter factors   value emin emax
      Q1 factor1      15   14   16
      Q2 factor1      25   22   29
      Q3 factor1       5    1    6
      Q4 factor1      10    8   11
      Q1 factor2      25   22   23
      Q2 factor2      20   19   21
      Q3 factor2      25   20   26
      Q4 factor2      10   10   11
 ')

someLessData <- subset(someMoreData, someMoreData$factors=='factor1')
someLessData$factors <- factor(someLessData$factors)

twoDimData  <-expand.grid(x=seq(from=-4, to=4, length.out=30), y=seq(from=-3, to=5, length.out=30))
twoDimData$z<-sin(twoDimData$x*twoDimData$x+twoDimData$y*twoDimData$y)
twoDimData$p<-sign(twoDimData$z)

# Transform a matrix into a data.frame
volcanoDF   <- stack(as.data.frame(volcano))
volcanoDF$x <- rep(seq_len(nrow(volcano)), ncol(volcano))
volcanoDF$y <- rep(seq_len(ncol(volcano)), rep(nrow(volcano), ncol(volcano)))

twoDimNormal <- data.frame(x=rnorm(10000, sd=2.0),
                           y=rnorm(10000))
twoDimNormal$y <- twoDimNormal$y+twoDimNormal$x*.75
  
################################################################################################################################################################
# Simple "scatterplot" of ONE variable.  i.e. do something like 'plot(someData$y)'
# Alternate with qplot: qplot(seq_along(someData$y), someData$y)

ggplot(data=someData, aes(y=y, x=seq_along(y))) + geom_point()

################################################################################################################################################################
# Scatterplot with a couple vectors.  i.e. do something like 'plot(x, y)'
# Alternate with qplot:  qplot(rnorm(100), y=rnorm(100))

ggplot() + geom_point(aes(x=rnorm(100), y=rnorm(100)))

################################################################################################################################################################
# And we can use other ascetics as well

ggplot() + geom_point(aes(x=rnorm(50),
                          y=rnorm(50),
                          shape=factor(as.integer(runif(50, min=0, max=5))),
                          size=rnorm(50),
                          col=rnorm(50)))

################################################################################################################################################################
# An x-y scatter plot from data in a data.frame -- ggplot is generally much less difficult to use when everything is in a data frame.

ggplot(data=someData, aes(x=xi, y=y)) +
  geom_point(col='red')                 ## Just make the points in the graph RED!! ;)
  #geom_point(aes(col='red')) +         ## This is very different -- see next example

################################################################################################################################################################
# An x-y scatter plot from data in a data.frame broken down by a factor (displayed as color).  
# Extra care is take to select date breaks for the x-axis and a date format for the tick labels.
# We also customize the legend with a title and different labels.
# Alternate with qplot: qplot(data=someData,x=xd, y=y, color=k1, main='Main Title', xlab='X Title', ylab='Y Title', geom='point')
        
ggplot(data=someData, aes(x=xd, y=y, col=k1)) +
  scale_x_datetime(breaks       = scales::date_breaks("1 year"),           ## Set major break lines to 1 year
                   minor_breaks = scales::date_breaks("1 month"),          ## Set minor break lines to 1 year
                   labels       = scales::date_format("%Y")) +             ## Set the date format
  geom_point(size=3, pch=21) +                                         ## BIG Circles for points
  #geom_line() +                                                       ## Use this to add lines
  labs(title='Main Title',
       x='X Title',
       y='Y Title') +
  scale_colour_discrete(name ="Legend Title\nLine 2",                  ## Control title and labels in the color legend
                        breaks=c("K1F_4", "K1F_1", "K1F_2", "K1F_3"),
                        labels=c("k1f_4", "k1f_1", "k1f_2", "k1f_3"))

################################################################################################################################################################
# Controlling the axes and labels

ggplot(data=someData, aes(x=xi, y=y, col=k1)) + geom_point() +
  #theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())  ## No grid lines, but have tick marks
  #theme(axis.ticks = element_blank(), axis.text.x = element_blank())         ## x-axis: grid lines, but no tick marks (or labels)
  scale_y_continuous(breaks=c(1, 1.5, 2.5, 5, 6,8))                           ## y-axis: Set location for ticks and grid lines
  #scale_y_continuous(breaks=NULL)                                            ## Zap y tick and grid lines
  #scale_x_continuous(breaks=NULL)                                            ## Zap x tick and grid lines
  #scale_x_continuous(limits=c(-1, 2))                                        ## Remove data outside the limits
  #xlim(NA, 2)                                                                ## Like above.  Note you can use NA if you don't want to set a limit
  #coord_cartesian(xlim = c(-1, 2))                                           ## Simply set the visible region.

################################################################################################################################################################
# Adding a 'ribbon' of color to a plot (useful for things like confidence intervals and ranges around a smoother). Also demo some crazy title and legend stuff.

ggplot(data=someData, aes(x=xd)) + 
  geom_ribbon(aes(ymin=y2-1, ymax=y2+1), alpha=.3, fill='pink', col='grey') +  ## alpha allows us to see the background grid
  geom_point(aes(y=y, col=k1)) +                                               ## Colorful points
  geom_line(aes(y=y2), col='black') +                                          ## Center line
  ggtitle("Title\nLine 2 Of Title") +                                          ## Add a title
  theme(plot.title = element_text(lineheight=.7, face="bold",                  ## title line space, face, color, size, & angle
                                  colour="#aa0000", size=20, angle=0)) +
  xlab('x-title\nLine 2 Of x-title') +                                         ## Add an x-axis title
  theme(axis.title.x = element_text(lineheight=0.7, face="plain",              ## x-axis line space, face, color, size, & angle
                                    colour="#ffa000", size=18, angle=0),
        axis.text.x  = element_text(angle=0, color='red',                      ## x-axis tick label angle, color, centering, & size
                                    vjust=0.5, size=10)) +
  ylab('y-title\nLine 2 Of y-title') +                                         ## Add an y-axis title
  theme(axis.title.y = element_text(lineheight=0.7, face="italic",             ## y-axis line space, face, color, size, & angle
                                    colour="#ff00a0", size=18, angle=90),
        axis.text.y  = element_text(angle=0, color='brown',                    ## y-axis tick label angle, color, centering, & size
                                    vjust=0.5, size=10)) +
  #theme(legend.title=element_blank())                                         ## This is how you nix a legend title entirely
  guides(col=guide_legend(title="HELLO")) +                                    ## Title for legend
  theme(legend.title = element_text(colour="blue",                             ## legend title size, color, face
                                    size=10, face="bold.italic"),
        legend.text  = element_text(colour="purple",                           ## legend elements size, color, face
                                    size=8, face="plain"),
        legend.position="right")                                               ## Put the legend on the right side of the graph

################################################################################################################################################################
# Simple x-y graphs with facets.  Points & Lines.  X data numeric or date.  Fancy strip options.

ggplot(data=someData, aes(x=xd, y=y, col=k1)) +              ## col=k1 has nothing to do with facets -- you could leave it off.
  geom_point() +
  #geom_line() +                                             ## lines work just as well -- could do both if you wanted!
  facet_wrap(~k1, ncol=2) +                                  ## Gets facets in a grid
  #facet_grid(.~k1)  +                                       ## Gets facets arranged left to right
  #facet_grid(k1~.)  +                                       ## Gets facets arranged top to bottom
  labs(title='Main Title',                                   ## Title at the top
       x='X Title',                                          ## On the bottom -- centered
       y='Y Title') +                                        ## On the left -- centered
  theme(strip.text.x = element_text(size=8, face="plain",    ## Strips along y-axis text
                                    colour='red', angle=0),
        strip.text.y = element_text(size=8, face="bold",     ## Strips along y-axis text
                                    colour='blue', angle=0),
        strip.background = element_rect(colour="blue",       ## Strip border and background
                                        fill="pink"))


################################################################################################################################################################
# Simple x-y graphs with facets and trends (lowess & linear regression)

ggplot(data=someData, aes(x=xd, y=y)) +
  geom_point(col='pink') +                             ## Leave off to just see the smoothed line and interval (pink shows off the line and inteval)
  #geom_smooth(method="loess") +                       ## DEFAULT. lowess.  Show confidence interval. 
  #geom_smooth(method="lm") +                          ## Linear model instead.  
  geom_smooth(method="lm", level=0.9999) +             ## Linear model instead.  Explicitly set confidence level for interval
  #geom_smooth(method="lm", se=FALSE) +                ## Don't show the confidence interval
  #geom_smooth(method="gam", formula = y~s(x)) +       ## Generalised additive model.  Needs library(mgcv)
  #geom_smooth(method="rlm") +                         ## Robust linear regression. Needs library(MASS)
  facet_wrap(~k1, ncol=2)

################################################################################################################################################################
# Simple x-y graphs with linear regression lines

ggplot(data=someData, aes(x=xd, y=y, col=k1)) +
  geom_point() +                                      ## Leave this off and you will just see the smoothed line and interval
  geom_smooth(method="lm", se=FALSE)                  ## Don't show the confidence interval
  #geom_smooth(method="loess", span=.2, se=FALSE)     ## lowess.  No confidence interval. Note the span=.2 -- lowess(x, f=.2)

################################################################################################################################################################
# Linear regression used for future prediction

expandedRange <- c(min(someData$xi),                     ## Range from min to max+1/2 the range.
                   max(someData$xi) +
                       diff(range(someData$xi))/2) 
ggplot(data=someData, aes(x=xi, y=y)) +
  scale_x_continuous(limits = expandedRange) +         ## Extend the x-axis on the right. NOTE: don't use coord_cartesian here!!
  #geom_line() +                                       ## Add this if you want to connect the dots. ;)
  geom_point() +                                       ## You can also use points!
  geom_smooth(method="lm", fullrange=TRUE, level=0.99) ## Linear model with a .99 confidence interval

################################################################################################################################################################
# Area graph
        
ggplot(data=subset(someData, k1=='K1F_1'), aes(x=xi, y=abs(y))) +
  geom_area(col = "black", fill='red')

################################################################################################################################################################
# Stacked area graph -- Yuck!

ggplot(data=someData, aes(x=xi, y=abs(y), fill=k1)) +
  geom_area(stat='identity', position="stack")

################################################################################################################################################################
# Basic Box-n-Wisker

ggplot(data=someData, aes(x=k1, y=y)) + 
  geom_boxplot(col='red', fill='pink')

################################################################################################################################################################
# Colorful Box-n-Wisker -- useful to match colors with other graphs in a report

ggplot(data=someData, aes(x=k1, y=y, fill=k1))+
  geom_boxplot(show.legend=FALSE)                  ## Suppress the legend
   
################################################################################################################################################################
# Exerting a little control over the axis and legend lables 

ggplot(data=someData, aes(x=k1, y=y, fill=k1)) +  geom_boxplot(col='black', alpha=.4) +
  #scale_x_discrete(labels=c("x1", "x2", "x3", "x4"))      ## x-axis: Set the tick marks
  #scale_x_discrete(breaks=NULL)                           ## x-axis: Nix the tick marks and grid lines
  scale_fill_discrete(name="Title\nSecond Line Of Title",  ## Set title of legend
                      labels=c("x1", "x2", "x3", "x4"))    ## Set labels of legend

################################################################################################################################################################
# A standard violin plot -- note the white borders.  I think this helps them stand out in the standard color scheme.

ggplot(data=someData, aes(x=k1, y=y, fill=k1)) + 
  geom_violin(col='white', show.legend=FALSE)            ## Suppress the legend and put a white border on the violins

################################################################################################################################################################
# Combine a violin and box-n-wisker plot -- This is how I normally do them.

ggplot(data=someData, aes(x=k1, y=y, fill=k1)) + 
  geom_boxplot(col='black', alpha=.4) +
  geom_violin(alpha=.25, col=NA) +
  theme(legend.position="none")

################################################################################################################################################################
# Grab a graph that is already created, save it as a PDF file, and load that PDF up with a PDF viewer.  See the margin annotation example below for another
# option.

#ggsave("x.pdf", width=10, height=8, dpi=300, units='in'); system('xpdf x.pdf &')

################################################################################################################################################################
# How to create a graphic object, and then apply a theme to it.

aGrob <- ggplot(data=someData, aes(x=xd, y=y)) + facet_wrap(~k1, ncol=2) + geom_line() + geom_smooth(method="lm")

# Here is the default -- remember, a grob is plotted when the REPL 'prints' it so the following generates a plot
aGrob + labs(title='Theme: DEFAULT (same as theme_grey)') 

# I sometimes use this for printed material
aGrob + labs(title='Theme: light  ') + theme_light()

# I have used this as a starting place for new themes
aGrob + labs(title='Theme: minimal') + theme_minimal()          

################################################################################################################################################################
# Annotate within the plot region

qplot(data=someData,x=xi, y=y, color=k1) +
  geom_abline(intercept=0, slope=1, col='blue', size=3) +                                                    ## A 'data' line
  geom_hline(yintercept=-.5, col='red') +                                                                    ## A horizontal line
  geom_vline(xintercept=.55, col='red') +                                                                    ## A vertical line
  annotate("text", x = .25, y = .75, label = "HI", col='red', size=14) +                                     ## Text
  annotate("rect", xmin=0, xmax=.5, ymin=0, ymax=1.5, alpha=.1, fill='red', col='black') +                   ## Rectangles
  annotate("segment", x = 0.0, xend = 0.25, y = -2.0, yend = -1.0) +                                         ## Segments
  annotate("segment", x = 1.0, xend = 0.75, y = -2.0, yend = -1.0, arrow = arrow(length = unit(0.5,"cm")))   ## Arrows!

################################################################################################################################################################
# Annotate the bottom margin on the right

aGrob <- arrangeGrob(qplot(data=someData,x=xi, y=y, color=k1),
                     bottom=textGrob("HELLO\ndog\nmy name is\nfoo", x = unit(1, 'npc')-unit(4, 'mm'), just = c("right", "bottom"),
                                    gp = gpar(fontface = "italic", fontsize = 10, col='blue')))
grid.newpage()
grid.draw(aGrob)

################################################################################################################################################################
# Annotate the bottom margin on the left
 
aGrob <- arrangeGrob(qplot(data=someData,x=xi, y=y, color=k1),
                     bottom=textGrob("HELLO\ndog\nmy name is\nfoo", x = unit(0, 'npc')+unit(4, 'mm'), just = c("left", "bottom"),
                                     gp = gpar(fontface = "italic", fontsize = 10, col='blue')))
grid.newpage()
grid.draw(aGrob)

################################################################################################################################################################
# Annotate the right (or left) margin with a significant quantity of text. My normal use case is a paragraph explaining some aspect of how the graph was created
# -- things like underlying assumptions required to understand the graph.  paragraph.  Put the text on the left instead of the right by replacing 'legend' with
# 'left'

ourText <- "HELLO\nWhere we discuss the\nmany tiny little\nthings we faked\nin the production of our data."
aGrob <- arrangeGrob(qplot(data=someData,x=xi, y=y, color=k1),
                     right=textGrob(ourText, x = unit(2, 'mm'), y = unit(1, 'npc')-unit(2, 'mm'), just = c("left", "top"),
                              gp = gpar(fontface = "italic", fontsize = 10, col='red')))
grid.newpage()
grid.draw(aGrob)

# You could save and display a PDF like so:

#ggsave("x.pdf", width=10, height=8, dpi=300, units='in', plot=aGrob); system('xpdf x.pdf &')

################################################################################################################################################################
# 2D data: An image

ggplot(data=twoDimData, aes(x=x, y=y, fill=z)) +
  geom_tile() 
  #geom_raster()                                ## When length of x and y are the same, we can use the faster geom_raster()

################################################################################################################################################################
# 2D data: Image with a dot in each cell scaled to 'abs(z)'

ggplot(data=twoDimData, aes(x=x, y=y)) +
  geom_tile(aes(fill=z)) +
  geom_point(aes(size=abs(z)), col='red')

################################################################################################################################################################
# 2D data: Image with white text in each cell

ggplot(data=twoDimData, aes(x=x, y=y)) +
  geom_tile(aes(fill=z)) +
  geom_text(aes(label=p), col='white', size=4)

################################################################################################################################################################
# 2D data: An image with text in each cell with a color set by the z value

ggplot(data=twoDimData, aes(x=x, y=y)) +
  geom_tile(aes(fill=z)) +
  geom_text(aes(label=p), col=c('red', 'black', 'green')[sign(twoDimData$z)+2], size=4)

################################################################################################################################################################
# 2D data: An Image with contours in white

ggplot(data=volcanoDF, aes(x=x, y=y)) +
  geom_raster(aes(fill=values), interpolate=TRUE) +  # tile has no "interpolate" option.
  geom_contour(aes(z=values), col='white', size=1)

################################################################################################################################################################
# 2D data: Using stat_contour with polygons to fill in contours (not geom_tile for an image).  White contour lines.

ggplot(data=volcanoDF, aes(x=x, y=y, z=values)) +
  stat_contour(geom="polygon", aes(fill=..level..))  + 
  stat_contour(col='white', size=1)

################################################################################################################################################################
# 2D data: Just contour lines without an image.  The contour line color determined by contour level

ggplot(data=volcanoDF, aes(x=x, y=y, z=values)) +
  geom_contour(aes(col=..level..), size=2)              ##  Fatten up the line so the color shows up

################################################################################################################################################################
# 2D histogram with rectangular or hexagon bins

ggplot(data=twoDimNormal, aes(x=x,y=y)) + 
  #geom_rug() +                           ## Add a rug (dot-plot) to each axis for lower density plots
  stat_bin2d(aes(fill=..count..))         ## Use this for rectangular bins!
  #stat_binhex(aes(fill=..count..))       ## Use this for hexagon bins -- fill artifacts on raster devices.  Use PDF.

################################################################################################################################################################
# 2D Density Estimation with semi-transparent data points visually indicating higher density

ggplot(data=twoDimNormal, aes(x=x,y=y)) +
  geom_point(alpha=.2, col='red') +           ## Alpha lets points visually producing darker shades in high density regions
  #geom_rug() +                               ## Add a rug (dot-plot) to each axis for lower density plots
  geom_density2d(col='black', size=1)         ## Put contour lines after points to make sure we can see them.

################################################################################################################################################################
# 2D Density Estimation shown as a classic, filled contour graph

ggplot(data=twoDimNormal, aes(x=x,y=y)) +
  geom_point(alpha=.5, col='black') +         ## Show outlier with dots -- must be first so contour covers most of them up
  #geom_rug() +                               ## Add a rug (dot-plot) to each axis for lower density plots
  stat_density2d(aes(fill = ..level..),       ## Fill in the contour graph -- covering up non-outlier points.
                 geom="polygon", col='white')

################################################################################################################################################################
# Scatter plot with marginal histograms

histTop <- ggplot(twoDimNormal) +                                     ## Create histogram that goes at the top
  geom_histogram(aes(x=x),
                 col='white',
                 fill='red',
                 binwidth=diff(range(twoDimNormal$x))/50) +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(margin=margin(0,0,0,0,"pt")),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0,0,0,0),"lines"),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(0,"null")) +
  scale_x_continuous(limits=range(twoDimNormal$x))

histRight <- ggplot(twoDimNormal) +                                   ## Create histogram that goes at the right
  geom_histogram(aes(x=y),
                 col='white',
                 fill='red',
                 binwidth=diff(range(twoDimNormal$y))/50) +
  coord_flip() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(margin=margin(0,0,0,0,"pt")),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,0),"lines"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(0,"null")) +
  scale_x_continuous(limits=range(twoDimNormal$y))

maxCount = max(c(max(ggplot_build(histTop)$data[[1]]$count),          ## Set the scale of both histograms to the same value
                 max(ggplot_build(histRight)$data[[1]]$count)))
histTop   <- histTop + scale_y_continuous(limits=c(0,maxCount+1))
histRight <- histRight + scale_y_continuous(limits=c(0,maxCount+1))

scatter <- ggplot(twoDimNormal)+                                      ## Create scatter plot in the center.  Use transparency to indicate density.
  geom_point(aes(x=x,y=y), col=rgb(1,0,0,.05)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,0),"lines"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(0,"null")) +
  scale_x_continuous(limits=range(twoDimNormal$x)) + 
  scale_y_continuous(limits=range(twoDimNormal$y))

aGrob <- arrangeGrob(histTop,                                         ## Put it all together
                     grob(),
                     scatter,
                     histRight,
                     ncol=2,
                     nrow=2,
                     widths=c(3, 1),
                     heights=c(1, 3))
grid.newpage()
grid.draw(aGrob)

################################################################################################################################################################
# Barcharts with pre-computed data -- i.e. you have 'labels' and values

ggplot(data=someLessData, aes(x=quarter, y=value)) +
  geom_bar(stat='identity', col='black', fill='red') +       ## Draws red bars with black borders
  geom_text(aes(label=value), vjust='top', nudge_y=-0.25) +  ## Adds the numerical label to each bar
  theme(panel.grid.minor.x=element_blank(),                  ## Get rid of the vertical grid lines
        panel.grid.major.x=element_blank())

################################################################################################################################################################
# Barcharts with pre-computed data with color filled bars -- i.e. you have 'labels' and values Why? To keep the bar colors consistent with other color coded
# plots in the same report (like the xy plots shown above).

ggplot(data=someLessData, aes(x=quarter, y=value, fill=quarter)) +
  geom_bar(stat='identity', col='black', show.legend=FALSE) +        ## col adds the black lines separating the bars
  geom_label(aes(label=value), fill='white', vjust='center') +       ## Adds the numerical label to each bar.  geom_label is geom_text with a box
  theme(panel.grid.minor.x=element_blank(),                          ## Get rid of the vertical grid lines
        panel.grid.major.x=element_blank())

################################################################################################################################################################
# The every-day stacked barchart

ggplot(data=someMoreData, aes(x=quarter, y=value, fill=factors)) +
  guides(fill=guide_legend(override.aes=list(colour=NA))) +       ## Get rid of the slash across the legend color boxes
  geom_bar(stat='identity', col='black', position="stack") +
  theme(panel.grid.minor.x=element_blank(),                       ## Get rid of the vertical grid lines
        panel.grid.major.x=element_blank())

################################################################################################################################################################
# Side by side barchart

ggplot(data=someMoreData, aes(x=quarter, y=value, fill=factors)) +
  guides(fill=guide_legend(override.aes=list(colour=NA))) +       ## Get rid of the slash across the legend color boxes
  geom_bar(stat='identity', col='black', position="dodge") +
  theme(panel.grid.minor.x=element_blank(),                       ## Get rid of the vertical grid lines
        panel.grid.major.x=element_blank())

################################################################################################################################################################
# Stacked with constant height -- normally used instead of a pie chart to compare proportions

ggplot(data=someMoreData, aes(x=quarter, y=value, fill=factors)) +
  guides(fill=guide_legend(override.aes=list(colour=NA))) +       ## Get rid of the slash across the legend color boxes
  geom_bar(stat='identity', col='black', position="fill") +
  theme(panel.grid.minor.x=element_blank(),                       ## Get rid of the vertical grid lines
        panel.grid.major.x=element_blank())

################################################################################################################################################################
# Circular bar chart -- yea, a pie.  Yuck.

ggplot(someLessData, aes(x=factor(1), y=value, fill=quarter)) +
  geom_bar(width=1, stat='identity', col='black') +         ## col puts black lines between slices
  guides(fill=guide_legend(override.aes=list(colour=NA))) + ## Get rid of the slash across the legend color boxes
  coord_polar(theta='y') +                                  ## This is how it gets round
  theme(axis.ticks=element_blank(),                         ## Get rid of axis ticks and labels
        axis.text.y=element_blank(),
        axis.text.x=element_text(colour='black')) +
  labs(title='Main Title',                                  ## You can set the title, but the x & y are ignored
       x='', y='')                                          ## Ignored!

################################################################################################################################################################
# Basic 1D histogram

ggplot(data=someData, aes(x=rn)) +
  geom_histogram(col = "black", fill = "red", binwidth=diff(range(someData$rn))/20)

################################################################################################################################################################
# Basic 1D histogram with percentage instead of count

ggplot(data=someData, aes(x=rn, y = (..count..)/sum(..count..))) +
  geom_histogram(col = "black", fill = "red", binwidth=diff(range(someData$rn))/20) + 
  scale_y_continuous(labels=percent) +
  labs(y='%')

################################################################################################################################################################
# 1D histogram with bin breaks hard coded -- inconsistent sizes in this case

ggplot(data=someData, aes(x=rn)) +
  geom_histogram(col = "black", fill = "red", breaks=c(-1.0,-0.5,-0.25,-0.1,0.0,0.1,0.25,0.5,1.0))

################################################################################################################################################################
# 1D histogram with fill color determined by bin count

ggplot(data=someData, aes(x=rn, fill=..count..)) + 
  geom_histogram(col = "black", binwidth = .1) +   ## 'col' gets us black lines separating bars
  theme(legend.position="none")                    ## legend provides no new information (y-axis shows bar height)

################################################################################################################################################################
# 1D histogram + 1D density curve

ggplot(data=someData, aes(x=rn, y=..density..)) +
  geom_histogram(col = "black", fill = "red", binwidth = .11) + ## First so we always see density line.  Leave off for just line.
  geom_density(col = "blue", size=2)                            ## Fatten up line so we can see it

################################################################################################################################################################
# Draw two graphs from the same data frame on the same set of axes

ggplot(someData) +
  geom_point(aes(x=xi,  y=y, col='dots')) + # Note use of col= in aes()
  geom_line(aes(x=xi2, y=y2, col='line')) +
  labs(colour="foo")

################################################################################################################################################################
# Draw two graphs potentially from different data frames on the same set of axes

ggplot() +
  geom_point(data=twoDimNormal, aes(x=x,  y=y, col='big'),    size=2, alpha=.25) +
  geom_point(data=someData,     aes(x=xi, y=y, col='little'), size=3, alpha=.25) +
  scale_colour_manual(name='foo', values=c('big'='red', 'little'='blue'), labels=c('little'='Da Small One', 'big'='Da Big One'))

################################################################################################################################################################
# Draw two graphs potentially from different data frames on the same set of axes with with fixed colors and no legend

ggplot() +
  geom_point(data=twoDimNormal, aes(x=x,  y=y), col='red',  size=2, alpha=.25) +
  geom_point(data=someData,     aes(x=xi, y=y), col='blue', size=3, alpha=.25)


################################################################################################################################################################
# A completely custom pallet.  For example, a color blind safe palette from http://jfly.iam.u-tokyo.ac.jp/color/
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(data=someData, aes(x=k1, y=y, fill=k1)) +
  geom_boxplot(col='black', alpha=.4) +
  scale_fill_manual(values=cbPalette)

ggplot(data=someData, aes(x=xd, y=y, col=k1)) + geom_point() +
  scale_colour_manual(values=cbPalette)

################################################################################################################################################################
# Use a colorbrewer pallet
ggplot(data=someData, aes(x=k1, y=y, fill=k1)) +
  geom_boxplot(col='black', alpha=.4) +
  scale_fill_brewer(palette="Set2")

ggplot(data=someData, aes(x=xd, y=y, col=k1)) + geom_point() +
  scale_colour_brewer(palette="Set1")
