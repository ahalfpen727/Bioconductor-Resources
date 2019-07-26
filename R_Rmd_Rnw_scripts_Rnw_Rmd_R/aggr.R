# -*- Mode:R; Coding:us-ascii-unix; fill-column:160 -*-

################################################################################################################################################################
##
# @file      aggr.R
# @author    Mitch Richling <https://www.mitchr.me>
# @Copyright Copyright 2015 by Mitch Richling.  All rights reserved.
# @brief     Simple aggregation using base R.@EOL
# @Keywords  base r aggregation tapply aggregate dcast by lapply
#
# Also check out the data.table and dplyr examples.  Note that the 'dcast' function in reshape2.R can do some aggregation as well.
#            

################################################################################################################################################################
# We will use the mtcars data set for all the demonstrations

mtcars

################################################################################################################################################################
# Apply a function to ONE column ('disp' as a vector) of the data frame broken up by ONE factor ('cyl')

tapply(mtcars$disp, mtcars$cyl, mean)

################################################################################################################################################################
# Apply a function to ONE column ('disp' as a vector) of the data frame broken up by TWO factors ('cyl' & 'gear').  In the resulting table the rows (4, 6 & 8)
# are 'cyl' while columns are (3, 4 & 5) are 'gear'.

tapply(mtcars$disp,list(mtcars$cyl, mtcars$gear), mean)

################################################################################################################################################################
# With a single column in the list for the second argument, aggregate works much like tapply but returns a data.frame.  Note that this preserves the type of the
# 'cyl' column -- it is retained as the 'names' part of the returned vector tapply.

aggregate(mtcars$disp,list(cyl=mtcars$cyl), mean)

################################################################################################################################################################
# Unlike tapply, aggregate can work with multiple data columns -- aggregated independently, and returns a data.frame.

aggregate(mtcars$disp,list(cyl=mtcars$cyl, gear=mtcars$gear), mean)

################################################################################################################################################################
# aggregate has a handy formula based interface too

aggregate(disp ~ cyl, mtcars, mean)

################################################################################################################################################################
# aggregate has a handy formula based interface too

aggregate(disp ~ cyl + gear, mtcars, mean)

################################################################################################################################################################
## dcast (from reshape2) can do some aggregation too. Note it sorts the resulting data.frame.

dcast(mtcars, cyl + gear ~ . , fun.aggregate=mean, value.var='disp')

################################################################################################################################################################
# by is tapply, but for data.frames instead of vectors -- i.e. the function gets a data.frame for each factor level.

byo <- by(mtcars, mtcars$cyl, function(x) { mean(x$disp); } )
byo

names(byo)

as.vector(byo)

################################################################################################################################################################
# by can return complex results

by(mtcars, mtcars$cyl, function(x) { c(mean(x$disp), sd(x$disp)); } )

################################################################################################################################################################
# split actually splits up a data.frame into a list of data.frames broken out by a factor

lapply(split(mtcars, mtcars$cyl), function(x) { length(x$cyl) })
