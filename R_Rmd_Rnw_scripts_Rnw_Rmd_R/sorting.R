# -*- Mode:R; Coding:us-ascii-unix; fill-column:160 -*-

################################################################################################################################################################
##
# @file      sorting.R
# @author    Mitch Richling <https://www.mitchr.me>
# @Copyright Copyright 2015 by Mitch Richling.  All rights reserved.
# @brief     Sorting things in R -- sometimes harder than it should be.@EOL
# @Keywords  sort base r
#
# See the 'arrange' and 'add_rownames' functions in the dplyr:: package
#            

################################################################################################################################################################
# sortong a vector is easy -- note that the right thing happens with names
aVector        <- as.integer(rnorm(25, mean=20, sd=10))
names(aVector) <- as.integer(rnorm(25, mean=20, sd=10))
aVector
sort(aVector)

################################################################################################################################################################
# sort by cyl
mtcars[order(mtcars$cyl),]

################################################################################################################################################################
# sort by cyl and carb
mtcars[order(mtcars$cyl, mtcars$carb),]

################################################################################################################################################################
#sort by cyl (ascending) and carb (descending)
mtcars[order(mtcars$cyl, -mtcars$carb),]


