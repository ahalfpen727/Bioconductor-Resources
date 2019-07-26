# -*- Mode:R; Coding:us-ascii-unix; fill-column:160 -*-

################################################################################################################################################################
##
# @file      longAndWide.R
# @author    Mitch Richling <https://www.mitchr.me>
# @Copyright Copyright 2015 by Mitch Richling.  All rights reserved.
# @brief     Transforming between long and wide data sets.@EOL
# @Keywords  r cran reshape2::melt reshape2::dcast utils::unstack utils::stack stats::aggregate tidyr::gather tidyr::spread
#
# We demonstrate two reshape2 commands:
#   * reshape2::melt -- transform from 'wide' to 'long'
#       * arg 1 is the data
#       * id.vars -- the variables to keep but not split apart on
#       * measure.vars -- source columns
#       * variable.name -- Name of the destination column that will identify the original
#       * value.name - column that the measurement came from
#   * reshape2::dcast --transform from 'long' to 'wide'
#       * Arg 1 is the data
#       * Arg 2 is a formula
#       * foo + bar ~ foobar + foobam
#         LHS are columns you want to keep as is (vertical)
#         RHS are columns who's values will be used to construct new columns (horizontal)
#       * value.var is the variable with the measurement you want in the table
#       * fun.aggregate is used when the table will have multiple value.var values for each cell.
#         This function is then used to aggregate the multiple values into one for the table

################################################################################################################################################################
# First, load the library.

library(reshape2)
library(tidyr)

################################################################################################################################################################
# This is the most trivial, and common, use case when converting from 'wide' to 'long'

someWideData <- read.table(header=T, text='
  subject sex breakfast dinner lunch
1     Dad   M     101.9  102.1  99.2
2     Mom   F      98.6   98.9  95.4
3     Son   M     120.3  120.4 110.1
')

reshape2::melt(someWideData, id.vars=c('subject', 'sex'), measure.vars=c('breakfast', 'dinner', 'lunch'), variable.name='meal', value.name='blgl')

tidyr::gather(someWideData, meal, blgl, breakfast:lunch)

################################################################################################################################################################
# If you don't have an 'x variable', then you can use 'aggregate' to make a wide table horizontal:

head(InsectSprays)

stats::aggregate(InsectSprays$count, list(spray=InsectSprays$spray), identity)

################################################################################################################################################################
# Now lets look at the most common wide to long use case:

someLongData <- read.table(header=T, text='
 subject sex        meal blgl
     Dad   M   breakfast  101.9
     Dad   M       lunch   99.2
     Dad   M      dinner  102.1
     Mom   F   breakfast   98.6
     Mom   F       lunch   95.4
     Mom   F      dinner   98.9
     Son   M   breakfast  120.3
     Son   M       lunch  110.1
     Son   M      dinner  120.4
 ')

reshape2::dcast(someLongData, sex ~ meal, value.var='blgl', fun.aggregate=mean)

################################################################################################################################################################
# Another common 'wide' to 'long' case is when we have more than one variable that should be on the vertical axis of the table.  For example, suppose we want
# 'subject' & 'sex' on the vertical, 'meal' on the horizontal, and the unique measurement for that combination in the table.  Then we might do this:
  
reshape2::dcast(someLongData, subject + sex ~ meal, value.var='blgl')

tidyr::spread(someLongData, meal, blgl)

################################################################################################################################################################
# Lets do another example with a different dataset -- mtcars. Suppose we want a table with 'cyl' & 'gear' on the vertical, 'carb' on the horizontal, and mean
# 'mpg' in the table.  We can do that like so:

head(mtcars, 20)
reshape2::dcast(mtcars, cyl + gear ~ carb, value.var='mpg', fun.aggregate=mean)

################################################################################################################################################################
# Lets do one more example with yet another data set -- tips. We want 'time' on vertical, and both 'sex' & 'smoker' on the horizontal, and average tip in the
# table:
  
head(tips)
reshape2::dcast(tips, time ~ sex + smoker, value.var='tip', fun.aggregate=mean)

################################################################################################################################################################
# Sometimes you have a group of vectors each containing measurements for different groups -- or a data frame with columns that are each measurements for
# different groups, and you simply want to stack them into a two column, 'long' data frame.

utils::stack(list(group1=c(0,1,2,3), group2=c(4,5,6,7)))

utils::stack(data.frame(group1=c(0,1,2,3), group2=c(4,5,6,7)))

################################################################################################################################################################
# Alternately, if you have a data.frame with two columns -- one for the group and one for the measurements, then you can get a wide

a <- data.frame(vals=c(0,1,2,3,4,5,6), group=c('group1', 'group1', 'group1', 'group2', 'group2', 'group2', 'group2'))
a

utils::unstack(a, vals ~ group)

################################################################################################################################################################
# In the case when the number of measurements in each group is equal, unstack will return a data.frame

a <- data.frame(vals=c(0,1,2,3,4,5), group=c('group1', 'group1', 'group1', 'group2', 'group2', 'group2'))
a

utils::unstack(a, vals ~ group)

