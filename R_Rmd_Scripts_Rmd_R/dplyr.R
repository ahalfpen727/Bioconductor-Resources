# -*- Mode:R; Coding:us-ascii-unix; fill-column:160 -*-

################################################################################################################################################################
##
# @file      dplyr.R
# @author    Mitch Richling <https://www.mitchr.me>
# @Copyright Copyright 2015 by Mitch Richling.  All rights reserved.
# @brief     Demo of some of my most used dplyr features.@EOL
# @Keywords  package cran dplyr R
#
# Notes on the design of dplyr:
#
#  1) Operations in dplyr don't modify the arguments passed to them, but rather return new results -- i.e. they don't have side effects.  This means that you
#     can use functional-style programming with dplyr.
#
#  2) dplyr tends to be faster than the equivalent base R functionality
#
#  3) dplyr works on data.tables, data.frames, and SQL connections.  This is very cool because it provides a uniform interface to all of those objects.  For
#     example, you can write your code for data.frames and it will still work later if your data gets bigger and you switch from data.frames to data.tables.

################################################################################################################################################################
# First we load up the library

library(dplyr)

################################################################################################################################################################
# The discussion below is example focused.  We will be using the same dataset for all of the examples that follow:

mtcars

################################################################################################################################################################
# First we need to have some data to work with.

daObjects <- read.table(header=T, text='
     hairs SEX  name
  brunette   M  rock
    blonde   M  paper
  brunette   F  scissors
    blonde   F  hammer
       red   M  gun
       red   F  knife
', as.is=TRUE)

daPeople <- read.table(header=T, text='
      hair sex  name
  brunette   M  bill
    blonde   M  ted
    blonde   M  fred
  brunette   F  wilma
    blonde   F  betzi
      pink   F  pinkey
', as.is=TRUE)

################################################################################################################################################################
# Most dplyr operations zap row names on data.frames.   -- so add a real column with the row.names

mtcarsR <- dplyr::add_rownames(mtcars)                         
mtcarsR

################################################################################################################################################################
# Grab rows by Boolean expression. You can actually replace that & with a comma, but I'm usually explicit about this sort of thing

dplyr::filter(mtcarsR, cyl>4 & carb==4)                        

################################################################################################################################################################
# Grab rows by index

dplyr::slice(mtcarsR, 15:20)

################################################################################################################################################################
# Sort things.  Could have used '-gear' here as 'gear' is numeric, but 'desc' works on strings as well.

dplyr::arrange(mtcarsR, cyl, desc(gear)) 

################################################################################################################################################################
# Select columns.  Love how you can use column names with the range operator.  Numbers work too, but that is boring.

dplyr::select(mtcarsR, cyl:drat)         

################################################################################################################################################################
# Pull out a column and rename it

dplyr::select(mtcarsR, displacement=disp)

################################################################################################################################################################
# You can separate selection rules with commas (works like and)

dplyr::select(mtcarsR, displacement=disp, cyl)

################################################################################################################################################################
# If you wanted all the cols, but just wanted to rename a few you can do this.

dplyr::rename(mtcarsR, displacement=disp, weight=wt)

################################################################################################################################################################
# Pull out a row for each unique value of cyl & carb.  More useful if dplyr::select is used to limit input to the column(s)

dplyr::distinct(mtcarsR, cyl, carb)           

################################################################################################################################################################
# Add new column(s)

dplyr::mutate(mtcarsR, mpc = mpg/cyl, impc=1/mpc) 

################################################################################################################################################################
# New table with JUST the new computed columns. No back references in formulas to other new columns.

dplyr::transmute(mtcarsR, mpc = mpg/cyl, mpd=mpg/disp)

################################################################################################################################################################
# Compute various aggregate functions on the columns

dplyr::summarise(mtcarsR, avgmpg=mean(mpg), sdmpc=sd(mpg/cyl))

################################################################################################################################################################
# Combining the previous functions (especially 'summarize') with a 'group_by' object is where the true utility of dplyr shines.

by_cyl <- dplyr::group_by(mtcarsR, cyl)
dplyr::summarise(by_cyl, cnt=n(), meanmpg=mean(mpg, na.rm = TRUE))

################################################################################################################################################################
# mutate and group_by are a magical combination!!!

dplyr::mutate(dplyr::group_by(mtcarsR, cyl), meanMPGbyCYL=mean(mpg))

################################################################################################################################################################
# Chaining is a technique using the %>% operator that can be used to "chain" together dplyr calls into sequential steps.  The idea is very much like how pipes
# are used on the UNIX command line.

mtcarsR %>%
  dplyr::group_by(cyl, carb) %>%
  dplyr::select(mpg, hp, cyl, carb) %>%
  dplyr::summarise(meanMPG = mean(mpg, na.rm = TRUE),
                   meanCYL = mean(hp, na.rm = TRUE))  %>%
  filter(meanMPG > 20)

################################################################################################################################################################
# Some people don't like chaining, and you are not forced to use it if you don't wish.  The previous example could have been done like this:

tmp <- dplyr::group_by(mtcarsR, cyl, carb)
tmp <- dplyr::select(tmp, mpg, hp, cyl, carb)
tmp <- dplyr::summarise(tmp,
                        meanMPG = mean(mpg, na.rm = TRUE),
                        meanCYL = mean(hp, na.rm = TRUE))
tmp <- dplyr::filter(tmp, meanMPG > 20)
tmp

################################################################################################################################################################
# dplyr has nice join (merge in R) capability too

dplyr::full_join( daPeople, daObjects, by=c('hair'='hairs', 'sex'='SEX'))

dplyr::left_join( daPeople, daObjects, by=c('hair'='hairs', 'sex'='SEX'))

dplyr::right_join(daPeople, daObjects, by=c('hair'='hairs', 'sex'='SEX'))

dplyr::inner_join(daPeople, daObjects, by=c('hair'='hairs', 'sex'='SEX'))

dplyr::semi_join( daPeople, daObjects, by=c('hair'='hairs', 'sex'='SEX'))

dplyr::anti_join( daPeople, daObjects, by=c('hair'='hairs', 'sex'='SEX'))
