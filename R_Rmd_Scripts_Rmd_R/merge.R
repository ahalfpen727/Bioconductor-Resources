# -*- Mode:R; Coding:us-ascii-unix; fill-column:160 -*-

################################################################################################################################################################
##
# @file      merge.R
# @author    Mitch Richling <https://www.mitchr.me>
# @Copyright Copyright 2015 by Mitch Richling.  All rights reserved.
# @brief     Demo of The R merge (join) command.@EOL
# @Keywords  base r merge
#

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
')

daPeople <- read.table(header=T, text='
      hair sex  name
  brunette   M  bill
    blonde   M  ted
    blonde   M  fred
  brunette   F  wilma
    blonde   F  betzi
      pink   F  pinkey
')

################################################################################################################################################################
# The most common case -- we expect every row in both tables to have a precise match in the output

base::merge(x=daPeople, y=daObjects, by.x=c('hair', 'sex'), by.y=c('hairs', 'SEX'))

################################################################################################################################################################
# It is quite common for the x-table to be the master data set, and a y-table be a set of "look-up" values related to some quality of each x-table row.  In this
# case, the y-table is frequently missing some values for the master data in the x-table, but it is important that every data row appear in the output -- even
# if the 'look-up' value is missing.

base::merge(x=daPeople, y=daObjects, by.x=c('hair', 'sex'), by.y=c('hairs', 'SEX'), all.x=TRUE)

################################################################################################################################################################
# Really, this is just the reverse of the previous case.

base::merge(x=daPeople, y=daObjects, by.x=c('hair', 'sex'), by.y=c('hairs', 'SEX'), all.y=TRUE)

################################################################################################################################################################
# Lastly we do a 'full outer join' -- in this case we want a row in the output to appear if it appears in the x-table or y-table.

base::merge(x=daPeople, y=daObjects, by.x=c('hair', 'sex'), by.y=c('hairs', 'SEX'), all=TRUE)
