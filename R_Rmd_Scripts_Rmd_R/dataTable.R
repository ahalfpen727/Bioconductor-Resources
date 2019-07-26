# -*- Mode:R; Coding:utf-8; fill-column:160 -*-

################################################################################################################################################################
# @file      dataTable.R
# @author    Mitch Richling <https://www.mitchr.me>
# @Copyright Copyright 2015 by Mitch Richling.  All rights reserved.
# @brief     Some data.table basics.@EOL
# @Keywords  r package cran data.table dataTable
#
# Some of the basics of data.table.  Note that dplyr works with data.table.  So if you have been using dplyr with data.frames, you can get a significant boost
# in performance by simply dropping in data.table -- no code changes required!
#

################################################################################################################################################################
# We will be working with the following data

mtcars

################################################################################################################################################################
# One can directly coerce a data.frame into a data.table; however, row names will be lost.  We can do something like this to keep them:

mtcarsDT <- data.table(mtcars, names=rownames(mtcars))
mtcarsDT

################################################################################################################################################################
# Grab row #2.  

mtcarsDT[2,]

################################################################################################################################################################
# Note that when we only have one index, we don't need the comma!!

mtcarsDT[2]

################################################################################################################################################################
# Get rows 2 through 5

mtcarsDT[2:5]

################################################################################################################################################################
# Get rows where cyl is equal to 6

mtcarsDT[cyl==6]

################################################################################################################################################################
# Grab column "cyl" the data.frame way.

mtcarsDT$cyl

################################################################################################################################################################
# The second index of a data.table is an expression on the rows selected by the first index.  In this example, the expression is "cyl".  As this will evaluate
# to the cyl column, that is what is returned.

mtcarsDT[,cyl]

################################################################################################################################################################
# When the second index is a list, then a new data.table will be returned instead of a simple vector. For example we can extract a new data.table with three
# columns (two of which we rename)

mtcarsDT[,list(cylinder=cyl, weight=wt, gear)]

################################################################################################################################################################
# The result of the computation in the second index can be stuffed back into the data.table.

mtcarsDT[,cylTimesGear:=cyl*gear]

################################################################################################################################################################
# That second expression can do more than transform columns into new columns; it aggregate things.  This is roughly equivalent to: with(mtcars, sum(wt))

mtcarsDT[,sum(wt)]

################################################################################################################################################################
# As before, if the second expression is a list we get a new data.table.  If the second expression computes aggregates, then we get more than one
# aggregation. This is roughly equivalent to: with(mtcars, list(daSum=sum(wt), daSd=sd(wt)))

mtcarsDT[,list(daSum=sum(wt), daSd=sd(wt))]

################################################################################################################################################################
# We can aggregate by groups too.  This is roughly equivalent to: tapply(mtcars$wt, mtcars$cyl, sum)

mtcarsDT[,sum(wt),cyl]

################################################################################################################################################################
# The groups can contain more than one factor if we use a list.  We can also name the aggregate column if we put that in a list.

mtcarsDT[,list(sumWt=sum(wt)),list(cyl, gear)]

################################################################################################################################################################
# We can aggregate by groups and put it back into the data.table!!

mtcarsDT[,sumWtByCyl := sum(wt),cyl]

################################################################################################################################################################
# If you have a "key column" set, you can index with strings -- note that this will change the order of the table to store on the key.

setkey(mtcarsDT,names)
mtcarsDT["Valiant"]

################################################################################################################################################################
# The values in the "key column" need not be unique:

mtcarsDT$cylNames <- paste("cyl", mtcarsDT$cyl, sep='=')
setkey(mtcarsDT,cylNames)
mtcarsDT["cyl=6"]

################################################################################################################################################################
# You can get just the first result from a key column lockup like so (get the last with "last")

mtcarsDT["cyl=6",mult="first"]

################################################################################################################################################################
# With things sorted on cyl, things like cumulative sums become interesting -- especially when added to the table.

mtcarsDT[,cSumWtByCyl := cumsum(wt),by=cyl]

################################################################################################################################################################
# If the key column represents experiment tags, then we can remove duplicate experiments like so

unique(mtcarsDT)

################################################################################################################################################################
# Merge is much like with data.frames.  One nice feature is that key columns will be used for a merge automatically if they are set.

crbn <- data.table(carb=c(1,2,3), carbs=c("one", "two", "three"))
setkey(mtcarsDT, carb)
setkey(crbn, carb)
merge(mtcarsDT, crbn, all.x=TRUE)

################################################################################################################################################################
# Read a CSV into a data.table (much faster than read.csv & read.table).  With CSV, things usually just work.

df1 <- fread("dataTable_f1.csv")
df1

################################################################################################################################################################
# Same as previous, but from a URL!!
df1u <- fread('https://www.mitchr.me/SS/exampleR/rcode/dataTable_f1.csv')
df1u

################################################################################################################################################################
# Read output from a command (in this case compressing with gunzip and filtering with awk)

df2 <- fread("gunzip < dataTable_f2.csv.gz | awk -F, 'NR==1 || $3<100 { print $0 }'")
df2

################################################################################################################################################################
# Read a colon (:) separated file with extra whitespace and no column names

df3 <- fread('dataTable_f3.txt', sep=':', header=FALSE, strip.white=TRUE, col.names=c('name', 'age', 'weight'))
df3
