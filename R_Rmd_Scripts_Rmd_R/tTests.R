# -*- Mode:R; Coding:us-ascii-unix; fill-column:160 -*-

################################################################################################################################################################
##
# @file      tTests.R
# @author    Mitch Richling <https://www.mitchr.me>
# @Copyright Copyright 2015 by Mitch Richling.  All rights reserved.
# @brief     t-test.@EOL
# @Keywords  t test ttest t-test base r
#

################################################################################################################################################################
popsz  <- 1000
mean12 <- 1
mean3  <- -1
group1 <- rnorm(popsz, mean=mean12)
group2 <- rnorm(popsz, mean=mean12)
group3 <- rnorm(popsz, mean=mean3)
allDat <- stack(list(group1=group1, group2=group2, group3=group3))
names(allDat) <- c('v', 'group')

################################################################################################################################################################
ggplot(data=allDat, aes(x=v, col=group)) + geom_density()

################################################################################################################################################################
# Welch Two Sample t-test -- i.e. when you don't know the variance of the two populations is equal

t.test(group1, group2)

################################################################################################################################################################
# Two Sample t-test -- i.e. when you DO know the variance of the two populations is equal

t.test(group1, group2, var.equal=TRUE)

################################################################################################################################################################
# Do a paired T-test -- i.e. when the measurements in each group are related pairwise.  For example, the data could be temperature measurements taken with two
# thermometers each hour.

t.test(group1, group2, paired=TRUE)

################################################################################################################################################################
# One Sample t-test for equality -- i.e. you want to know if the sample mean is equal to a hypothesized population mean

t.test(group1, mu=mean12)
t.test(group2, mu=mean12)
t.test(group3, mu=mean12)

################################################################################################################################################################
# One Sample t-test inequality -- i.e. you want to know if the sample mean is less than a hypothesized population mean

t.test(group3, mu=mean12, alternative="greater")

################################################################################################################################################################
# Wilcoxon single sample signed rank test -- i.e. you want to know if the sample mean is equal to a hypothesized population mean

wilcox.test(group1, mu=mean12)

################################################################################################################################################################
# Wilcoxon two sample independent signed rank test -- i.e. when the measurements in each group are related pairwise.  This test is also known as the
# "independent 2-group Mann-Whitney U Test".  T-test above.

wilcox.test(group1, group2)

################################################################################################################################################################
# Wilcoxon two sample paired signed rank test -- i.e. when the measurements in each group are related pairwise.  See the paired T-test above.

wilcox.test(group1, group2, paired=TRUE)

################################################################################################################################################################
## The test functions have a formula interface.  Be careful to make sure that your grouping factor has precisely 2 levels.

t.test(v ~ group, subset(allDat, group %in% c('group1', 'group2')))


