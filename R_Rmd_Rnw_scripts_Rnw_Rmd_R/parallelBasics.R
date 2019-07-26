# -*- Mode:R; Coding:utf-8; fill-column:160 -*-

################################################################################################################################################################
# @file      parallelBasics.R
# @author    Mitch Richling <https://www.mitchr.me>
# @Copyright Copyright 2015 by Mitch Richling.  All rights reserved.
# @brief     Basic use of the parallel package bundled with R.@EOL
# @Keywords  base r package cran parallel
#

################################################################################################################################################################
# Load up the library
library(parallel)
 
################################################################################################################################################################
# Calculate the number instances based on core count
instancesInMyComputeCluster <- detectCores() - 1
 
################################################################################################################################################################
# Start up the cluster (FORK only works on UNIX'ish OSes)
myComputeCluster <- makeCluster(instancesInMyComputeCluster, type="FORK")

################################################################################################################################################################
# Run a function in each instance (you can add arguments after the function name)
unlist(clusterCall(myComputeCluster, 'Sys.getpid'))

################################################################################################################################################################
# Run an expression in each instance (not necessary for the code below -- just showing how to run something)
unlist(clusterEvalQ(myComputeCluster, Sys.getpid()))

################################################################################################################################################################
# Create some data, and compute the sin of each element.  Notice aVar need not be exported.
aVar <- 1:10
parSapply(myComputeCluster, aVar, sin)
parLapply(myComputeCluster, aVar, sin)

################################################################################################################################################################
# Create another vector.  Note that it must be exported as it is in the last arg of parSapply.
bVar <- 10
clusterExport(myComputeCluster, 'bVar')
parSapply(myComputeCluster, aVar, function (x) bVar*x)

################################################################################################################################################################
# Create some big data data and put it in 'cVar'
cVar <- runif(1000)

################################################################################################################################################################
# Export 'cVar' to to each instance
clusterExport(myComputeCluster, "cVar")

################################################################################################################################################################
# Load microbenchmark in the 'master' instance -- it will not be in the other instances
library(microbenchmark)

################################################################################################################################################################
# Benchmark computing 'sin' on our vector.  We do it on the cluster vs just the master instance
microbenchmark(parSapply(myComputeCluster, cVar, sin),
               sapply(cVar, sin))

################################################################################################################################################################
# Shut down cluster
stopCluster(myComputeCluster)

