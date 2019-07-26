


setwd("/project/umb_triley/cpct/rna-seq/urine1/cuffdiff_results_GRCh38_gtf_only/LUTS-over-CTRL")

library(glmnet)
nested.cvGlmnet <- function(x.train, y.train, num.iter = 10, nfold = 4, verbose = TRUE){
  
  n <- length(y.train)
  pred <- matrix(NA, nrow = n, ncol = num.iter)
  best.lambdas <- matrix(0, nrow = (nfold+1), ncol = 2*num.iter)
  
  size = as.integer(floor(n/(nfold+1)))
  ind.split.outer = c(seq.int(0L,by=size,len=(nfold+1)),n)
  outer.indecies <- lapply(1:(nfold+1), function(x) list()) # predicted responses 		
  outer.indecies <- lapply(1:num.iter, function(x) outer.indecies) ## repeate the list for all itterations
  
  
  labs <- matrix(NA, nrow = n, ncol = num.iter)
  nonzero.genes = {}
  nonzero.coeffs = {}
  ## generating CV folds.
  for(iter in 1:num.iter){
    if(verbose){
      cat('\n\n')
      cat(paste('running outer iteration', iter))
      cat('\n')
    }
    
    
    while(TRUE){
      Flag = TRUE
      ind <- sample(1:n)
      for (o in 1:(nfold+1))
      {
        ind.outer <- ind[(ind.split.outer[o]+1):ind.split.outer[o+1]]
        ind.inner <- setdiff(ind,ind.outer)
        if(sum(y.train[ind.inner]==0) <= 3 || sum(y.train[ind.inner]==1) <= 3){
          Flag = FALSE
          break
        }
      }
      
      if(Flag)
        break
    }
    for (o in 1:(nfold+1))
    {
      if (verbose)
        cat("\n*** OUTER NFOLD", o, "***")
      ind.outer <- ind[(ind.split.outer[o]+1):ind.split.outer[o+1]]
      ind.inner <- setdiff(ind,ind.outer)
      outer.indecies[[iter]][[o]] = list(ind.outer)
      if (verbose)
        cat("\n*** Running CRE ***")		
      inner.data <- list(x = x.train[ind.inner,], y = y.train[ind.inner])
      
      ## Cross-validation
      if (verbose)
        cat("\n*** Running Cross-Validation ***")
      ## looping through the inner folds ans selecting the best model (lambda)
      cv.fit <- NULL
      while(is.null(cv.fit)){
        try(cv.fit<-cv.glmnet(inner.data$x, inner.data$y, family="binomial", type.measure = "deviance", nfolds = nfold))
      }
      #cv.fit <- cv.glmnet(inner.data$x, inner.data$y, family="binomial", type.measure = "deviance", nfolds = nfold)
      lambdas <- cv.fit$lambda
      lambda.min <- cv.fit$lambda.min
      best.lam.ind <- which(lambdas == lambda.min)
      best.lambdas[o, (2*iter-1):(2*iter)] <- c(lambda.min, best.lam.ind)
      
      ## fir on all folds 
      inner.fit <- glmnet(inner.data$x, inner.data$y, family="binomial", alpha = 0.47,lambda = lambdas)
      ##x.outer = standardize(x.train[ind.outer, ,drop = F])$x
      x.outer = x.train[ind.outer, ,drop = F]
      pred[ind.outer, iter] <- predict(inner.fit, newx = x.outer, s = lambda.min, type = "response", mode = "lambda")
      labs[ind.outer, iter] <- ifelse(pred[ind.outer, iter] > 0.5, 1, 0)
      
      
      ##
      nonzero.genes = c(nonzero.genes, colnames(x.train[,which(inner.fit$beta[, best.lam.ind] != 0)]))
      nonzero.coeffs = c(nonzero.coeffs, inner.fit$beta[, best.lam.ind][which(inner.fit$beta[, best.lam.ind] != 0)])
    }
    ##nonzero.genes = unique(nonzero.genes)
  }
  
  L <- list(pred = pred, labs = labs, best.lambdas = best.lambdas, nonzero.genes = nonzero.genes, 
            nonzero.coeffs = nonzero.coeffs, outer.indecies = outer.indecies)
  
  class(L) = c("nested.cvGlmnet")
  
  return(L)
  
}
load("SampleByGene-gtf-only.rda")
result=nested.cvGlmnet(sigGenes2,c(rep(1,9),rep(0,9)) )
