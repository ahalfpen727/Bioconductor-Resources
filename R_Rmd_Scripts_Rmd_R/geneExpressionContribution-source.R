VlsVt <- function(data,rescale=TRUE){
  tol <- 10^-14
  k   <- dim(data)[2] # number of individuals
  n   <- dim(data)[1] # number of transcripts
  unC <- (k-1) / k    # to uncorrect
  if(rescale)
    data = sqrt(data)
  sigma <- unC * cov(t(data)) # corrected covariance
  ## Vt: Total variance
  Vt = sum(diag(sigma))
  ## Vls: variance in the model of constant splicing
  conjecture    <- svd(data,1,0)$u[,1]   
  indexP   <- conjecture > tol
  indexN   <- conjecture < -tol    
  nPos     <- sum(indexP)
  nNeg     <- sum(indexN)    
  if ((nPos > 0) & (nNeg>0)) { # Perron-Frobenius fails?
    require(NMFN)
    matrixFactor <- nnmf(data,1,method = "nnmf_als")
    bestDirection <- matrixFactor$W /  sqrt((t(matrixFactor$W) %*% matrixFactor$W)[1,1])
    Vls = ( t(bestDirection) %*% sigma ) %*% bestDirection
  }
  else {
    if (nNeg > 0) conjecture <- -1*conjecture 
    Vls = ( t(conjecture) %*% sigma ) %*% conjecture
  }
  ## Vls2: variance in the model of constant gene expression
  vecModel = vec = as.numeric(unlist(data))       # vectorize
  ## directament com OLS a veure si hi ha sort
  bigSum <- sum(vec)  
  for (i in 1:k) {
    littleSum                   <- k*sum(vec[(1+(i-1)*n):(i*n)])  
    vecModel[(1+(i-1)*n):(i*n)] <- vec[(1+(i-1)*n):(i*n)] + (bigSum - littleSum)/(n*k)
  }  
  nCompNeg    <- sum(vecModel < 0)
  if (nCompNeg == 0) data2g <- matrix(vecModel,nrow=n,ncol=k) 
  if (nCompNeg > 0) {
    Dg    <- matrix(0,nrow=(n*k),ncol=(n*k))
    Amatg <- matrix(0,nrow=(n*k+k-1),ncol=n*k)
    diag(Dg) = 1
    Amatg[1:(k-1),1:n] <- 1
    for (i in 1:(n*k)) {
      Amatg[(k-1+i),i] <- 1
      }
    for (i in 2:k) {
      Amatg[(i-1),(((i-1)*n)+1):(i*n)] <- -1
    }
    Amatg <- t(Amatg)
    dg    <- vec
    require(quadprog)
    sQP  <- solve.QP(Dg,dg,Amatg,meq=(k-1))
    Y2g         <- sQP$solution
    data2g       <- matrix(Y2g,nrow=n,ncol=k)
  }
  sigma2g    <- unC * cov(t(data2g)) ## Uncorrected variance
  trSigma2g  <- sum(diag(sigma2g))
  Vls2 = trSigma2g
  return(list(vls.vt=Vls/Vt,vls2.vt=Vls2/Vt))
}

VlsVt.bootstrap <- function(data,rescale=TRUE,nb.bs=100,prop.bs=.5){
  if(!is.null(nb.bs)){
    require(plyr)
    vls.bs = ldply(1:nb.bs,function(ii){
      dat = data[,sample(1:ncol(data),ncol(data)*prop.bs,TRUE)]
      vlsvt.g = VlsVt(dat,rescale=rescale)
      data.frame(vls.vt = vlsvt.g$vls.vt,
                 vls2.vt = vlsvt.g$vls2.vt)
    })
    return(list(vls.vt=median(vls.bs$vls.vt),vls2.vt=median(vls.bs$vls2.vt)))
  } else {
    return(VlsVt(data,rescale=rescale))
  }
}

## Compute VlsVt estimator on from RPKM values across a data.frame.
## Input:
##   data: a data.frame with:
##            col1 = transcript ID; col2 = gene ID; col3 = RPKM sample 1; col4 = RPKM sample 2; ...
##            row = transcripts
##   min.iso2.RPKM: minimum expression for the second major transcript. Particularly relevant to filter out genes expressing almost only one transcript. Eventually post-computation filtering is also possible.
##   rescale: recent improvement to counter-balance a bias toward gene expression contribution in some situations.
##   nb.bs: number of bootstrap replicates. Bootstrapping is particularly useful when outlier are present(heterogeneous samples). If 0 or NULL, bootstrap is not used.
##   prop.bs: proportion of sample used in the bootstrap replicates. If enough samples, 0.5 is good.
##   verbose: verbose mode, displaying gene ID.
## Output:
##   vls.vt: a data.frame with gene name, vls.vt(i.e. gene expression contribution estimate), vls2.vt(i.e. splicing contribution extimate), iso2.RPKM(i.e. the average RPKM expression of the second major isoform).
##   rescale, nb.bs, prop.bs: parameters used for the computation.
VlsVt.comp <- function(data,min.iso2.RPKM=.01,rescale=TRUE,nb.bs=100,prop.bs=.5,verbose=FALSE){
  require(plyr)
  genes = unique(as.character(data[,2]))
  res = ldply(genes,function(gene){
    if(verbose) cat(gene,"...\n")
    if(length((gene.ids = which(data[,2]==gene)))>1){ ## At least two isoforms
      dat.g = data[gene.ids,-(1:2)]
      iso.RPKM.mean = rowMeans(dat.g,na.rm=TRUE)
      iso2.RPKM = sort(iso.RPKM.mean,decreasing=TRUE)[2]
      vls.vt = vls2.vt = NA
      if(iso2.RPKM >= min.iso2.RPKM){ ## Second major isoform expressed
        vls.vt.o = VlsVt.bootstrap(dat.g,rescale=rescale,nb.bs=nb.bs,prop.bs=prop.bs)
        vls.vt = vls.vt.o$vls.vt
        vls2.vt = vls.vt.o$vls2.vt
      }
      return(data.frame(gene=gene,
                        vls.vt=vls.vt,
                        vls2.vt=vls2.vt,
                        iso2.RPKM=iso2.RPKM))
    }
  })
  return(list(vls.vt=res,rescale=rescale,nb.bs=nb.bs,prop.bs=prop.bs))
}
