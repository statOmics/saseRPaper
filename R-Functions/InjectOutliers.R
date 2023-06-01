
# Injection of corrupted counts based on a negative binomial distribution
# for the search of the optimal number of latent dimensions. 
# Adapted from Outrider - Brechtmann et al.
.injectOutliersHyperparam <- function(se, freq, zScore, inj, lnorm, sdlog){
  
  # copy true counts to be able to acces them later
  assay(se, 'trueCounts', withDimnames=FALSE) <- counts(se)
  
  # generate index of injected corrupted counts
  size <- prod(dim(se))
  index <- sample(c(0,1,-1), size, prob = c(1 - freq, freq/2, freq/2), 
                  replace = TRUE)
  index <- matrix(index, nrow = nrow(se))
  switch(inj,
         low = { index <- -abs(index) },
         high = { index <- abs(index) }
  )
  
  # Generate z-values for corrupted counts according to a lognormal  
  # distribution if these are not specified
  tmpzScore <- matrix(0, ncol=ncol(se), nrow=nrow(se))
  if(isTRUE(lnorm)){
    tmpzScore[index!=0] <- rlnorm(sum(index!=0), log(zScore), sdlog=sdlog)
  } else {
    tmpzScore[index!=0] <- zScore
  }
  zScore <- tmpzScore
  
  # Define maximal count that should be injected
  max_out <- min(10*max(counts(se), na.rm=TRUE), .Machine$integer.max)
  
  # compute size factor normalized counts.
  # don't use it on the se to not influence the later calculation.
  se <- .estimateSizeFactors(se)
  sf_r <- getSizeFactors(se)
  # extract counts
  counts <- counts(se)
  # list of locations where corrupted counts have to be injected 
  list_index <- which(index != 0, arr.ind = TRUE)
  # Loop over the list_index to inject corrupted counts
  
  DGE <- DGEList(counts=counts(se))
  DGE <- calcNormFactors(DGE)
  DGE <- estimateDisp(DGE) # estimate dispersion estimates
  fit_DGE <- glmFit(y = DGE)
  
  for(i in seq_len(nrow(list_index))){
    idxCol <- list_index[i,'col']
    idxRow <- list_index[i,'row']
    
    # Extract dispersion of negative binomial regression and define
    # the count that will be injected as corrupted count 
    theta <- 1/fit_DGE$dispersion[idxRow]
  
    prob_z <- pnorm(zScore[idxRow, idxCol])
    
    if (index[idxRow, idxCol]==1){
      art_out <- qnbinom(p=prob_z, mu=counts(se)[idxRow,idxCol], size=theta)            
    }
    else {
      art_out <- qnbinom(p=1-prob_z, mu=counts(se)[idxRow,idxCol], size=theta)
    }
    
    # only insert outliers if they are different from before 
    # and not too large
    if(art_out < max_out & counts[idxRow, idxCol] != art_out){
      counts[idxRow, idxCol] <- art_out
    }else{
      index[idxRow, idxCol] <- 0
    }
  }
  
  # save the new count matrix, the index of corrupted counts and the 
  # magnitude of the corrupted count via the z-value
  assay(se, 'counts', withDimnames=FALSE) <- matrix(as.integer(counts), 
                                                     nrow=nrow(se))
  assay(se, 'trueCorruptions', withDimnames=FALSE) <- index
  assay(se, 'injectedZscore', withDimnames=FALSE) <- zScore
  return(se)
}

# Adapted from Outrider - Brechtmann et al.
injectOutliers <- function(se, freq, zScore, inj, lnorm, sdlog){
  # copy true counts to be able to acces them in the loss later
  assay(se, 'Counts_before_outliers', withDimnames=FALSE) <- counts(se)
  
  # generate index of injected counts
  size <- prod(dim(se))
  index <- sample(c(0,1,-1), size, prob = c(1 - freq, freq/2, freq/2), 
                  replace = TRUE)
  index <- matrix(index, nrow = nrow(se))
  switch(inj,
         low = { index <- -abs(index) },
         high = { index <- abs(index) }
  )
  
  tmpzScore <- matrix(0, ncol=ncol(se), nrow=nrow(se))
  if(isTRUE(lnorm)){
    tmpzScore[index!=0] <- rlnorm(sum(index!=0), log(zScore), sdlog=sdlog)
  } else {
    tmpzScore[index!=0] <- zScore
  }
  zScore <- tmpzScore
  
  #inject counts
  max_out <- min(10*max(counts(se), na.rm=TRUE), .Machine$integer.max)
  
  # compute size factor normalized counts.
  # don't use it on the se to not influence the later calculation.
  se <- .estimateSizeFactors(se)
  sf_r <- getSizeFactors(se)
  
  counts <- counts(se)
  list_index <- which(index != 0, arr.ind = TRUE)
  for(i in seq_len(nrow(list_index))){
    idxCol <- list_index[i,'col']
    idxRow <- list_index[i,'row']
    nb <-  try(glm.nb(counts(se)[idxRow,]~1+offset(log(sf_r))))
    if("try-error" %in% class(nb)) {
      index[idxRow, idxCol] <- 0
    }
    else{
      mu <- exp(nb$coefficients[[1]])
      theta <- nb$theta
      
      prob_z <- pnorm(zScore[idxRow, idxCol])
      if (index[idxRow, idxCol]==1){
        art_out <- qnbinom(p=prob_z, mu=mu*sf_r[idxCol], size=theta)            
      }
      else {
        art_out <- qnbinom(p=1-prob_z, mu=mu*sf_r[idxCol], size=theta)
      }
      
      # only insert outliers if they are different from before 
      # and not too large
      if(art_out < max_out & counts[idxRow, idxCol] != art_out){
        counts[idxRow, idxCol] <- art_out
      }else{
        index[idxRow, idxCol] <- 0
      }
    }
  }
  # save coruppted counts and index of corruption into se
  assay(se, 'counts', withDimnames=FALSE) <- matrix(as.integer(counts), 
                                                     nrow=nrow(se))
  assay(se, 'Outliers', withDimnames=FALSE) <- index
  assay(se, 'injectedZscore_Outliers', withDimnames=FALSE) <- zScore
  return(se)
}

# Adapted from Outrider - Brechtmann et al.
injectOutliers_per_confounder <- function(se, freq, zScore, inj, lnorm, sdlog, confounder_vector){
  # copy true counts to be able to acces them in the loss later
  assay(se, 'trueCounts', withDimnames=FALSE) <- counts(se)
  
  # generate index of injected counts
  size <- prod(dim(se))
  index <- sample(c(0,1,-1), size, prob = c(1 - freq, freq/2, freq/2), 
                  replace = TRUE)
  index <- matrix(index, nrow = nrow(se))
  switch(inj,
         low = { index <- -abs(index) },
         high = { index <- abs(index) }
  )
  
  tmpzScore <- matrix(0, ncol=ncol(se), nrow=nrow(se))
  if(isTRUE(lnorm)){
    tmpzScore[index!=0] <- rlnorm(sum(index!=0), log(zScore), sdlog=sdlog)
  } else {
    tmpzScore[index!=0] <- zScore
  }
  zScore <- tmpzScore
  
  #inject counts
  max_out <- min(10*max(counts(se), na.rm=TRUE), .Machine$integer.max)
  
  # compute size factor normalized counts.
  # don't use it on the se to not influence the later calculation.
  se <- .estimateSizeFactors(se)
  sf_r <- getSizeFactors(se)

  counts <- counts(se)
  list_index <- which(index != 0, arr.ind = TRUE)
  for(i in seq_len(nrow(list_index))){
    idxCol <- list_index[i,'col']
    idxRow <- list_index[i,'row']
    
    nb <-  try(glm.nb(counts(se)[idxRow,]~1+ confounder_vector + offset(log(sf_r))))
    if("try-error" %in% class(nb)) {
      index[idxRow, idxCol] <- 0
    }
    else{
      mu <- exp(predict.glm(object = nb, 
                            newdata = data.frame("confounder_vector" = confounder_vector,"sf_r" = sf_r)))
      theta <- nb$theta
      
      prob_z <- pnorm(zScore[idxRow, idxCol])
      if (index[idxRow, idxCol]==1){
        art_out <- qnbinom(p=prob_z, mu=mu[idxCol], size=theta)            
      }
      else {
        art_out <- qnbinom(p=1-prob_z, mu=mu[idxCol], size=theta)
      }
      
      # only insert outliers if they are different from before 
      # and not too large
      if(art_out < max_out & counts[idxRow, idxCol] != art_out){
        counts[idxRow, idxCol] <- art_out
      }else{
        index[idxRow, idxCol] <- 0
      }
    }
  }
  # save coruppted counts and index of corruption into se
  assay(se, 'counts', withDimnames=FALSE) <- matrix(as.integer(counts), 
                                                     nrow=nrow(se))
  assay(se, 'Outliers', withDimnames=FALSE) <- index
  assay(se, 'injectedZscore', withDimnames=FALSE) <- zScore
  return(se)
}

# Adapted from Outrider - Brechtmann et al.
injectOutliers_lognormal <- function(se, freq, zScore, inj, lnorm, sdlog){
  # copy true counts to be able to acces them in the loss later
  assay(se, 'trueCounts', withDimnames=FALSE) <- counts(se)
  
  # generate index of injected counts
  size <- prod(dim(se))
  index <- sample(c(0,1,-1), size, prob = c(1 - freq, freq/2, freq/2), 
                  replace = TRUE)
  index <- matrix(index, nrow = nrow(se))
  switch(inj,
         low = { index <- -abs(index) },
         high = { index <- abs(index) }
  )
  
  tmpzScore <- matrix(0, ncol=ncol(se), nrow=nrow(se))
  if(isTRUE(lnorm)){
    tmpzScore[index!=0] <- rlnorm(sum(index!=0), log(zScore), sdlog=sdlog)
  } else {
    tmpzScore[index!=0] <- zScore
  }
  zScore <- tmpzScore
  
  #inject counts
  max_out <- min(10*max(counts(se), na.rm=TRUE), .Machine$integer.max)
  
  # compute size factor normalized counts.
  # don't use it on the se to not influence the later calculation.
  se <- .estimateSizeFactors(se)
  sf_r <- getSizeFactors(se)
  normtable <- log2(t(t(counts(se))/sf_r)+1)
  counts <- counts(se)
  list_index <- which(index != 0, arr.ind = TRUE)
  for(i in seq_len(nrow(list_index))){
    idxCol <- list_index[i,'col']
    idxRow <- list_index[i,'row']
    
    mu <- mean(normtable[idxRow,])
    sd <- sd(normtable[idxRow,])
    if (index[idxRow, idxCol]==1){
      art_out <- round(sf_r[idxCol]*2^(mu+zScore[idxRow,idxCol]*sd),digits=0)           
    }
    else {
      art_out <- round(sf_r[idxCol]*2^(mu-zScore[idxRow,idxCol]*sd),digits=0)           
    }
    # only insert outliers if they are different from before 
    # and not too large
    if(art_out < max_out & counts[idxRow, idxCol] != art_out){
      counts[idxRow, idxCol] <- art_out
    }else{
      index[idxRow, idxCol] <- 0
    }
  }
  
  # save coruppted counts and index of corruption into se
  assay(se, 'counts', withDimnames=FALSE) <- matrix(as.integer(counts), 
                                                     nrow=nrow(se))
  assay(se, 'Outliers', withDimnames=FALSE) <- index
  assay(se, 'injectedZscore', withDimnames=FALSE) <- zScore
  return(se)
}



