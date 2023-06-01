library(BiocParallel)
library(data.table)
library(PRROC)
library(pracma)


findEncodingDim_RUV_edgeR <- function(se,
                                dimensions=seq(
                                  2, min(100, ncol(se) - 2, nrow(se) - 1), 2),
                                freq=1E-2,
                                zScore=3,
                                sdlog=log(1.6),
                                lnorm=TRUE,
                                inj='both',
                                BPPARAM=bpparam(),
                                ...){
  # Find the optimal number of hyperparameters when using RUV with 
  # deviance residuals.
  # Estimate size factors
  if (is.null(getSizeFactors(se))){
    se <- .estimateSizeFactors(se)
  }
  # Inject corrupted counts
  se <- .injectOutliersHyperparam(se, 
                                 freq=freq,
                                 zScore=zScore, 
                                 inj=inj, 
                                 lnorm=lnorm,
                                 sdlog=sdlog)
  
  # Fit edgeR model on intercept and known covariates. This is
  # further used to calculate the deviances.
  DGE <- .fitEdgeRDisp(se = se, design = getDesign(se))
  fit_DGE <- glmFit(y = DGE)
  
  # Principal component analysis on the deviance residuals
  svd_latent_factors <- .svd_deviances(se = se,
                                       fit_DGE = fit_DGE)[,1:max(dimensions)]
  
  eval <- bplapply(X=dimensions, ..., BPPARAM=BPPARAM, 
                   FUN=function(i, ..., .evalAucPRLoss=NA){
                     .evalAutoCorrection(se,
                                         svd_latent_factors = svd_latent_factors,
                                         dimensions=i, 
                                         fit_DGE = fit_DGE,
                                         BPPARAM=BPPARAM)})
  
  # Input the grid search results in the OutriderDataSet
  metadata(se)[['encDimTable']] <- data.table(
    encodingDimension= dimensions,
    evaluationLoss= unlist(eval), 
    evalMethod='aucPR')
  
  #Obtain the optimal number of latent factors
  metadata(se)[['optimalEncDim']] <- NULL
  metadata(se)[['optimalEncDim']] <- .getBestQ(se)
  
  # Input the original counts again as count matrix, as some values were
  # changed with corrupted counts.
  assay(se,'counts', withDimnames = F) <- assay(se, 'trueCounts')
  return(se)
}

fitRUVMD_edgeR <- function(se, dimensions = NULL, padjust = "BH", BPPARAM = bpparam()){
  
  if(is.null(dimensions)){
    if(is.null(metadata(se)[['optimalEncDim']])){
      dimensions <- round(ncol(se)/4,0)
    }
    else{
      dimensions <- metadata(se)[['optimalEncDim']]
    }
  }
  
  DGE <- .fitEdgeRDisp(se = se, design = getDesign(se))
  fit_DGE <- glmFit(y = DGE)
  
  # Principal component analysis on the deviance residuals
  svd_latent_factors <- .svd_deviances(se = se,
                                       fit_DGE = fit_DGE)[,1:dimensions]
  

  
  
  # Calculate for each number of latent factors (= dimensions) that is
  # given as input the area under the curve of all merged p-values when 
  # doing a RUV analysis with as true positives the corrupted counts.
  # This grid search is done to obtain the optimal number of latent factors
  # to include in the model.
  

  se <- .fitRUV_edgeR(se = se, svd_latent_factors = svd_latent_factors,
                dimensions = dimensions, 
                fit_DGE = fit_DGE)
  
  se <- .adjustPValues(se, method = padjust)
  
  return(se)
  
}

.fitEdgeRDisp <- function(se, offset=NULL, design=~1, ...){
  # Fit edgeR model on intercept and known covariates. This is
  # further used to calculate the deviances.
  
  DGE <- DGEList(counts=counts(se))
  if(is.null(offset)){
    DGE <- calcNormFactors(DGE)
  } else {
    DGE$offset <- log(offset)
  }
  
  if(!(design == ~1)){
    if(!("matrix" %in% class(model.matrix(~1, colData(se))))){
      design <- model.matrix(design, data = colData(se))
    }
  } else{
    design <- model.matrix(~1, data = colData(se))
  }
  
  DGE <- estimateDisp(DGE, design = design) 
  return(DGE)
}





.svd_deviances <- function(se, fit_DGE){
  # Returns the latent factors of the deviance residuals to be used in 
  # parameter estimation of RUV
  
  deviance <- nbinomUnitDeviance(counts(se),
                                 mean = fit_DGE$fitted.values,
                                 dispersion = fit_DGE$dispersion)
  deviance_scaled <- sqrt(abs(deviance))*sign(counts(se)-fit_DGE$fitted.values)
  
  svd_latent <- svd(t(deviance_scaled))
  return(svd_latent$u)
}



.evalAutoCorrection <- function(se, svd_latent_factors,
                                dimensions, beta_initial, fit_DGE, BPPARAM,...){
  # Search for corrupted counts by RUV with simultaneous parameter estimation
  se <- .fitRUV_edgeR(se = se, svd_latent_factors = svd_latent_factors,
                dimensions = dimensions,
                fit_DGE = fit_DGE, ...)
  eloss <- .evalAucPRLoss(se)
  print(paste0('Evaluation loss: ', eloss,' for q=',dimensions))
  return(eloss)
}



.fitRUV_edgeR <- function(se, svd_latent_factors, dimensions, fit_DGE){
  # Extract the number of latent factors that is estimated to be optimal
  # by prior knowledge or hyperparameter optimisation.
  total_dimensions <- dimensions
  
  reduced_svd_latent_factors <- svd_latent_factors[,1:total_dimensions]

  # Form a design matrix with the known covariates and latent factors.
  # Then, an edgeR regression is performed with these known covariates
  # and latent factors.

  if (getDesign(se) != ~1){
      env <- environment()
      formula <- update(getDesign(se), ~ . + reduced_svd_latent_factors, env = env)
      environment(formula) <- env
      design <- model.matrix(formula, data = colData(se))
  } else {
       design <- model.matrix(~1+reduced_svd_latent_factors, data = colData(se))
  }

  DGE <- .fitEdgeRDisp(se=se, design=design)
  fit_DGE <- glmFit(DGE, design = design)
  
  mu <- fit_DGE$fitted.values
  theta <- 1/fit_DGE$dispersion
  
  assay(se,"mu",withDimnames = FALSE) <- mu
  rowData(se)$theta <- theta
  
  
  #Calculate 2-sided p-values for every count based on previous estimates.
  se <- .calcPValues(se, mu = mu, theta = theta)
  
  # Store the p-values, estimated means and dispersion to the
  # OutriderDataSet.
  
  return(se)
  
  
}


.adapted_IRLS <- function(se, initial_beta, reduced_orthonormal_vector_space, fit_DGE, maxiter=50){ #offset slot van maken; maxiter toevoegen(standaard 50); bekijken of wel time gain met deviance residuals of andere stop
  
  mu <- exp(t(reduced_orthonormal_vector_space%*%t(initial_beta)) + fit_DGE$offset)
  deviance <- sum(nbinomUnitDeviance(counts(se),
                                     mean=mu,
                                     dispersion=fit_DGE$dispersion))
  
  for(i in c(1:maxiter)){
    z <- log(mu+0.001) - fit_DGE$offset + 
      pmax(
        pmin(counts(se)/mu,matrix(10,ncol = ncol(se),nrow = nrow(se))),
        matrix(0.1,ncol = ncol(se),nrow = nrow(se))
      ) - 1
    
    beta <- t(reduced_orthonormal_vector_space)%*%t(z)
    mu_previous <- mu 
    mu <- exp(t(reduced_orthonormal_vector_space%*%beta) + fit_DGE$offset)*0.2 +
      0.8* mu_previous
    deviance_previous <- deviance
    deviance <- sum(nbinomUnitDeviance(counts(se),
                                       mean = mu, 
                                       dispersion = fit_DGE$dispersion))
    print(deviance)
    if(deviance_previous - deviance < 0){
      mu <- mu_previous
      return(mu)
    }
  }
  return(mu)
}

.calcPValues <- function(se, mu, theta){
  left_quantile <- pnbinom(q = counts(se), mu = mu, size = theta)
  right_quantile <- pnbinom(q = counts(se)-1, mu = mu, size = theta)
  PValues <- pmin(left_quantile*2, (1-right_quantile)*2,1)
  assay(se, 'pValue', withDimnames=FALSE) <- PValues
  return(se)
}


.adjustPValues <- function(se, method = "BH"){
  assay(se, 'pValueAdjust', withDimnames=FALSE) <- 
    apply(X = assay(se, 'pValue', withDimnames=FALSE), MARGIN = 2, FUN = p.adjust, method = method)
  metadata(se)[['pValueAdjustMethod']] <- method
  return(se)
}

# Adapted from OUTRIDER - Brechtmann et al.
.estimateSizeFactors <- function(se, design = ~1){
  dse <- DESeqDataSet(se = se, design = design)
  dse <- DESeq2::estimateSizeFactors(dse)
  colData(se, 'sizeFactor', withDimnames=FALSE) <- 
    colData(dse, 'sizeFactor', withDimnames=FALSE)
  return(se)
}


# Adapted from OUTRIDER - Brechtmann et al.
.getBestQ <- function(ods){
  if('optimalEncDim' %in% names(metadata(ods))){
    return(metadata(ods)[['optimalEncDim']])
  }
  
  if('encDimTable' %in% names(metadata(ods))){
    encTable <- metadata(ods)[['encDimTable']]
    return(.getBestQDT(encTable, 'aucPR'))
  }
  # warning('Please find the optimal encoding dimension by running. ')
  return(NA_integer_)
}

# Adapted from OUTRIDER - Brechtmann et al.
.getBestQDT <- function(dt, usedEvalMethod='aucPR', digits=10){
  # if('evalMethod' %in% colnames(dt)){
  #     testFun <- ifelse(all(dt[,evalMethod == usedEvalMethod]), 
  #             which.max, which.min)
  # } else {
  testFun <- which.max
  #}
  
  dt[,encodingDimension[
    seq_len(.N) == testFun(round(evaluationLoss, digits))]]
}

# Adapted from OUTRIDER - Brechtmann et al.
.evalAucPRLoss <- function(se){
  # Calculate the precision-recall curve of all p-values, using corrupted
  # counts as true positives.
  scores <- -as.vector(assay(se, 'pValue'))
  labels <- as.vector(assay(se, 'trueCorruptions') != 0) + 0
  
  if(any(is.na(scores))){
    warning(sum(is.na(scores)), " P-values where NAs.")
    scores[is.na(scores)] <- min(scores, na.rm=TRUE)-1
  }
  pr <- pr.curve(scores, weights.class0=labels)
  return(max(0, pr$auc.integral, na.rm=TRUE))
}

