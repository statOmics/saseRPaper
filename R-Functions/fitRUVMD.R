
#' @title The grid search for the optimal number of latent factors

#'

#' @description

#'

#' findEncodingDimRUV is used to optimise the number of latent factors included

#' in the final regression to search for aberrant expression. It introduces

#' corrupted counts, after which multiple analyses are done at different

#' number of latent factors included. The number of latent factors for which

#' the highest area under the precision-recall curve is obtained, is defined

#' as the optimal number of latent factors.

#' @param se A `SummarizedExperiment` instance generated with the
#' SummarizedExperiment function of the SummarizedExperiment package.
#' In the assay slot, provide the expression counts as an
#' ordinary `matrix`, `DataFrame`, a `sparseMatrix` or a `DelayedMatrix`.
#' `colData` is a `DataFrame` describing the samples in the experiment.
#' Finally, specify the experimental design as a formula in the metadata slot
#' with design as name. This formula must be based on the colData.

#' @param dimensions A vector containing the number of latent factors that are
#' compared in the grid search for the optimal number of latent factors.

#' @param freq the frequency of corrupted counts injected in the dateset.

#' @param zScore The mean magnitude of the corrupted counts in the dataset.

#' @param sdlog Standard deviation of the distribution of the corrupted counts
#' on the log-scale.

#' @param lnorm If TRUE, the corrupted counts are simulated from a log-normal
#' distribution with a mean of log(zScore) and standard deviation of sdlog.

#' @param inj Injection of overexpression ('high'), underexpression ('low') or
#' both ('both') types of corrupted counts.

#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#' back-end to be used for computations. See
#' \code{bpparam} in \code{BiocParallel} package for details.


#' @return An updated 'SummarizedExperiment' instance, now including a new
#' data.table ('encDimTable') and number ('optimalEncDim') in the metadata slot.
#' @import edgeR
#' @import MASS
#' @import pracma
#' @import SummarizedExperiment
#' @import precrec
#' @import PRROC
#' @import BiocGenerics
#' @import S4Vectors
#' @importFrom limma lmFit
#' @importFrom data.table data.table
#' @importFrom DESeq2 DESeqDataSet
#' @importFrom BiocParallel bplapply bpparam
#'
#' @examples
#' CountMatrix <- matrix(rnbinom(n = 1000, mu = 10, size = 5), ncol = 20, nrow = 50)
#' Metadata <- c(rep(c("Male","Female"), each = 10))
#' se <- SummarizedExperiment::SummarizedExperiment(
#' assays = list("counts" = CountMatrix),
#'  colData = data.frame("Sex" = Metadata),
#'  metadata = list(design = ~Sex))
#' se <- findEncodingDimRUV(se = se, dimensions = c(0:10))
#' se <- fitRUVMD(se = se)
#' SummarizedExperiment::assays(se)[["mu"]][1:5,1:5]
#' SummarizedExperiment::assays(se)[["pValue"]][1:5,1:5]
#' SummarizedExperiment::assays(se)[["pValueAdjust"]][1:5,1:5]
#' SummarizedExperiment::rowData(se)[["theta"]][1:5]
#' S4Vectors::metadata(se)[["optimalEncDim"]]
#'
#'
#' @export




findEncodingDimRUV <- function(se,
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


  #Orthonormalisation

  orthonormal_vector_space <- .orthonormalisation(se=se,
                                        svd_latent_factors=svd_latent_factors)$Q

  # Calculate for each number of latent factors (= dimensions) that is
  # given as input the area under the curve of all merged p-values when
  # doing a RUV analysis with as true positives the corrupted counts.
  # This grid search is done to obtain the optimal number of latent factors
  # to include in the model.

  beta_initial <- .initial_beta(se=se,
                              orthonormal_vector_space=orthonormal_vector_space,
                              offset=fit_DGE$offset)

  eval <- bplapply(X=dimensions, ..., BPPARAM=BPPARAM,
                   FUN=function(i, ..., .evalAucPRLoss=NA){
                     .evalAutoCorrection(se,
                            orthonormal_vector_space = orthonormal_vector_space,
                            dimensions=i,
                            beta_initial,
                            fit_DGE,
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





#' @title Searching for aberrant expression

#'

#' @description

#'

#' fitRUVMD is used to perform a negative binomial regression and to search
#' for aberrant expression in the context of rare Mendelian disorders.
#' It uses a predefined or an optimised with the grid search number of latent
#' factors in the regression. It returns the p-values for the prioritisation
#' of aberrant expression.

#' @param se A `SummarizedExperiment` instance generated with the
#' SummarizedExperiment function of the SummarizedExperiment package.
#' In the assay slot, provide the expression counts as an
#' ordinary `matrix`, `DataFrame`, a `sparseMatrix` or a `DelayedMatrix`.
#' `colData` is a `DataFrame` describing the samples in the experiment.
#' Finally, specify the experimental design as a formula in the metadata slot.
#' This formula must be based on the colData. If 'optimalEncDim' is available
#' in the metadata slot, this number will be used as the number of latent
#' factors.

#' @param dimensions A single integer representing the number of latent factors
#' that are used in the regression framework. If not specified, the
#' 'optimalEncDim' in the metadata will be used. If nor 'optimalEncDim' is
#' specified, the number of samples divided by 4 will be used.

#' @param padjust the method used to compute FDR-adjusted p-values. Both
#' Benjamini-Hochberg ('BH') or Benjamini-Yekutieli ('BY') are recommended.

#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#' back-end to be used for computations. See
#' \code{bpparam} in \code{BiocParallel} package for details.

#' @return An updated 'SummarizedExperiment' instance, now including a matrix
#' of p-values ('pValue') and a matrix of FDR-adjusted p-values
#' ('pValueAdjust') in the assay slot. Also 'pValueAdjustMethod' is available
#' in the metadata slot, stating which method was used to obtain FDR-adjusted
#' pvalues.
#'
#' @import edgeR
#' @import MASS
#' @import pracma
#' @import SummarizedExperiment
#' @import precrec
#' @import PRROC
#' @import BiocGenerics
#' @import S4Vectors
#' @importFrom limma lmFit
#' @importFrom data.table data.table
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom DESeq2 DESeqDataSet
#'
#' @examples
#'
#' CountMatrix <- matrix(rnbinom(n = 1000, mu = 10, size = 5), ncol = 20, nrow = 50)
#' Metadata <- c(rep(c("Male","Female"), each = 10))
#' se <- SummarizedExperiment::SummarizedExperiment(
#' assays = list("counts" = CountMatrix),
#'  colData = data.frame("Sex" = Metadata),
#'  metadata = list(design = ~Sex))
#' se <- findEncodingDimRUV(se = se, dimensions = c(0:10))
#' se <- fitRUVMD(se = se)
#' SummarizedExperiment::assays(se)[["mu"]][1:5,1:5]
#' SummarizedExperiment::assays(se)[["pValue"]][1:5,1:5]
#' SummarizedExperiment::assays(se)[["pValueAdjust"]][1:5,1:5]
#' SummarizedExperiment::rowData(se)[["theta"]][1:5]
#' S4Vectors::metadata(se)[["optimalEncDim"]]
#'
#' @export

fitRUVMD <- function(se, dimensions = NULL, padjust = "BH", BPPARAM = bpparam()){

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


  #Orthonormalisation
  QR <-.orthonormalisation(se=se,
                           svd_latent_factors=svd_latent_factors)
  metadata(se)$Q <- QR$Q
  metadata(se)$R <- QR$R

  orthonormal_vector_space <- QR$Q


  # Calculate for each number of latent factors (= dimensions) that is
  # given as input the area under the curve of all merged p-values when
  # doing a RUV analysis with as true positives the corrupted counts.
  # This grid search is done to obtain the optimal number of latent factors
  # to include in the model.

  beta_initial <- .initial_beta(se=se,
                              orthonormal_vector_space=orthonormal_vector_space,
                              offset=fit_DGE$offset)

  se <- .fitRUV(se = se, orthonormal_vector_space = orthonormal_vector_space,
                 dimensions = dimensions, beta_initial = beta_initial,
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

  if(!(is.null(design) | design == ~1)){
    design <- model.matrix(getDesign(se), data = colData(se))
  } else{
    design <- model.matrix(~1, data = colData(se))
  }

  DGE <- estimateDisp(DGE, design = design)
  return(DGE)
}

.orthonormalisation <- function(se, svd_latent_factors){
  if (getDesign(se) != ~1){
    env <- environment()
    formula <- update(getDesign(se), ~ . + svd_latent_factors, env = env)
    environment(formula) <- env
    design <- model.matrix(formula, data = colData(se))
  } else {
    design <- model.matrix(~1+svd_latent_factors, data = colData(se))
  }
  qr <- gramSchmidt(A = design)

  return(qr)

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

.initial_beta <- function(se, orthonormal_vector_space, offset = NULL){
  if(is.null(offset)){
    offset <- matrix(0,ncol = ncol(se), nrow = nrow(se))
  }

  log_data <- log((counts(se)+1))-offset
  beta_initial<- lmFit(log_data, design = orthonormal_vector_space)$coefficients

  return(beta_initial)
}

# Adapted from OUTRIDER - Brechtmann et al.
.evalAutoCorrection <- function(se, orthonormal_vector_space,
                                dimensions, beta_initial, fit_DGE, BPPARAM,...){
  # Search for corrupted counts by RUV with simultaneous parameter estimation
  se <- .fitRUV(se = se, orthonormal_vector_space = orthonormal_vector_space,
                       dimensions = dimensions, beta_initial = beta_initial,
                fit_DGE = fit_DGE, ...)
  eloss <- .evalAucPRLoss(se)
  print(paste0('Evaluation loss: ', eloss,' for q=',dimensions))
  return(eloss)
}



.fitRUV <- function(se, orthonormal_vector_space, dimensions, beta_initial, fit_DGE){

  # Extract the number of latent factors that is estimated to be optimal
  # by prior knowledge or hyperparameter optimisation.
  number_of_known_confounders <- dim(model.matrix(getDesign(se), data = colData(se)))[2]
  total_dimensions <- dimensions + number_of_known_confounders

  reduced_orthonormal_vector_space <- orthonormal_vector_space[,1:total_dimensions]
  beta_initial_reduced_dimensions <- beta_initial[,1:total_dimensions]

  # Form a design matrix with the known covariates and latent factors.
  # Then, an edgeR regression is performed with these known covariates
  # and latent factors.


  mu <- .adapted_IRLS(se = se,
                      initial_beta=beta_initial_reduced_dimensions,
                      reduced_orthonormal_vector_space =
                        reduced_orthonormal_vector_space,
                      fit_DGE = fit_DGE)
  assay(se, 'mu', withDimnames=FALSE) <- mu


  DGE <- .fitEdgeRDisp(se=se, offset=mu, design=~1)
  theta <- 1/DGE$tagwise.dispersion
  rowData(se)$theta <- theta


  #Calculate 2-sided p-values for every count based on previous estimates.
  se <- .calcPValues(se, mu = mu, theta = theta)

  # Store the p-values, estimated means and dispersion to the
  # OutriderDataSet.

  return(se)


}


.adapted_IRLS <- function(se, initial_beta, reduced_orthonormal_vector_space, fit_DGE, maxiter=50){

  mu <- exp(t(reduced_orthonormal_vector_space%*%t(initial_beta)) + fit_DGE$offset)

  for(i in c(1:maxiter)){
    z <- log(mu) - fit_DGE$offset + counts(se)/mu - 1

    beta <- t(reduced_orthonormal_vector_space)%*%t(z)
    mu <- pmin(exp(t(reduced_orthonormal_vector_space%*%beta) + fit_DGE$offset)*0.2 + 0.8* mu,
               matrix(.Machine$double.xmax,ncol = ncol(se),nrow = nrow(se)))
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

