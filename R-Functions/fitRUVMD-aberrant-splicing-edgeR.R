findEncodingDimRUVSplicing <- function(se,
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
  # Inject corrupted counts
  se <- .injectOutliersHyperparamSplicing(se,
                                          freq=freq,
                                          zScore=zScore,
                                          inj=inj,
                                          lnorm=lnorm,
                                          sdlog=sdlog)
  
  real_offsets <- assays(se)$offsets
  
  se <- OffsetSplicing(se)
  
  # Fit edgeR model on intercept and known covariates. This is
  # further used to calculate the deviances.
  DGE <- .fitEdgeRDisp(se = se, design = getDesign(se), assays(se)$offsets)
  fit_DGE <- glmFit(y = DGE)
  
  # Principal component analysis on the deviance residuals
  svd_latent_factors <- .svd_deviances(se = se,
                                       fit_DGE = fit_DGE)[,1:max(dimensions)]
  
  

  
  # Calculate for each number of latent factors (= dimensions) that is
  # given as input the area under the curve of all merged p-values when
  # doing a RUV analysis with as true positives the corrupted counts.
  # This grid search is done to obtain the optimal number of latent factors
  # to include in the model.
  

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
  assay(se,'offsets', withDimnames = F) <- real_offsets 
  
  return(se)
}


# Adapted from FRASER - Mertes et al.
.injectOutliersHyperparamSplicing <- function(se, freq, zScore, inj, lnorm, sdlog, deltaMin = 0.2){
  # copy true counts to be able to acces them later
  assay(se, 'trueCounts', withDimnames=FALSE) <- counts(se)

  size <- length(unique(rowData(se)$locus))*ncol(se)
  index <- sample(c(0,1,-1), size, prob = c(1 - freq, freq/2, freq/2), 
                  replace = TRUE)
  index <- matrix(index, ncol = ncol(se))
  rownames(index) <- unique(rowData(se)$locus)
  colnames(index) <- colnames(se)
  
  switch(inj,
         low = { index <- -abs(index) },
         high = { index <- abs(index) }
  )
  
  list_index <- which(index != 0, arr.ind = TRUE)
  binOutliers <- matrix(0,ncol = ncol(se),nrow = nrow(se))
  rownames(binOutliers) <- rownames(se)
  colnames(binOutliers) <- colnames(se)
  
  locusOutliers <- matrix(0, ncol = ncol(se), nrow = nrow(se))
  rownames(locusOutliers) <- rowData(se)$locus
  colnames(locusOutliers) <- colnames(se)
  
  for(i in seq_len(nrow(list_index))){
    idxCol <- list_index[i,'col']
    idxRow <- list_index[i,'row']
    
    locus_outlier <- rownames(index)[idxRow]
    
    binnames <- rownames(se)[rowData(se)$locus == locus_outlier]
    
    outlier_bin <- sample(binnames, 1)
    other_bins <- binnames[!(binnames %in% outlier_bin)]
    
    currentProportionOutlier <- assays(se)$counts[outlier_bin,idxCol] / assays(se)$offsets[outlier_bin,idxCol]
    ProportionOthers <- assays(se)$counts[other_bins,idxCol] / assays(se)$offsets[other_bins,idxCol]
    
    
    if(is.nan(currentProportionOutlier) | length(ProportionOthers) == 0){
      index[idxRow,idxCol] <- 0
    } else{ 
      
      if(index[idxRow,idxCol] == 1){
        deltaMax <- 1-currentProportionOutlier
      } else {
        deltaMax <- currentProportionOutlier
      }
      
      
      if(deltaMin > deltaMax){
        index[idxRow,idxCol] <- -index[idxRow,idxCol]
        
        if(index[idxRow,idxCol] == 1){
          deltaMax <- 1-currentProportionOutlier
        } else {
          deltaMax <- currentProportionOutlier
        }
      }
      
      DeltaProportionOutlier <- runif(n = 1, min = deltaMin, max = deltaMax)
      
      #outlier_count <- round((currentProportion + DeltaProportionOutlier*index[idxRow,idxCol])*(assays(se)$offsets[outlier_bin,idxCol]+2)-1,0)
      if(currentProportionOutlier != 1){
        DeltaOthers <- -DeltaProportionOutlier*index[idxRow,idxCol] * (ProportionOthers/(1-currentProportionOutlier))
      } else {
        DeltaOthers <- -DeltaProportionOutlier*index[idxRow,idxCol] / rep(1/length(other_bins), times = length(other_bins))
        
      }

      #other_counts <- round((ProportionOthers + DeltaOthers)*(assays(se)$offsets[other_bins,idxCol]+2)-1,0)
      if((currentProportionOutlier + DeltaProportionOutlier*index[idxRow,idxCol] < 0) ||
         ( ProportionOthers + DeltaOthers < 0)){
      }
      injectedCounts <- rmultinom(n = 1, 
                                  size = assays(se)$offsets[other_bins,idxCol],
                                  prob = c(currentProportionOutlier + DeltaProportionOutlier*index[idxRow,idxCol],
                                           ProportionOthers + DeltaOthers))
      rownames(injectedCounts) <- c(outlier_bin,other_bins)
      
      assays(se)$counts[rownames(injectedCounts),idxCol] <- injectedCounts
      binOutliers[outlier_bin,idxCol] <- index[idxRow,idxCol]
      locusOutliers[rownames(locusOutliers) == rowData(se)[outlier_bin,"locus"],idxCol] <- index[idxRow,idxCol]
      
    }
  }
  assay(se, 'binCorruptedCounts', withDimnames=FALSE) <- binOutliers
  assay(se, 'locusCorruptedCounts', withDimnames=FALSE) <- locusOutliers
  metadata(se)$locusCorruptedCounts <- index
  
  return(se)
}


.evalAutoCorrection <- function(se, svd_latent_factors,
                                dimensions, fit_DGE, BPPARAM, ...){
  # Search for corrupted counts by RUV with simultaneous parameter estimation
    se <- .fitRUV_edgeR(se = se, svd_latent_factors = svd_latent_factors,
                dimensions = dimensions,
                fit_DGE = fit_DGE, ...)

  eloss <- .evalAucPRLossSplicing(se)
  print(paste0('Evaluation loss: ', eloss,' for q=',dimensions))
  return(eloss)
}

fitRUVMD <- function(se, dimensions = NULL, padjust = "BH", BPPARAM = bpparam()){
  if(is.null(dimensions)){
    if(is.null(metadata(se)[['optimalEncDim']])){
      dimensions <- round(ncol(se)/4,0)
    }
    else{
      dimensions <- metadata(se)[['optimalEncDim']]
    }
  }
  
  #se <- OffsetSplicing(se)
  
  DGE <- .fitEdgeRDisp(se = se, design = getDesign(se), assays(se)$offsets)
  fit_DGE <- glmFit(y = DGE)
  
  # Principal component analysis on the deviance residuals
  svd_latent_factors <- .svd_deviances(se = se,
                                       fit_DGE = fit_DGE)[,1:dimensions]
  
  
  se <- .fitRUV_edgeR(se = se, svd_latent_factors = svd_latent_factors,
                      dimensions = dimensions, 
                      fit_DGE = fit_DGE)
  
  se <- .adjustPValues(se, method = padjust)
  
  return(se)
  
}

.fitRUV_edgeR <- function(se, svd_latent_factors, dimensions, fit_DGE){
  # Extract the number of latent factors that is estimated to be optimal
  # by prior knowledge or hyperparameter optimisation.
  total_dimensions <- dimensions
  if (total_dimensions == 1){
  	reduced_svd_latent_factors <- svd_latent_factors
  } else{
  	reduced_svd_latent_factors <- svd_latent_factors[,1:total_dimensions]
  }
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
  
  DGE <- .fitEdgeRDisp(se = se, design = getDesign(se), assays(se)$offsets)
  fit_DGE <- glmFit(DGE, design = design)
  
  mu <- fit_DGE$fitted.values
  theta <- 1/fit_DGE$dispersion
  
  assay(se,"mu",withDimnames = FALSE) <- mu
  rowData(se)$theta <- theta
  
  
  #Calculate 2-sided p-values for every count based on previous estimates.
  se <- .calcPValues(se, mu = mu, theta = theta)
  se <- .locusPValues(se)
  
  
  # Store the p-values, estimated means and dispersion to the
  # OutriderDataSet.
  
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
    design <- model.matrix(getDesign(se), data = colData(se))
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








.calcPValues <- function(se, mu, theta){
  left_quantile <- pnbinom(q = counts(se), mu = mu, size = theta)
  right_quantile <- pnbinom(q = counts(se)-1, mu = mu, size = theta)
  PValues <- pmin(left_quantile*2, (1-right_quantile)*2,1)
  assay(se, 'pValue', withDimnames=FALSE) <- PValues
  return(se)
}

library(poolr)



.locusPValues <- function(se){

  pvalues_before_group <- data.frame(assays(se)$pValue,"locus" = rowData(se)$locus)
  
  pvalues_grouped <- pvalues_before_group %>% group_by(locus) %>% 
    summarise(across(colnames(se),function(i) min(i)))
  
  rownames <- pvalues_grouped$locus
  pvalues_grouped[,1] <- NULL
  rownames(pvalues_grouped) <- rownames
  metadata(se)$pValuesLocus <- as.matrix(pvalues_grouped)
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
.evalAucPRLossSplicing <- function(se){

  # Calculate the precision-recall curve of all p-values, using corrupted
  # counts as true positives.
  order <- match(rownames(metadata(se)$locusCorruptedCounts),rownames(metadata(se)$pValuesLocus))
  scores <- -as.vector( metadata(se)$pValuesLocus)
  labels <- as.vector(metadata(se)$locusCorruptedCounts != 0) + 0
  
  if(any(is.na(scores))){
    warning(sum(is.na(scores)), " P-values where NAs.")
    scores[is.na(scores)] <- min(scores, na.rm=TRUE)-1
  }
  pr <- pr.curve(scores, weights.class0=labels)
  return(max(0, pr$auc.integral, na.rm=TRUE))
}


# Adapted from ASpli - Mancini et al.
wrapper_counts <- function(counts = NULL, 
                           features,
                           targets, 
                           minReadLength, 
                           maxISize,
                           libType="SE",
                           strandMode=0,
                           minAnchor = 10,
                           threshold = 5,
                           BPPARAM = bpparam()){
  
  if(is.null(counts)){
    counts <- gbCounts(features = features,
                       targets = targets, 
                       minReadLength = minReadLength, 
                       maxISize = maxISize, 
                       minAnchor = minAnchor,
                       libType=libType,
                       strandMode=strandMode,
                       BPPARAM = BPPARAM)
  }
  
  jcounts <- ASpli:::jCounts(counts= counts, 
                             features = features, 
                             minReadLength = minReadLength, 
                             threshold = 5,
                             minAnchor = minAnchor,
                             libType=libType,
                             strandMode=strandMode)
  
  
  
  
  SEcounts <- SummarizedExperiment::SummarizedExperiment(metadata = 
                                                           list("binCounts" = counts@exon.intron.counts,
                                                                "junctionCounts" = jcounts@junctionsPJU,
                                                                "IRCounts" = jcounts@junctionsPIR),
                                                         colData = list("Names" = rownames(targets)))
  return(SEcounts)
  
}

# Adapted from ASpli - Mancini et al.
FilterBinCounts <- function(SEcounts, 
                            ignoreExternal = FALSE,
                            ignoreIo = TRUE,
                            ignoreI = FALSE,
                            min.count = 10,
                            min.total.count = 15,
                            large.n = 10,
                            min.prop = 0.7){
  SEBinCounts <- metadata(SEcounts)$binCounts
  SEBinCounts = SEBinCounts[ ! ignoreExternal | SEBinCounts$event != "external" ,] 
  SEBinCounts = SEBinCounts[ ! ignoreIo | SEBinCounts$feature != "Io" ,] 
  SEBinCounts = SEBinCounts[ ! ignoreI | SEBinCounts$feature != "I" ,] 
  counts <- SEBinCounts[,colData(SEcounts)$Names]
  
  return(SEcounts)
}

# Adapted from ASpli - Mancini et al.
FilterJunctionCounts <- function(SEcounts, 
                                 ignoreExternal = FALSE,
                                 ignoreIo = TRUE,
                                 ignoreI = FALSE,
                                 min.count = 10,
                                 min.total.count = 15,
                                 large.n = 10,
                                 min.prop = 0.7,
                                 start_J3 = 9,
                                 number_of_samples = dim(SEcounts)[2]){
  
  
  
  SEJunctionCounts <- metadata(SEcounts)$junctionCounts
  
  
  start_J1              <- grep("StartHit", colnames(SEJunctionCounts)) + 1
  start_J2              <- grep("EndHit", colnames(SEJunctionCounts)) + 1
  start_J3              <- start_J3
  end_J3                <- start_J3 + number_of_samples - 1
  
  junctions_of_interest <- rep(TRUE, dim(SEJunctionCounts)[1])
  names(junctions_of_interest) <- rownames(SEJunctionCounts)
  
  # junctions_of_interest <- filterByExpr(SEJunctionCounts[,start_J3:end_J3],
  #                                       min.count = min.count,
  #                                       min.total.count = min.total.count,
  #                                       large.n = large.n,
  #                                       min.prop = min.prop)
  # 
  
  if(sum(junctions_of_interest) == 0) stop("No junctions to analyze! Is minReadLength less than bam read length?")
  
  J1                    <- as.character(SEJunctionCounts$StartHit[rownames(SEJunctionCounts) %in% names(junctions_of_interest)])
  J2                    <- as.character(SEJunctionCounts$EndHit[rownames(SEJunctionCounts) %in% names(junctions_of_interest)])
  J3                    <- names(junctions_of_interest)
  
  if(length(J1) == 0 | length(J2) == 0 | length(J3) == 0) stop("No junctions to analyze! Is minReadLength less than bam read length? Maybe set strongFilter to FALSE?")
  clusters              <- ASpli:::.makeClusters(J1, J2, J3, bStrongFilter = FALSE)
  
  if(clusters$no == 0) stop("No junctions to analyze! Is minReadLength less than bam read length? Maybe set strongFilter to FALSE?")
  
  countData             <- ASpli:::.makeCountDataWithClusters(SEJunctionCounts[names(clusters$membership),start_J3:end_J3], 
                                                              clusters)
  
  
  
  metadata(SEcounts)$junctionCounts <- countData
  return(SEcounts)
}



# Adapted from ASpli - Mancini et al.
gbCounts <- function( features, targets,  minReadLength,  
                      maxISize, minAnchor = 10,
                      libType="SE",
                      strandMode=0,
                      BPPARAM = bpparam()) {
  
  counts <- readCounts( features = features,
                        bam = NULL, 
                        targets = targets, 
                        readLength = minReadLength, 
                        maxISize = maxISize, 
                        minAnchor = minAnchor,
                        libType=libType,
                        strandMode=strandMode,
                        BPPARAM = BPPARAM)
  
  counts@.ASpliVersion = "2" #Marks ASpliCounts object with the ASpli update 2.0.0
  
  return(counts)
}

# Adapted from ASpli - Mancini et al.
readCounts <- function( features, bam, targets, cores = 1,
                        readLength,  
                        maxISize, 
                        minAnchor = 10,
                        libType=libType,
                        strandMode=strandMode,
                        BPPARAM = bpparam()) {
  
  if(!is.null(bam)){
    .Deprecated("gbCounts")
  }
  minReadLength <- readLength
  cores <- 1 #Allways use 1 core.
  
  #Create result object
  counts <- new(Class="ASpliCounts")
  counts@.ASpliVersion = "1" #Last version before 2.0.0 was 1.14.0.
  
  #Generates sample names in case there arent any
  targets <- ASpli:::.generateSamplesNames(targets)
  counts@targets <- ASpli:::.condenseTargetsConditions(targets) #ACH
  group                  <- counts@targets$condition
  counts@condition.order <- levels(factor( group, unique( group ), ordered = TRUE ))
  
  #Minimal anchors
  minAnchor <- if ( ! is.null(minAnchor) ) minAnchor else 10
  minA <- round( minAnchor * minReadLength / 100 )
  ptm <- proc.time()
  if(is.null(bam)) {
    ntargets <- nrow(targets)
  }else{
    ntargets <- 1
  }
  
  # Adapted from ASpli - Mancini et al.
  counts_list <- bplapply(X = c(1:ntargets), BPPARAM = BPPARAM,
                          FUN = function(i){
                            bamToCounts(counts = counts,
                                        features = features,
                                        targets = targets[i,], 
                                        minReadLength = minReadLength, 
                                        maxISize = maxISize, 
                                        minA = minA,
                                        libType=libType,
                                        strandMode=strandMode)})
  
  for (target in c(1:ntargets)){
    
    if(ncol(counts@gene.counts) == 0){
      counts@gene.counts <- counts_list[[target]]$gene.hits
    } else{
      counts@gene.counts <- cbind(counts@gene.counts, 
                                  ASpli:::.extractCountColumns(counts_list[[target]]$gene.hits, targets[target, ]))
      colnames(counts@gene.counts)[ncol(counts@gene.counts)] <- rownames(targets)[target]
    }
    
    
    if(ncol(counts@exon.intron.counts) == 0){
      counts@exon.intron.counts <- counts_list[[target]]$exons.hits
    } else{
      counts@exon.intron.counts <- cbind(counts@exon.intron.counts, 
                                         ASpli:::.extractCountColumns(counts_list[[target]]$exons.hits, targets[target, ]))
      colnames(counts@exon.intron.counts)[ncol(counts@exon.intron.counts)] <- rownames(targets)[target]      
    }
    
    if(ncol(counts@e1i.counts) == 0){
      counts@e1i.counts <- counts_list[[target]]$e1i.hits
    }else{
      counts@e1i.counts <- cbind(counts@e1i.counts, 
                                 ASpli:::.extractCountColumns(counts_list[[target]]$e1i.hits, targets[target, ]))
      colnames(counts@e1i.counts)[ncol(counts@e1i.counts)] <- rownames(targets)[target]         
    }
    
    if(ncol(counts@ie2.counts) == 0){
      counts@ie2.counts <- counts_list[[target]]$ie2.hits
    }else{
      counts@ie2.counts <- cbind(counts@ie2.counts, 
                                 ASpli:::.extractCountColumns(counts_list[[target]]$ie2.hits, targets[target, ]))
      colnames(counts@ie2.counts)[ncol(counts@ie2.counts)] <- rownames(targets)[target]   
    }
    
    if(ncol(counts@junction.counts) == 0){
      counts@junction.counts <- counts_list[[target]]$junction.hits
      junction.hits <- counts_list[[target]]$junction.hits
    }else{
      dt1                    <- data.table(counts@junction.counts, keep.rownames = TRUE)
      dt2                    <- data.table(ASpli:::.extractCountColumns(counts_list[[target]]$junction.hits, 
                                                                        targets[target, ]),
                                           keep.rownames = T)
      dt3                    <- data.frame(merge(dt1, dt2, by="rn", all.x=T, all.y=TRUE))
      junction.hits <- counts_list[[target]]$junction.hits
      for(s in c("junction", "gene", "strand", "multipleHit", "symbol", "gene_coordinates", "bin_spanned", "j_within_bin")){
        dt3[, s]           <- as.character(dt3[, s])
        junction.hits[, s] <- as.character(junction.hits[, s])
      }
      rownames(dt3)          <- dt3[, "rn"]
      dt3                    <- dt3[, -1]
      dt3[dt2$rn, 1:8]       <- ASpli:::.extractDataColumns(junction.hits, targets[target, ])
      counts@junction.counts <- dt3
      counts@junction.counts[is.na(counts@junction.counts)] <- 0
    }
    
    if(length(grep("NA", rownames(counts@junction.counts))) > 0){
      print(target)
      break
    }
    if(length(grep("NA", rownames(junction.hits ))) > 0){
      print(target)
    }
    gc()
    
    
  }
  
  
  for(s in c("junction", "gene", "strand", "multipleHit", "symbol", "gene_coordinates", "bin_spanned", "j_within_bin")){
    counts@junction.counts[, s] <- as.factor(counts@junction.counts[, s])
  }
  colnames(counts@junction.counts)[9:ncol(counts@junction.counts)] <- rownames(targets)
  junctions.order <- sort(rownames(counts@junction.counts))
  junctions.order <- strsplit2(junctions.order, "[.]")
  junctions.order <- GRanges(seqnames=junctions.order[, 1], IRanges(start=as.numeric(junctions.order[, 2]), end=as.numeric(junctions.order[, 3])))
  junctions.order <- sort(junctions.order)
  junctions.order <- paste(junctions.order@seqnames, junctions.order@ranges@start, (junctions.order@ranges@start+junctions.order@ranges@width-1), sep=".")
  counts@junction.counts <- counts@junction.counts[junctions.order, ]
  
  # Create result object
  counts <- ASpli:::rds( counts, targets )
  gc()
  return(counts)
  
  
}

# Adapted from ASpli - Mancini et al.
bamToCounts <- function(  counts = counts,
                          features = features,
                          targets = targets[i,], 
                          minReadLength = minReadLength, 
                          maxISize = maxISize, 
                          minA = minA,
                          libType=libType,
                          strandMode=strandMode){
  
  bam <- ASpli:::loadBAM(targets, cores = NULL,
                         libType=libType, 
                         strandMode=strandMode)
  
  # Count Genes
  gene.hits <- ASpli:::.counterGenes( bam, featuresg( features ))
  counts@gene.counts <- gene.hits
  
  # Count exons 
  bins <- ASpli:::featuresb( features )
  exons.hits <- ASpli:::.counterBin( bam, bins, gene.hits)
  counts@exon.intron.counts <- exons.hits
  
  
  
  # Count introns
  introns <- c( bins[ mcols(bins)$feature == "I" ], 
                bins[ mcols(bins)$feature == "Io"],
                bins[ mcols(bins)$eventJ  == "IR"])
  
  # Count exon1 - intron regions
  e1i <- introns
  start( e1i ) <- start( introns ) - ( minReadLength - minA )
  end( e1i )   <- start( introns ) + ( minReadLength - minA )
  e1i.hits     <- ASpli:::.counterJbin(bam, e1i, gene.hits, minReadLength)
  counts@e1i.counts <- e1i.hits
  
  
  # Count intron - exon2 regions
  ie2 <- introns
  start( ie2 ) <- end( introns ) - ( minReadLength - minA )
  end( ie2 )   <- end( introns ) + ( minReadLength - minA )
  ie2.hits     <- ASpli:::.counterJbin( bam, ie2, gene.hits, minReadLength )
  counts@ie2.counts <- ie2.hits
  
  # Count junctions
  junction.hits    <- ASpli:::.counterJunctions( features, bam, maxISize )
  counts@junction.counts <- junction.hits
  
  
  return(list("gene.hits" = as.data.frame(gene.hits), 
              "exons.hits" = as.data.frame(exons.hits), 
              "e1i.hits" = as.data.frame(e1i.hits),
              "ie2.hits" = as.data.frame(ie2.hits),
              "junction.hits" = as.data.frame(junction.hits)))
  
}


library(tidyr)

OffsetSplicing <- function(se){
  
  offsets <- aggregate(assays(se)$counts,
                       by = list("locus" = as.factor(rowData(se)$locus)),
                       FUN = sum) 
  
  offsets <- merge(list("locus" = rowData(se)$locus),offsets, by = "locus") #%>% mutate("locus" = NULL)
  offsets <- offsets[match(rowData(se)$locus,offsets$locus),] %>% mutate("locus" = NULL)
  
  colnames(offsets) <- colnames(se)
  rownames(offsets) <- rownames(se)
  assays(se)$counts[offsets == 0] <- 1
  offsets[offsets == 0] <- 1
 
  assays(se)$offsets <- as.matrix(offsets)
  
  return(se)
  
}


getCountsSplicing <- function(SECounts, type = c("bin", "junction","IR")){
  
  counts <- as.matrix(metadata(SECounts)[[paste0(type,"Counts")]][,colData(SECounts)$Names])
  
  rownames(counts) <- rownames(metadata(SECounts)[[paste0(type,"Counts")]])
  
  
  se <- SummarizedExperiment::SummarizedExperiment(assays = SimpleList("counts" = counts),
                                                   colData = colData(SECounts),
                                                   rowData = list("locus" = as.factor(metadata(SECounts)[[paste0(type,"Counts")]]$locus)))
  
  
  return(se)
}



