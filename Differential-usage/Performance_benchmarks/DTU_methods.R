# Load all 5 DTU methods

satuRn_DTU <- function(countData, tx2gene, sampleData, quiet = FALSE){
    
    colnames(tx2gene)[1:2] <- c("isoform_id", "gene_id")
    tx2gene <- tx2gene[match(rownames(countData), tx2gene$isoform_id),]
    condition <- factor(sampleData$condition)
    rownames(sampleData) <- colnames(countData)
    design <- model.matrix(~0+condition)
    
    sumExp <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=countData), colData = sampleData, rowData = tx2gene)

    sumExp <- satuRn::fitDTU(object = sumExp,
                            formula = ~0+condition,
                            parallel = FALSE,
                            BPPARAM = BiocParallel::bpparam(),
                            verbose = TRUE
    )

    L <- matrix(0, ncol = 1, nrow = ncol(design))
    rownames(L) <- colnames(design)
    colnames(L) <- c("AvsB")
    L[1:2,] <- c(-1,1)
    
    sumExp <- satuRn::testDTU(object = sumExp, 
                              contrasts = L)
    
    output <- data.frame(rownames(countData))
    rownames(output) <- rownames(countData)
    colnames(output) <- "TXNAME"
    
    output$GENEID <- tx2gene[match(output$TXNAME,tx2gene$TXNAME),"GENEID"]
    
    result <- rowData(sumExp)[["fitDTUResult_AvsB"]]
    
    output$tvalue <- result$t
    output$p_value_raw <- result$pval
    output$FDR_raw <- result$regular_FDR
    output$p_value <- result$empirical_pval
    output$FDR <- result$empirical_FDR
    
    localRes <- output
    
    return(localRes)
}

## edgeR_diffsplice
edgeR_diffsplice_DTU <- function(countData, tx2gene, sampleData) {
    
    y <- DGEList(counts = countData, group=as.data.frame(sampleData)$condition)
    y <- calcNormFactors(y)
    
    design <- model.matrix(~condition, data=sampleData)
    
    y <- estimateDisp(y, design = design)
    fit <- glmFit(y, design = design) ## glmQLfit --> very poor results
    
    tx2gene <- tx2gene[match(rownames(fit$counts), tx2gene$TXNAME),]
    
    dtuEdgeR <- edgeR::diffSpliceDGE(
      glmfit = fit,
      exonid = tx2gene$TXNAME,
      geneid = tx2gene$GENEID,
      coef = "conditionb",
      verbose = FALSE
    )
    
    localRes <- topSpliceDGE(
      dtuEdgeR,
      test = 'exon',
      number = Inf
    )
    
    rownames(localRes) <- NULL
    colnames(localRes)[which(colnames(localRes) == 'GeneID')] <- 'GENEID'
    colnames(localRes)[which(colnames(localRes) == 'ExonID')] <- 'TXNAME'
    
    localRes <- localRes[,c(
      which( colnames(localRes) == 'TXNAME'),
      which( colnames(localRes) == 'GENEID'),
      which( ! colnames(localRes) %in% c('TXNAME','GENEID'))
    )]
    
    colnames(localRes)[which(colnames(localRes) == 'P.Value')] <- 'p_value'
    
    return(localRes)
}

### DEXSeq
DEXSeq_DTU <- function(countData, tx2gene, sampleData) {
    # get tx2gene in better format
    geneForEachTx <- tx2gene$GENEID[match(rownames(countData),tx2gene$TXNAME)]
    
    countData <- round(countData)
    
    dxd <- DEXSeqDataSet(countData = countData,
                         sampleData = sampleData,
                         design = ~ sample + exon + condition:exon,
                         featureID = rownames(countData),
                         groupID = as.character(geneForEachTx))
    
    dxd <- estimateSizeFactors(dxd)
    ## note that estimating dispersions takes some time!
    print("estimating the dispersions")
    dxd <- estimateDispersions(dxd)
    
    print("testing for DEU")
    dxd <- testForDEU(dxd, reducedModel=~ sample + exon)
    
    print("Getting the results")
    dxr <- DEXSeqResults(dxd)
    
    ### Massage result
    localRes <- as.data.frame(dxr)
    localRes <- localRes[,c('featureID','groupID','exonBaseMean','dispersion','stat','pvalue','padj')]
    colnames(localRes)[1:2] <- c('TXNAME','GENEID')
    rownames(localRes) <- NULL
    
    colnames(localRes)[which(colnames(localRes) == 'pvalue')] <- 'p_value'
    colnames(localRes)[which(colnames(localRes) == 'padj')] <- 'FDR'
    
    return(localRes)
}

edgeR_offset_DTU <- function(countData, 
                             tx2gene, 
                             sampleData,
                             quiet = FALSE){
  
  colnames(tx2gene)[1:2] <- c("isoform_id", "gene_id")
  tx2gene <- tx2gene[match(rownames(countData), tx2gene$isoform_id),]
  condition <- factor(sampleData$condition)
  design <- model.matrix(~condition)
  
  # get totalCounts
  totalCounts <- satuRn:::getTotalCount(countData = countData, 
                                        tx2gene = tx2gene)
  
  # Run edgeR analysis with offset
  y <- DGEList(counts = countData)
  
  countData[totalCounts==0] <- 1
  y <- DGEList(counts = countData)
  totalCounts[totalCounts==0] <- 1
  y$offset <- log(totalCounts)

  yEst <- estimateDisp(y, design = design)
  yFit <- glmFit(yEst, design = design)
  
  yTest <- glmLRT(yFit, coef = 2)
  topRes <- topTags(yTest, n = nrow(yTest$table))
  
  output <- topRes$table
  output$TXNAME <- rownames(output)
  output$GENEID <- tx2gene[match(output$TXNAME,tx2gene$TXNAME),"GENEID"]
  colnames(output)[which(colnames(output) == 'PValue')] <- 'p_value'
  
  localRes <- output
  
  return(localRes)
}

DESeq2_offset_DTU <- function(countData, 
                              tx2gene, 
                              sampleData,
                              quiet = FALSE){
  
  colnames(tx2gene)[1:2] <- c("isoform_id", "gene_id")
  tx2gene <- tx2gene[match(rownames(countData), tx2gene$isoform_id),]
  sampleData$condition <- as.factor(sampleData$condition)
  
  # DESeq2 needs integers. I will do the rounding before computing totalCounts
  countData <- round(countData)
  
  # get totalCounts
  totalCounts <- satuRn:::getTotalCount(countData = countData, 
                                        tx2gene = tx2gene)
  
  # Run DESeq2 analysis with offset
  if("sample_id" %in% colnames(sampleData)){
    rownames(sampleData) <- sampleData$sample_id
  }
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = sampleData,
                                design= ~ condition)
  
  countData[totalCounts==0] <- 1
  dds <- DESeqDataSetFromMatrix(countData = countData + 1,
                                colData = sampleData,
                                design= ~ condition)
  totalCounts[totalCounts==0] <- 1
  normalizationFactors(dds) <- totalCounts  

  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  
  output <- results(dds, name="condition_b_vs_a")

  output$TXNAME <- rownames(output)
  output$GENEID <- tx2gene[match(output$TXNAME,tx2gene$TXNAME),"GENEID"]
  colnames(output)[which(colnames(output) == 'pvalue')] <- 'p_value'
  colnames(output)[which(colnames(output) == 'padj')] <- 'FDR'
  
  localRes <- output
  
  return(localRes)
}
