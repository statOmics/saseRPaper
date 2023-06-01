suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(satuRn))

## run edgeR_offset
run_DESeq2_offset <- function(L,countData,tx2gene) {
  message("DESeq2_offset")
  session_info <- sessionInfo()
  timing <- system.time({
    
    tx2gene <- tx2gene_init
    quantsf_counts <- countData_init
    removed <- "none"
    
    ## select cells
    set.seed(123)
    selected_cells <- sample(colnames(quantsf_counts), 2*i)
    
    quantsf_counts <- quantsf_counts[,selected_cells]
    
    ## filter TXs with zero rowSum after cell selection
    quantsf_counts <- quantsf_counts[which(rowSums(quantsf_counts)!=0),]
    tx2gene <- tx2gene[tx2gene$TXNAME %in% rownames(quantsf_counts),]
    
    if(dim(quantsf_counts)[1] > j){
      
      ## select transcripts
      k <- 1
      randomized <- sample(unique(tx2gene$GENEID))
      TXnames <- c()
      
      while (length(TXnames) < j){
        TXnames <- c(TXnames, tx2gene$TXNAME[tx2gene$GENEID == randomized[k]])
        k <- k + 1
      }
      
      quantsf_counts <- quantsf_counts[TXnames,]
      quantsf_counts <- as.data.frame(quantsf_counts)
      
      sampleData <- as.data.frame(matrix(nrow = ncol(quantsf_counts), ncol=2))
      colnames(sampleData) <- c("sampleId", "group")
      sampleData$sampleId <- colnames(quantsf_counts)
      sampleData$group <- as.factor(sample(rep(c("A","B"),each=i)))
      
      ## check if there are cells without counts due to filtering/subsampling
      ## if so (unlikely); remove them, but store the fact they are removed as it may impact speed
      if(length(which(colSums(quantsf_counts) == 0)) > 0){
        removed <- length(which(colSums(quantsf_counts) == 0))
        sampleData <- sampleData[-which(colSums(quantsf_counts) == 0),]
        quantsf_counts <- quantsf_counts[,-which(colSums(quantsf_counts) == 0)]
      }
      
      colnames(tx2gene)[1:2] <- c("isoform_id", "gene_id")
      tx2gene <- tx2gene[which(tx2gene$isoform_id %in% TXnames),]
      tx2gene <- tx2gene[match(rownames(quantsf_counts), tx2gene$isoform_id),]
      
      # get totalCounts
      totalCounts <- satuRn:::getTotalCount(countData = quantsf_counts, 
                                            tx2gene = tx2gene)
      
      # handle zero offsets
      quantsf_counts[quantsf_counts==0] <- 1
      totalCounts[totalCounts==0] <- 1
      
      # Run DESeq2 analysis with offset
      dds <- DESeqDataSetFromMatrix(countData = round(quantsf_counts),
                                    colData = sampleData,
                                    design= ~ group)
      normalizationFactors(dds) <- totalCounts
      dds <- estimateDispersions(dds)
      dds <- nbinomWaldTest(dds)
      
      res <- results(dds, name="group_B_vs_A")
      
      finalResult <- nrow(res[which(res$padj < 0.05),])
      mem <- gc()
    } else {
      ## make sure to have an output for result and info
      quantsf_counts <- matrix(data=NA, nrow=1,ncol=1)
      finalResult <- "More TXs asked than there were present"
      mem <- gc()
    }
  })
  
  output_DESeq2_offset <- list(session_info = session_info,
                              info = list(size = i, 
                                          TXs = nrow(quantsf_counts)),
                              timing = timing,
                              removed = removed,
                              memory = mem,
                              result = finalResult)
  
  assign("result", output_DESeq2_offset, envir=globalenv())
  rm(list = ls())
  gc()
}
