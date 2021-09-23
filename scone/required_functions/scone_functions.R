suppressPackageStartupMessages(library(scone))
suppressPackageStartupMessages(library(RColorBrewer))
data(housekeeping) #housekeeping genes for FNR curves in scone
mmhousekeeping = read.csv("required_functions/mm_housekeeping.csv", header=F)
housekeeping = rbind(housekeeping, mmhousekeeping)

## main scone functions

scone_qc = function(sce, batchvar=NA){
  qcvars = c("total_features_by_counts","log10_total_features_by_counts","total_counts",
             "log10_total_counts","pct_counts_in_top_50_features",
             "pct_counts_in_top_100_features","pct_counts_in_top_200_features",
             "pct_counts_in_top_500_features","pct_mito","pct_ribo",
             "total_features","log10_total_features",    
             "pct_counts_top_50_features","pct_counts_top_100_features",
             "pct_counts_top_200_features","pct_counts_top_500_features",
             "nGene","nUMI")
  metadata(sce)$which_qc = intersect(colnames(colData(sce)),qcvars)
  print(metadata(sce)$which_qc)
  
  cc <- c(brewer.pal(9, "Set1"))
  
  # One batch per Biological Condition
  if(is.na(batchvar)){
    batch = factor(rep(1,nrow(sce)))
  }else{
    batchcol = which(colnames(colData(sce))==batchvar)
    unique(colData(sce)[,batchcol])
    batch = factor(colData(sce)[,batchcol])
  }
  
  # Alignment Quality Metrics
  if(length(metadata(sce)$which_qc)>0){
    qc = colData(sce)[,metadata(sce)$which_qc]
    colsums = colSums(as.matrix(qc))
    if(any(colsums==0)){
      qc = qc[,-which(colsums==0)] 
    }
    # Barplot of total read number
    if("total_counts" %in% names(qc)){
      nreads = qc$total_counts
      o = order(nreads)[order(batch[order(nreads)])] # Order by batch, then value
      barplot(nreads[o], col=cc[batch][o],
              border=cc[batch][o], main="Total number of reads")
      legend("topright", legend=levels(batch), fill=cc, cex=0.4)
    }
    
    ## ----- PCA of QC matrix -----
    qpc = prcomp(qc,center = TRUE,scale. = TRUE)
    barplot((qpc$sdev^2)/sum(qpc$sdev^2), border="gray", 
            xlab="PC", ylab="Proportion of Variance", main="Quality PCA")
    
    
    ## ----qpc_view--------------------------------------------------------------
    
    # Barplot of PC1 of the QC matrix
    qc1 = as.vector(qpc$x[,1])
    o = order(qc1)[order(batch[order(qc1)])]
    
    barplot(qc1[o], col=cc[batch][o], 
            border=cc[batch][o], main="Quality PC1")
    legend("bottomright", legend=levels(batch), 
           fill=cc, cex=0.8)
  }else{
    "No QC columns found"
  }
  ## ----fnr_fit---------------------------------------------------------------
  
  # Extract Housekeeping Genes
  hk = intersect(housekeeping$V1,rownames(assay(sce)))
  if(length(hk)>0){
    # Mean log10(x+1) expression
    mu_obs = rowMeans(log10(assay(sce)[hk,]+1))
    
    # Assumed False Negatives
    drop_outs = assay(sce)[hk,] == 0
    
    # Logistic Regression Model of Failure
    ref.glms = list()
    for (si in 1:dim(drop_outs)[2]){
      fit = glm(cbind(drop_outs[,si],1 - drop_outs[,si]) ~ mu_obs,
                family=binomial(logit))
      ref.glms[[si]] = fit$coefficients
    }
    
    ## ----fnr_vis,fig.width=8,fig.height=4,out.width="800px",out.height="400px"----
    
    par(mfrow=c(1,2))
    
    # Plot Failure Curves and Calculate AUC
    plot(NULL, main = "False Negative Rate Curves",
         ylim = c(0,1),xlim = c(0,6),
         ylab = "Failure Probability", xlab = "Mean log10 Expression")
    x = (0:60)/10
    AUC = NULL
    for(si in 1:ncol(assay(sce))){
      y = 1/(exp(-ref.glms[[si]][1] - ref.glms[[si]][2] * x) + 1)
      AUC[si] = sum(y)/10
      lines(x, 1/(exp(-ref.glms[[si]][1] - ref.glms[[si]][2] * x) + 1),
            type = 'l', lwd = 2, col = cc[batch][si])
    }
    
    # Barplot of FNR AUC
    o = order(AUC)[order(batch[order(AUC)])]
    
    barplot(AUC[o], col=cc[batch][o], border=cc[batch][o], main="FNR AUC")
    legend("topright", legend=levels(batch), fill=cc, cex=0.4)
    
  }
  return(sce)
  
}

scone_init = function(goodDat, biovar=NA, batchvar=NA, poscongenes=NA){
  
  ## ----scone_init------------------------------------------------------------
  
  ## Expression Data (Required); if n_genes is not NA, take only top n_genes most variable genes
  expr = assay(goodDat)
  
  ## Biological Origin - Variation to be preserved (Optional)
  if(is.na(biovar)){
    bio = NULL
  }else{
    bio = factor(colData(goodDat)[,biovar])
  }

  ## Batch
  if(is.na(batchvar)){
    batch = NULL
  }else{
    batch = factor(colData(goodDat)[,batchvar])
  }
  
  # Processed Alignment Metrics - Variation to be removed (Optional)
  qc = colData(goodDat)[,metadata(goodDat)$which_qc]
  print(dim(qc))
  ppq = scale(qc[,apply(qc,2,sd) > 0],center = TRUE,scale = TRUE)
  print(dim(ppq))
  
  # Positive Control Genes - Prior knowledge of DE (Optional)
  if(is.na(poscongenes)){
    poscon = NULL
  }else{
    poscon = intersect(rownames(expr),poscongenes)
  }
  
  # Negative Control Genes - Uniformly expressed transcripts (Optional)
  negcon = intersect(rownames(expr),housekeeping$V1)
  if(length(negcon) == 0){
    negcon = NULL
  }else{
    negcon = rownames(expr) %in% negcon
  }
  
  # Creating a SconeExperiment Object
  if(is.null(batch) & is.null(bio) & is.null(poscon)){
    print("Initializing scone experiment with: no bio, no batch, no poscon")
    scone <- SconeExperiment(expr,
                             qc=ppq, 
                             negcon_ruv = negcon) #,poscon = rownames(expr) %in% poscon # bio = bio, batch = batch,
  }else if(is.null(bio) & is.null(poscon)){
    print("Initializing scone experiment with: no bio, yes batch, no poscon")
    print(length(batch))
    scone <- SconeExperiment(expr,
                             qc=ppq, batch = batch,
                             negcon_ruv = negcon)
  }else if(is.null(poscon)){
    print("Initializing scone experiment with: yes bio, yes batch, no poscon")
    scone <- SconeExperiment(expr,
                             qc=ppq, bio = bio, batch = batch,
                             negcon_ruv = negcon)
  }else{
    print("Initializing scone experiment with: yes bio, yes batch, yes poscon")
    #not tested
    scone <- SconeExperiment(expr,
                             qc=ppq, bio = bio, batch = batch,
                             poscon = rownames(expr) %in% poscon,
                             negcon_ruv = negcon)
  }
  
  
  return(scone)

}

scone_main = function(scone, scaling, batchvar=NA, biovar=NA){

  BiocParallel::register(
    BiocParallel::SerialParam()
  ) # Register BiocParallel Serial Execution
  
  if(length(scone@which_negconruv) == 0 & is.na(batchvar) & is.na(biovar)){
    scone <- scone(scone,
                   scaling=scaling,
                   k_qc=3, k_ruv = 0,
                   adjust_bio="no",
                   adjust_batch="no",
                   run=TRUE,
                   eval_kclust = 2:6,stratified_pam = FALSE,
                   return_norm = "no",
                   zero = "postadjust")
  }else if(length(scone@which_negconruv) == 0 & is.na(biovar)){
    scone <- scone(scone,
                   scaling=scaling,
                   k_qc=3, k_ruv = 0,
                   adjust_bio="no",
                   adjust_batch="yes",
                   run=TRUE,
                   eval_kclust = 2:6,stratified_pam = TRUE,
                   return_norm = "no",
                   zero = "postadjust")
  }else if(length(scone@which_negconruv) == 0 & is.na(batchvar) & !is.na(biovar)){
    scone <- scone(scone,
                   scaling=scaling,
                   k_qc=3, k_ruv = 0,
                   adjust_bio="yes",
                   adjust_batch="no",
                   run=TRUE,
                   eval_kclust = 2:6,stratified_pam = FALSE,
                   return_norm = "no",
                   zero = "postadjust")
  }else if(length(scone@which_negconruv) == 0){
    scone <- scone(scone,
                   scaling=scaling,
                   k_qc=3, k_ruv = 0,
                   adjust_bio="yes",
                   adjust_batch="yes",
                   run=TRUE,
                   eval_kclust = 2:6,stratified_pam = TRUE,
                   return_norm = "no",
                   zero = "postadjust")
  }else if(is.na(batchvar) & is.na(biovar)){
    scone <- scone(scone,
                   scaling=scaling,
                   k_qc=3, k_ruv = 3,
                   adjust_bio="no",
                   adjust_batch="no",
                   run=TRUE,
                   eval_kclust = 2:6,stratified_pam = FALSE,
                   return_norm = "no",
                   zero = "postadjust")
  }else if(is.na(biovar)){
    scone <- scone(scone,
                   scaling=scaling,
                   k_qc=3, k_ruv = 3,
                   adjust_bio="no",
                   adjust_batch="yes",
                   run=TRUE,
                   eval_kclust = 2:6,stratified_pam = TRUE,
                   return_norm = "no",
                   zero = "postadjust")
  }else if(is.na(batchvar) & !is.na(biovar)){
    scone <- scone(scone,
                   scaling=scaling,
                   k_qc=3, k_ruv = 3,
                   adjust_bio="yes",
                   adjust_batch="no",
                   run=TRUE,
                   eval_kclust = 2:6,stratified_pam = FALSE,
                   return_norm = "no",
                   zero = "postadjust")
  }else{
    scone <- scone(scone,
                   scaling=scaling,
                   k_qc=3, k_ruv = 3,
                   adjust_bio="yes",
                   adjust_batch="yes",
                   run=TRUE,
                   eval_kclust = 2:6,stratified_pam = TRUE,
                   return_norm = "no",
                   zero = "postadjust")
  }

  print(get_params(scone))
  
  print(apply(get_params(scone),2,unique))
  return(scone)
}






