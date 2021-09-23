## v2: remove overall score

library(dplyr)
library(reshape2)
library(ggforce)
library(RColorBrewer)
library(tidyverse)

source("scIBknitTable.R")

heatmap_setup = function(metricsfile, sconescorefile, scalabilityfile, ruv, sim, platform, bio, batch, dataset){
  metricsdf = read.csv(metricsfile, row.names=1)
  sconescoredf = read.csv(sconescorefile, row.names=1)
  scalabilitydf = read.csv(scalabilityfile)
  
  ## subset above files per input filters
  if(!is.null(ruv)){
    metricsdf = subset(metricsdf, RUV_Adjustment == ruv)
    sconescoredf = subset(sconescoredf, RUV_Adjustment == ruv)
  }
  if(!is.null(sim)){
    metricsdf = subset(metricsdf, Simulation == sim)
    sconescoredf = subset(sconescoredf, Simulation == sim)
    # scalabilitydf = subset(scalabilitydf, simulation == sim)
  }
  if(!is.null(platform)){
    metricsdf = subset(metricsdf, Platform == platform)
    sconescoredf = subset(sconescoredf, Platform == platform)
    # scalabilitydf = subset(scalabilitydf, platform == platform)
  }
  if(!is.null(bio)){
    metricsdf = subset(metricsdf, Known_Bio == bio)
    sconescoredf = subset(sconescoredf, Known_Bio == bio)
  }
  if(!is.null(batch)){
    metricsdf = subset(metricsdf, Known_Batch == batch)
    sconescoredf = subset(sconescoredf,  Known_Batch == batch)
  }
  if(!is.null(dataset)){
    metricsdf = subset(metricsdf, Data == dataset)
    sconescoredf = subset(sconescoredf,  Data == dataset)
  }
  
  print(dim(metricsdf))
  print(dim(sconescoredf))
  print(dim(scalabilitydf))
  
  ## Define signatures that need to be interpreted inversely (neg sigs) and need a log scale (scalability)
  ## Changed this to exclude scone metrics since they have already been transformed within scone (-1*Metric)
  # negSig <- c("BATCH_SIL", "EXP_QC_COR", "EXP_UV_COR", "RLE_MED", "RLE_IQR", "RunTime.S", "MaxVMSize_G", "Robust")
  negSig <- c("RunTime.S", "MaxVMSize_G", "Robust")
  logSig <- c("RunTime.S", "MaxVMSize_G")
  
  ## Format metrics df 
  metricsdf = metricsdf[, c("Scaling", "Data", "Scone_Metric", "Z_CollapsedMean")] ## use datasubset normalized zscore instead of the raw score that was being used previously
  
  ## Get mean across datasets for each metric for final plot
  dfplot <- aggregate(metricsdf[, 4], list(metricsdf$Scaling, metricsdf$Scone_Metric), mean, na.rm=TRUE)
  to_plot <- dcast(dfplot, value.var = "x", Group.1 ~ Group.2)
  
  # add scone score and robustness calculation from scorescoredf 
  to_plot$Scone_Avg <- aggregate(sconescoredf[, "Z_CollapsedMean"], list(sconescoredf$Scaling), mean, na.rm=TRUE)$x
  to_plot$Robust <- aggregate(sconescoredf[, "Z_CollapsedSD"], list(sconescoredf$Scaling), mean, na.rm=TRUE)$x
  
  # add scalability results
  scalabilitydf = scalabilitydf[,c("normalization", "RunTime.S", "MaxVMSize_G")]
  dfplot2 <- aggregate(scalabilitydf[, c(2,3)], list(scalabilitydf$normalization), mean, na.rm=TRUE)
  to_plotfinal <- merge(dfplot2, to_plot)
  
  # convert time and memory to log scale
  to_plotfinal[logSig] = log10(to_plotfinal[logSig] + 1)
  
  print(colnames(to_plotfinal))
  # normalize scores to be in range 0 - 1
  # data <- as.data.frame(apply(to_plotfinal[1:13,2:12], MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X))))
  data <- as.data.frame(apply(to_plotfinal[,2:ncol(to_plotfinal)], MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))) ## removed row subsetting, changed col range since with adding "Scone_Avg" there are 13 total instead of 12
  data <- cbind(Method=to_plotfinal[,1], data)
  
  # 1-score for negative signatures
  data[negSig] = 1-data[negSig]
  
  # calculate summary of scalability (row mean)
  data$Resource_Avg = rowMeans(data[,c("RunTime.S", "MaxVMSize_G")])
  
  # calculated a WEIGHTED overall score, account for the fact that robustness calculation is only available for the large datasets 
  if(all(is.na(data$Robust)) & "BIO_SIL" %in% colnames(data)){
    #remove `Robust` from Overall_Score calculation and final plot
    # data$Overall_Score = (95*data$Scone_Avg + 5*data$Resource_Avg)/100
    # rename columns for plotting and sort by overall score
    data <- data[c("Method","Scone_Avg", "Resource_Avg", "BIO_SIL","PAM_SIL","EXP_WV_COR","BATCH_SIL","EXP_QC_COR","EXP_UV_COR","RLE_IQR","RLE_MED", "RunTime.S", "MaxVMSize_G")]
    colnames(data) <- c("Method","Scone_Score", "Scalability", "BIO_SIL","PAM_SIL","EXP_WV_COR","BATCH_SIL","EXP_QC_COR","EXP_UV_COR","RLE_IQR","RLE_MED", "1-log(RunTime)", "1-log(Memory)") 
  }else if(all(is.na(data$Robust)) & !"BIO_SIL" %in% colnames(data)){
    # no robustness and other metrics missing (simulations)
    #remove `Robust` from Overall_Score calculation and final plot
    # data$Overall_Score = (95*data$Scone_Avg + 5*data$Resource_Avg)/100
    # rename columns for plotting and sort by overall score
    data <- data[c("Method","Scone_Avg", "Resource_Avg", "PAM_SIL","BATCH_SIL","EXP_QC_COR","RLE_IQR","RLE_MED", "RunTime.S", "MaxVMSize_G")]
    colnames(data) <- c("Method", "Scone_Score", "Scalability", "PAM_SIL","BATCH_SIL","EXP_QC_COR","RLE_IQR","RLE_MED", "1-log(RunTime)", "1-log(Memory)") 
  }else if(!all(is.na(data$Robust)) & !"BIO_SIL" %in% colnames(data)){
    #robustness included but metrics missing, i don't think this case is present.. 
    # data$Overall_Score = (90*data$Scone_Avg + 5*data$Resource_Avg + 5*data$Robust)/100
    # rename columns for plotting and sort by overall score
    data <- data[c("Method","Scone_Avg", "Resource_Avg", "PAM_SIL","BATCH_SIL","EXP_QC_COR","RLE_IQR","RLE_MED", "RunTime.S", "MaxVMSize_G", "Robust")]
    # colnames(data) <- c("Method","Overall_Score", "Scone_Score", "Resource_Score", "BIO_SIL","PAM_SIL","EXP_WV_COR","1-BATCH_SIL","1-EXP_QC_COR","1-EXP_UV_COR","1-RLE_IQR","1-RLE_MED", "1-log(RunTime)", "1-log(Memory)", "1-Robustness")
    ## changed column names to remove the "1-" from all scone metrics and robustness
    colnames(data) <- c("Method", "Scone_Score", "Scalability", "PAM_SIL","BATCH_SIL","EXP_QC_COR","RLE_IQR","RLE_MED", "1-log(RunTime)", "1-log(Memory)", "Robustness")
  }else{
    # data$Overall_Score = (90*data$Scone_Avg + 5*data$Resource_Avg + 5*data$Robust)/100
    # rename columns for plotting and sort by overall score
    data <- data[c("Method","Scone_Avg", "Resource_Avg", "BIO_SIL","PAM_SIL","EXP_WV_COR","BATCH_SIL","EXP_QC_COR","EXP_UV_COR","RLE_IQR","RLE_MED", "RunTime.S", "MaxVMSize_G", "Robust")]
    # colnames(data) <- c("Method","Overall_Score", "Scone_Score", "Resource_Score", "BIO_SIL","PAM_SIL","EXP_WV_COR","1-BATCH_SIL","1-EXP_QC_COR","1-EXP_UV_COR","1-RLE_IQR","1-RLE_MED", "1-log(RunTime)", "1-log(Memory)", "1-Robustness")
    ## changed column names to remove the "1-" from all scone metrics and robustness
    colnames(data) <- c("Method", "Scone_Score", "Scalability", "BIO_SIL","PAM_SIL","EXP_WV_COR","BATCH_SIL","EXP_QC_COR","EXP_UV_COR","RLE_IQR","RLE_MED", "1-log(RunTime)", "1-log(Memory)", "Robustness") 
  }
  # data <- data[order(data$Overall_Score, decreasing = T), ]
  data <- data[order(data$Scone_Score, decreasing = T), ]
  return(data)
}

heatmap_main = function(ruv = NULL, sim = NULL, bio = NULL, batch = NULL, platform = NULL, dataset = NULL, plottitle="All datasets", plotfile="plots/allData_alluv.png"){
  metricsfile = "datafiles/public_final_individual_normalized_collapsed_scores.csv" ## this is the file with scores for the individual scone metrics (ex. BATCH_SIL), collapsed across dataset subsamplings
  sconescorefile = "datafiles/public_final_aggregated_normalized_collapsed_scores.csv" ## this is the file with the aggregated scone scores, collapsed across dataset subsamplings
  scalabilityfile = "datafiles/benchmarking_completed_validmatrix_public.csv" ## runtime and memory usage stats file
  
  data = heatmap_setup(metricsfile, sconescorefile, scalabilityfile, ruv, sim, platform, bio, batch, dataset) ## function to put it all together for the input slice/dice of data
  
  # plot
  # row_info <- data.frame(id = data$Method, group="All datasets")
  row_info <- data.frame(id = data$Method, group=plottitle)
  
  if("Robustness" %in% colnames(data) & "BIO_SIL" %in% colnames(data)){
    # all possible columns present
    column_info <- data.frame(id = colnames(data),
                              group = c("Text", rep("Overall",2), rep("Biological Accuracy",3), "Batch Removal", rep("Artifact Removal", 2), rep("Differential Reduction", 2), rep("Scalability", 3)), 
                              geom = c("text", rep("bar", 2), rep("circle", 11)),
                              width = c(5, rep(2.5,2), rep(1.5, 11)),
                              overlay = F) 
  }else if(!"Robustness" %in% colnames(data) & "BIO_SIL" %in% colnames(data)){
    # robustness missing, all other metrics present
    column_info <- data.frame(id = colnames(data),
                              group = c("Text", rep("Overall",2), rep("Biological Accuracy",3), "Batch Removal", rep("Artifact Removal", 2), rep("Differential Reduction", 2), rep("Scalability", 2)), 
                              geom = c("text", rep("bar", 2), rep("circle", 10)),
                              width = c(5, rep(2.5,2), rep(1.5, 10)),
                              overlay = F)
  }else if(!"Robustness" %in% colnames(data) & !"BIO_SIL" %in% colnames(data)){
    # robustness missing, other metrics also missing. present: "PAM_SIL","BATCH_SIL","EXP_QC_COR","RLE_IQR","RLE_MED",
    column_info <- data.frame(id = colnames(data),
                              group = c("Text", rep("Overall",2), rep("Biological Accuracy",1), "Batch Removal", rep("Artifact Removal", 1), rep("Differential Reduction", 2), rep("Scalability", 2)), 
                              geom = c("text", rep("bar", 2), rep("circle", 7)),
                              width = c(5, rep(2.5,2), rep(1.5, 7)),
                              overlay = F)
  }else if("Robustness" %in% colnames(data) & !"BIO_SIL" %in% colnames(data)){
    # robustness present, other metrics missing
    column_info <- data.frame(id = colnames(data),
                              group = c("Text", rep("Overall",2), rep("Biological Accuracy",1), "Batch Removal", rep("Artifact Removal", 1), rep("Differential Reduction", 2), rep("Scalability", 3)), 
                              geom = c("text", rep("bar", 2), rep("circle", 8)),
                              width = c(5, rep(2.5,2), rep(1.5, 8)),
                              overlay = F)
  }
  
  palettes <- list("Overall" = "Blues", 
                   "Biological Accuracy" = "Oranges",
                   "Batch Removal" = "Greens",
                   "Artifact Removal" = "Reds",
                   "Differential Reduction" = "Greys",
                   "Scalability" = "Purples")
  
  g <- scIB_knit_table(data = data, column_info = column_info, row_info = row_info, palettes = palettes, usability = F)
  # ggsave("allDatasets.png", width=10,height=6)
  ggsave(plotfile, width=10,height=6)
}

###### MAKE PLOTS########

# all data, all UV adjustments, all bio, all batch
heatmap_main(ruv=NULL, sim=NULL, platform=NULL, bio=NULL, batch=NULL, plottitle = "All datasets", plotfile="plots/allData.png")

# all data, no UV adjustment, no bio, no batch
heatmap_main(ruv="no_uv", sim=NULL, platform=NULL, bio="no_bio", batch="no_batch", plottitle = "All datasets", plotfile="plots/allData_nouv_nobio_nobatch.png")

# real 10x data, no UV adjustment, no bio, no batch
heatmap_main(ruv="no_uv", sim=FALSE, platform="10x", bio="no_bio", batch="no_batch", plottitle = "Real 10x data", plotfile="plots/real10x_nouv_nobio_nobatch.png")

# real 10x data, RUV k=3 adjustment, no bio, no batch
heatmap_main(ruv="ruv_k=3", sim=FALSE, platform="10x", bio="no_bio", batch="no_batch", plottitle = "Real 10x data", plotfile="plots/real10x_ruvk3_nobio_nobatch.png")

# real full length data, no UV adjustment, no bio, no batch
heatmap_main(ruv="no_uv", sim=FALSE, platform="full_length", bio="no_bio", batch="no_batch", plottitle = "Real full-length data", plotfile="plots/realfulllength_nouv_nobio_nobatch.png")

# real full length data, RUV k=3 adjustment, no bio, no batch
heatmap_main(ruv="ruv_k=3", sim=FALSE, platform="full_length", bio="no_bio", batch="no_batch", plottitle = "Real full-length data", plotfile="plots/realfulllength_ruvk3_nobio_nobatch.png")

# simulations, no UV adjustment, no bio, no batch
heatmap_main(ruv="no_uv", sim=TRUE, platform=NULL, bio="no_bio", batch="no_batch", plottitle = "Simulated data", plotfile="plots/sim_nouv_nobio_nobatch.png")

# simulations, qc adjustment, no bio, no batch
heatmap_main(ruv="qc_k=3", sim=TRUE, platform=NULL, bio="no_bio", batch="no_batch", plottitle = "Simulated data", plotfile="plots/sim_qck3_nobio_nobatch.png")

# cellmix, ruv k=3, no batch, no bio
heatmap_main(ruv="ruv_k=3", dataset="cellmix", bio="no_bio", batch="no_batch", plottitle = "cellmix", plotfile="plots/cellmix_ruvk3_nobio_nobatch.png")

# cellmix, ruv k=3, yes batch, no bio
heatmap_main(ruv="ruv_k=3", dataset="cellmix", bio="no_bio", batch="batch", plottitle = "cellmix", plotfile="plots/cellmix_ruvk3_nobio_yesbatch.png")

# cellmix, ruv k=3, yes batch, yes bio
heatmap_main(ruv="ruv_k=3", dataset="cellmix", bio="bio", batch="batch", plottitle = "cellmix", plotfile="plots/cellmix_ruvk3_yesbio_yesbatch.png")

# cellmix, no uv, no batch, no bio
heatmap_main(ruv="no_uv", dataset="cellmix", bio="no_bio", batch="no_batch", plottitle = "cellmix", plotfile="plots/cellmix_nouv_nobio_nobatch.png")

