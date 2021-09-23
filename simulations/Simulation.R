###########################################################################################
###  Jen Mollon, Astrid Wachter
###  Create simulated data sets from parameters estimated from all real data sets 
###########################################################################################

library(SingleCellExperiment)
library(Matrix)
library(splatter)
library(scater)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(Rtsne)
library(magrittr)

options(stringsAsFactors = FALSE)
mySeed = 190806
set.seed(mySeed)

#############################################################################################
##1 group simulations

#read in estimated parameters
pars <- read.csv2(file = paste0(path, "Estimated_parameters_from_real_datasets.csv"))
dataSets <- names(pars)[-c(1,2)]
print(pars)
pars <- pars[-1,]

for(k in 1: length(datasets))
{
  path     = "~/sim_datasets_AW/"
  inPath   = "~/scone_input_sce/"
  outpath  = paste0(path,"simData/")
  load(file = paste(inPath, datasets[k], sep = ""))
  
  #subsample large data set
  if(datasets[k] == "HCA_BM_sconeinput.RData" | datasets[k] == "HCA_CB_sconeinput.RData" | 
     datasets[k] == "pbmc_sconeinput.RData")
  { set.seed(mySeed)
    sce <- sce[ ,sample(dim(assay(sce))[2], 10000)]}
  
  #group/batch probabilities
  groupProbList = list(c(1))
  batchCellList = list(c(dim(sce)[2]), c(round(dim(sce)[2]/3),round(dim(sce)[2]/3),floor(dim(sce)[2]/3)))   
  scesim = SingleCellExperiment(assays = list(counts = sce_counts))
  
  #simulate data sets without dropout
  for(i in batchCellList) {
    
    paramListCurrent <- split(as.numeric(pars[, i]), pars[, 1]) %>% .[names(.) != "nCells"]
    paramListCurrent <- paramListCurrent[ !is.na(paramListCurrent)]
    
    localParams <- newSplatParams()
    localParams <- setParams(localParams, paramListCurrent)
    localParams <- setParams(params, update = list(nGenes = 30000, batchCells = i, group.prob=1, seed = mySeed))
    
    simRes <- splatSimulate(localParams, method = "groups")
    simRes <- normalize(simRes)
    simRes <- runTSNE(simRes, perplexity = 20, rand_seed = mySeed)
    
    save(simRes,  file = paste(outpath, Sys.Date(), "_", gsub("_sconeinput.RData", "", datasets[k]), 
                               "_scdatSim3_", length(i),"bat1groups.RData", sep = ""))
    
  }
}


#############################################################################################
##10 group simulations

path    = "~/parameter_estimation/"
outPath = "~/simsRealData/"
setwd(path)

#read in estimated parameters
pars <- read.csv2(file = paste0(path, "Estimated_parameters_from_real_datasets.csv"))
dataSets <- names(pars)[-c(1,2)]
print(pars)
pars <- pars[-1,]

#group/batch probabilities
groupProbs <- c(0.40, 0.17, 0.11, 0.08, 0.06, 0.05, 0.04, 0.03, 0.03)
batchCells <- list(c(10000), c(3333,3333,3334))
deScale <- c(0.6, 0.6, 0.7, 0.7, 0.7, 0.9, 0.9, 1.1, 1.2)

#simulate data sets 
for(i in dataSets) {
  for(j in batchCells) {

    paramListCurrent <- split(as.numeric(pars[, i]), pars[, 1]) %>% .[names(.) != "nCells"]
    paramListCurrent <- paramListCurrent[ !is.na(paramListCurrent)]
    
    localParams <- newSplatParams()
    localParams <- setParams(localParams, paramListCurrent)
    localParams <- setParams(localParams, update = list(nGenes = 30000, batchCells = j, group.prob=groupProbs, de.facScale = deScale, seed = mySeed))
    
    simRes <- splatSimulate(localParams, method = "groups")
    simRes <- normalize(simRes)
    simRes <- runTSNE(simRes, perplexity = 20, rand_seed = mySeed)
    
    save(simRes,  file = paste(outPath, "scdatSim_",i, "_",length(j),"bat_", "10groups.RData", sep = ""))

  }

}



