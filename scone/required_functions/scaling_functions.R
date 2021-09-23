##required variables: data (aml, appps1. etc.), set (full or 1k), normdir

CPM_FN = function(ei){
  cat("\nCPM\n")
  file=paste0(normdir,"/cpm_normalization_",set,"/",data,"_cpm_",set,".csv")
  if(file.exists(file)){
    norm = as.matrix(read.csv(file, row.names=1))
    print(dim(norm))
    print(norm[1:3,1:3])
    return(norm)
  }else{
    return(NA)
  }
}

ZINBWAVE_FN <- function(x){
  cat("\nZINBWAVE\n")
  file=paste0(normdir,"/zinbwave_normalization_",set,"/",data,"_zinbwave_",set,".csv")
	if(file.exists(file)){
	  norm = as.matrix(read.csv(file, row.names=1))
	  print(dim(norm))
	  print(norm[1:3,1:3])
	  return(norm)
	}else{
	  return(NA)
	}
}

SCNORM_FN <- function(ei){
  cat("\nSC NORM\n")
  file=paste0(normdir,"/scnorm_normalization_",set,"/",data,"_scnorm_",set,".csv")
  if(file.exists(file)){
    norm = as.matrix(read.csv(file, row.names=1))
    print(dim(norm))
    print(norm[1:3,1:3])
    return(norm)
  }else{
    return(NA)
  }
}

SEURAT_FN = function(ei){
  cat("\nSEURAT\n")
  file=paste0(normdir,"/seurat_normalization_",set,"/",data,"_seurat_",set,".csv")
  if(file.exists(file)){
    norm = as.matrix(read.csv(file, row.names=1))
    print(dim(norm))
    print(norm[1:3,1:3])
    return(norm)
  }else{
    return(NA)
  }
}


LINNORM_FN <- function(x){
  cat("\nLINNORM\n")
  file=paste0(normdir,"/linnorm_normalization_",set,"/",data,"_linnorm_",set,".csv")
  if(file.exists(file)){
    norm = as.matrix(read.csv(file, row.names=1))
    print(dim(norm))
    print(norm[1:3,1:3])
    return(norm)
  }else{
    return(NA)
  }
}


#full-length specific:
TPM_FN = function(ei){
  cat("\nTPM\n")
  file=paste0(normdir,"/tpm_normalization_",set,"/",data,"_tpm_",set,".csv")
  if(file.exists(file)){
    norm = as.matrix(read.csv(file, row.names=1))
    print(dim(norm))
    print(norm[1:3,1:3])
    return(norm)
  }else{
    return(NA)
  }
}

CENSUS_FN = function(ei){
  cat("\nCENSUS\n")
  file=paste0(normdir,"/census_normalization_",set,"/",data,"_census_",set,".csv")
  if(file.exists(file)){
    norm = as.matrix(read.csv(file, row.names=1))
    print(dim(norm))
    print(norm[1:3,1:3])
    return(norm)
  }else{
    return(NA)
  }
}