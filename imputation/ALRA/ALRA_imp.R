args = commandArgs(trailingOnly=TRUE)
library(SingleCellExperiment)

outdir = "~/lowrank_results/"
datadir  = "~/lowrank_input/"
source("~/lowrank_run/ALRA_function.R")

filen = args[1]
filepath = paste0(datadir, args[2], "/", filen)
sc_mat = read.csv(file.path(filepath))
rownames(sc_mat) = sc_mat[,1]
sc_mat[,1] = NULL
sc_mat = as.matrix(sc_mat)

outname = gsub("_full", "", gsub(".csv", "",  filen))

#assess processing time
start_time   <- Sys.time()
A_norm       <- t(sc_mat)
k_choice     <- choose_k(A_norm)
out_imp      <- alra(A_norm, k = k_choice$k)
end_time     <- Sys.time()
time_proc    <- round(difftime(end_time, start_time, units = "mins"),3)
print(paste("Time used for low-rank imputation of ", args[2], "normalized dataset: ", outname, " is ",
               round(difftime(end_time, start_time, units = "mins"),3), " min",  sep = ""))
sc_imp = t(out_imp[[3]])
cap  = data.frame(norm_method = args[2], dataset = outname, proc_time_min = time_proc)

#write results
write.csv(sc_imp, file = paste(outdir, args[2], "/", Sys.Date(),"_", outname, "_lowrank_imp.csv", sep = "") )
write.csv(cap, file = paste(outdir, args[2], "/metrics/", Sys.Date(),"_", outname, "_lowrank_impmetrics.csv", sep = ""), row.names = FALSE)

