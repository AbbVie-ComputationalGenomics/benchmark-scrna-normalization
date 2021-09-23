args = commandArgs(trailingOnly=TRUE)
library(SingleCellExperiment)
library(scImpute)

#paths
scratchdir = "/mnt/nvme/scratch/wachtax/"
outdir = "/mnt/nvme/users/wachtax/scImpute_results/"
datadir  = "/mnt/nvme/users/wachtax/scImpute_input/"
metricsdir = "/mnt/nvme/users/wachtax/scImpute_results/metrics/"
#read in cluster numbers
nclustall = read.csv("/mnt/nvme/users/wachtax/Seurat_clusters/2019-08-08_Seuratres0.8_clustno.csv", sep = ",")
#set imputation seed no
seed_no = 190806

#load data
filen = args[1]
filepath = paste0(datadir, filen, ".RData")
load(file.path(filepath))

#create scratch output directory
outname = gsub("_filtered", "", gsub("_sconeinput", "", gsub("_gl", "", gsub("2019-06-14_", "", gsub("2019-01-17_", "", gsub("2019-01-16_", "", gsub(".RData", "", filen)))))))
dir.create(paste0(scratchdir, "scRNAseq_methods/", outname, "/"))

#input cluster number
nclus = nclustall[which(nclustall$dataset == gsub(".RData", "", filen)) ,"clust_no"]
print(paste0("Dataset: ", outname, ", Input cluster number: ", nclus))
set.seed(seed_no)

#write input matrices
if(exists("sce_filt"))
{sce = sce_filt
rm(sce_filt)}
write.csv(as.matrix(counts(sce)), file = paste(outdir, "matrices/", outname , ".csv", sep = ""))
rm(sce)

#track time
start_time   <- Sys.time()
scimpute(count_path = paste0(outdir, "matrices/", outname, ".csv"),
			 infile = "csv", outfile = "csv",
			 out_dir = paste0(scratchdir, "scRNAseq_methods/", outname, "/"), Kcluster = nclus,labeled = FALSE, drop_thre = 0.5, ncores = 5)
end_time     <- Sys.time()
time_proc    <- round(difftime(end_time, start_time, units = "mins"),3)
print(paste("scImpute imputation of ", outname, " dataset, input cluster number: ", nclus, ", time used: ", 
               round(difftime(end_time, start_time, units = "mins"),3), " min",  sep = ""))

cap  = data.frame(dataset = outname, input_nclust = nclus, seed = seed_no, proc_time_min = time_proc)
print("\n")

#write metrics and results	
write.csv(cap, file = paste(metricsdir, Sys.Date(),"_", outname, "_scImputed_seed", seed_no,".csv", sep = ""), row.names = FALSE)
file.copy(from=paste0(scratchdir, "scRNAseq_methods/", outname, "/scimpute_count.csv"), to=paste0(outdir, Sys.Date(),"_", outname, "_scImputed.csv") )


