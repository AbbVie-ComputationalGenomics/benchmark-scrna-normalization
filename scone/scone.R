#!/usr/bin/env Rscript
# Priyanka Vijay
# Updated December 6, 2018

cat("\n\n")

############################################################
################### ARGUMENT PARSING #######################
############################################################

suppressPackageStartupMessages(library("argparser"))
p = arg_parser("Run scone on a given dataset to evaluate various normalization methods.

 Package dependencies: 
 argparser
 scone
 RColorBrewer")

p = add_argument(p,"input",
                 help="Path to scone input .RData file")
p = add_argument(p, "data",
                 help="Name of data set (must match prefix of normalized matrix file names) ")
p = add_argument(p, "normdir",
                 help="Path to directory with normalized files")
p = add_argument(p, "--scaling", 
                 short = "-s",
                 default = "fq,tmm,sum",
                 help = "comma separated list of scaling functions")
p = add_argument(p,"--outdir",
                 short = "-o",
                 default = "./",
                 help="Output directory. Will be created if it doesn't exist")
p = add_argument(p, "--subset",
                 flag=T,
                 help="If flag is specified, use 1K hvg subset datasets instead of full gene lists")
p = add_argument(p, "--batch",
                  help="Column name for batch variable. Should be present in colData(sce).")
p = add_argument(p, "--bio",
                  help="Column name for biological variable. Should be present in colData(sce) or 
                        sobj@meta.data if --seurat is specified.")
# p = add_argument(p, "--seurat",
#                  help="Path to RDS file with seurat object") ##untested


args <- parse_args(p)
print(args)

############################################################
################### REQUIRED FUNCTIONS #####################
############################################################
print("Loading required packages and functions ... ")
source("required_functions/scaling_functions.R")
source("required_functions/scone_functions.R")

############################################################
################### MAIN ###################################
############################################################

# Create output directory
dir.create(args$outdir, showWarnings = FALSE)
print(paste0("Created output directory at: ",args$outdir))

# Save arguments used to outdir directory
saveRDS(args, file.path(args$outdir,"arguments.rds"))
print(paste0("saved run arguments to ",file.path(args$outdir,"arguments.rds")))

# set variables
# (data, normdir, and set are required variables for scaling functions)
data = args$data
normdir = args$normdir 

# Check arguments - 1k subset
if(args$subset){
  set = "1k"
}else{
  set = "full"
}

# Get scaling functions
scaling = trimws(unlist(strsplit(args$scaling, ",")))
scaling_map = list(none=identity, # Identity - do nothing
               seurat = SEURAT_FN, # User-defined function
               scnorm = SCNORM_FN, 
               linnorm = LINNORM_FN,
               zinbwave = ZINBWAVE_FN,
               census = CENSUS_FN, #full-length specific
               sum = SUM_FN, # SCONE library wrappers...
               tmm = TMM_FN,
               tpm = TPM_FN,
               uq = UQ_FN,
               fq = FQT_FN,
               deseq = DESEQ_FN)
scalinglist = scaling_map[unique(c("none",scaling))]

# Load input
load(args$input)
print("Input object:")
print(sce)

# Check arguments - Seurat
# if(!is.na(args$seurat)){
#   print("Merging seurat data with existing sce object")
#   sobj = readRDS(args$seurat)
#   print(sobj)
#   df = sobj@meta.data
#   if(!all(rownames(df)==colnames(sce))){
#     stop("Cell names of seurat object don't match sce scone input object")
#   }
#   sce = mutate(sce, df)
# }


# Check arguments - bio and batch
if(!is.na(args$batch)){
  if(!args$batch %in% colnames(colData(sce))){
    stop("Specified batch does not exist in colnames(colData(sce))")
  }
}
if(!is.na(args$bio)){
  if(!args$bio %in% colnames(colData(sce))){
    stop("Specified bio does not exist in colnames(colData(sce))")
  }
}

# Run scone QC (generate QC plots but do not filter)
print("Generating QC plots ... ")
pdf(file.path(args$outdir,"qc_plots.pdf"), width=8,height=6)
scone_input = scone_qc(sce, batchvar = args$batch)
dev.off()

# Initialize scone run
print("Initializing scone experiment ... ")
scone_obj = scone_init(scone_input, batchvar = args$batch, biovar = args$bio)

# Run scone 
print("Running scone experiment ... ")

scone_result_obj = scone_main(scone_obj,
                              scaling = scalinglist,
                              batchvar = args$batch, 
                              biovar = args$bio)
print("Done.")

path = file.path(args$outdir, "scone_results_obj.RDS")
print(paste0("Saving scone results object to: ", path))
saveRDS(scone_result_obj, file=path)


print(scone_result_obj)
