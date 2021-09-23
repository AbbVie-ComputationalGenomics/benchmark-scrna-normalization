#!/usr/bin/env Rscript

cat("\n--------------------------SEURAT--------------------------------------\n\n")
print(Sys.time())
cat("\n")

# required args: input, output
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=2){
  stop("\nRequired inputs: [input sce .RData object] [name of output file]")
}
args = commandArgs(trailingOnly=TRUE)
infile = args[1] ## should be sce object
outfile = args[2]

cat("\nLoading required packages: Seurat, SingleCellExperiment ... \n")
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(library(SingleCellExperiment))

cat("\nReading input scone sce object ... \n")
print(infile)
load(infile)
print(dim(sce))

cat("\nRunning Seurat NormalizeData, FindVariableFeatures, and ScaleData functions ... \n")
sobj = CreateSeuratObject(counts(sce))
sobj = NormalizeData(object = sobj, normalization.method = "LogNormalize",
                     scale.factor = 10000)
sobj = FindVariableFeatures(object = sobj)
sobj = ScaleData(object = sobj)
mat = as.matrix(GetAssayData(sobj))

cat("\nWriting output to ",outfile," ... \n")
write.csv(mat, file=outfile)

cat("\nDone.\n")
print(Sys.time())
cat("\n----------------------------------------------------------------\n")
