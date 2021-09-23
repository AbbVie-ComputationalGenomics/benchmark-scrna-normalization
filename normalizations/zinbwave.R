#!/usr/bin/env Rscript

cat("\n--------------------------ZINBWAVE--------------------------------------\n\n")
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

cat("\nLoading required packages: zinbwave, SingleCellExperiment, Matrix\n")
suppressPackageStartupMessages(require(zinbwave))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(Matrix))

cat("\nReading input scone sce object ... \n")
load(infile)
print(dim(sce))

cat("\nFiltering cells with only 0 counts ... \n")
sums = colSums(counts(sce))
sce = sce[,which(sums>0)]
print(dim(sce))

cat("\nRunning zinbwave ... \n")
counts(sce) = as.matrix(counts(sce))
z = zinbwave(sce)
print(z)
mat = assay(z, "normalizedValues")
print(dim(mat))
print(mat[1:3,1:3])

cat("\nWriting output to ",outfile," ... \n")
write.csv(mat, file=outfile)

cat("\nDone.\n")
print(Sys.time())
cat("\n----------------------------------------------------------------\n")
