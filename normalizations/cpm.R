#!/usr/bin/env Rscript

cat("\n--------------------------CPM--------------------------------------\n\n")
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

cat("\nLoading required packages: SingleCellExperiment ... \n")
suppressPackageStartupMessages(library(SingleCellExperiment))

cat("\nReading input scone sce object ... \n")
print(infile)
load(infile)
print(dim(sce))

cat("\nCalculating CPM ... \n")
eo = 1000000*counts(sce)/colSums(counts(sce))
print(dim(eo))
print(eo[1:3,1:3])

cat("\nWriting output to ",outfile," ... \n")
write.csv(eo, file=outfile)

cat("\nDone.\n")
print(Sys.time())
cat("\n----------------------------------------------------------------\n")
