#!/usr/bin/env Rscript

cat("\n--------------------------LINNORM--------------------------------------\n\n")
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

cat("\nLoading required packages: Linnorm, SingleCellExperiment ... \n")
suppressPackageStartupMessages(require(Linnorm))
suppressPackageStartupMessages(library(SingleCellExperiment))

cat("\nReading input scone sce object ... \n")
print(infile)
load(infile)
print(dim(sce))

cat("\nRunning Linnorm ... \n")
r_linnorm = Linnorm(counts(sce))
r_linnorm = Linnorm.Norm(r_linnorm)

cat("\nWriting output to ",outfile," ... \n")
write.csv(r_linnorm, file=outfile)

cat("\nDone.\n")
print(Sys.time())
cat("\n----------------------------------------------------------------\n")
