#!/usr/bin/env Rscript

cat("\n--------------------------TPM--------------------------------------\n\n")
print(Sys.time())
cat("\n")

# required args: input, output, gene lengths file
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=3){
  stop("\nRequired inputs: [input sce .RData object] [name of output file] [gene lengths csv]")
}
infile = args[1] ## should be sce object
outfile = args[2] 

cat("\nLoading required packages: SingleCellExperiment\n")
suppressPackageStartupMessages(library(SingleCellExperiment))

cat("\nLoading input sce ... \n")
print(infile)
load(infile)
print(sce)

cat("\nReading gene lengths ... \n")
genelengths = read.csv(args[3],row.names=1)
print(dim(genelengths))
print(head(genelengths))

cat("\nCalculating TPM ... \n")

getTPM = function(counts, genelengths){
  ##helper function for CENSUS_FN
  ##input 1 is a matrix of raw counts where
  ##rows are genes with gene names specified as row names and
  ##samples are column with sample names as the header
  ##input 2 has row names as gene names and a column of gene lengths (bp)
  
  ##1. divide read counts by length of gene
  geneswdata = rownames(counts)[which(rownames(counts) %in% rownames(genelengths))]
  lengthfilt = data.frame(row.names = geneswdata,
                          length = genelengths[geneswdata,1])
  lengthkb = lengthfilt/1000
  RPK = counts[geneswdata,]/lengthkb$length
  
  ##2. calculate per million scaling factor
  pm = colSums(RPK)/1000000
  
  ##3. divide by scaling factor
  TPM = t(t(RPK)/pm)
  return(TPM)
}
tpm = getTPM(counts(sce), genelengths = genelengths)
print(dim(tpm))
print(tpm[1:3,1:3])

cat("\nWriting output to ",outfile,"\n")
write.csv(tpm, file=outfile)

cat("\nDone.\n")
print(Sys.time())
cat("\n----------------------------------------------------------------\n")