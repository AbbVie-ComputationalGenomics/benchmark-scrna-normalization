#!/usr/bin/env Rscript

cat("\n--------------------------CENSUS--------------------------------------\n\n")
print(Sys.time())
cat("\n")

# required args: input, output, gene lengths file
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=2){
  stop("\nRequired inputs: [input TPM matrix] [name of output file]")
}
infile = args[1] ## should be TPM matrix
outfile = args[2] 

suppressPackageStartupMessages(require(monocle))

cat("\nReading input ... \n")
print(infile)
tpm = read.csv(infile, row.names=1)
print(dim(tpm))

cat("\nCreating monocle object ... \n")
sampleann = data.frame(row.names = colnames(tpm), Cell = colnames(tpm))
geneinfo = data.frame(row.names = rownames(tpm), gene_short_name = rownames(tpm))
pd = new("AnnotatedDataFrame", data = sampleann)
fd = new("AnnotatedDataFrame", data = geneinfo)
cds = newCellDataSet(as.matrix(tpm),
                     phenoData = pd,
                     featureData = fd,
                     lowerDetectionLimit = 0,
                     expressionFamily = tobit(Lower = 0.1))
print(cds)

cat("\nPerforming census normalization ... \n")
rpc_matrix = relative2abs(cds, method = "num_genes")
print(dim(rpc_matrix))
print(rpc_matrix[1:3,1:3])

cat("\nWriting output to ",outfile,"\n")
write.csv(rpc_matrix, file=outfile)

cat("\nDone.\n")
print(Sys.time())
cat("\n----------------------------------------------------------------\n")