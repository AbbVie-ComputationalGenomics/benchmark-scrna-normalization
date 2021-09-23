#!/usr/bin/env Rscript

##########################################################################################
### scTransform
### Author: Sunantha Sethuraman
### Date: August 2019
### Usage: Rscript sct.R <input file> <output path>
### Description: Runs scTransform on single cell data. Input should be a SCE object file.
##########################################################################################

##########################################################################################
### Load libraries
##########################################################################################

suppressMessages(library(sctransform))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(data.table))

##########################################################################################
### Parse arguments
##########################################################################################

args = commandArgs(trailingOnly=TRUE)

scfile <- args[1]
outdir <- args[2]
fn <- tools::file_path_sans_ext(basename(scfile))

##########################################################################################
### Set up options and log
##########################################################################################

date <- gsub("-", "", as.character(Sys.Date()))
options(future.globals.maxSize= 891289600)
con <- file(paste0(outdir, "/", fn, "_", date, "_sct.log"))

sink(con, append=TRUE)
sink(con, append=TRUE, type="message")

##########################################################################################
### scTransform
##########################################################################################

print(fn)
print(Sys.time())

## Load and retrieve counts
load(scfile)
sccounts <- tryCatch( {counts(sce)}, error = function(e) {return(counts(sce_filt))})
roword <- row.names(sccounts) # save row order

## Run scTransform and convert back to counts
if(class(sccounts) == "data.frame") {sccounts <- as.matrix(sccounts)}
sccounts_vst <- sctransform::vst(sccounts, return_gene_attr = F)
norm <- sccounts_vst[[1]]
normcounts <- 10^norm # convert to normalized counts

## Add back filtered genes with 0 to keep consistent dimensions
genes_dropped <- rownames(sccounts)[! rownames(sccounts) %in% rownames(normcounts)]
if(length(genes_dropped) != 0) {
    zerodf <- data.frame(matrix(0, nrow = length(genes_dropped), ncol = ncol(normcounts)))
    rownames(zerodf) <- genes_dropped
    colnames(zerodf) <- colnames(normcounts)
    normcounts <- rbind(normcounts, zerodf)
  } 

## Write to file
normcounts <- normcounts[roword,] # return matrix in the original row order
normcounts <- data.frame(normcounts)
fwrite(normcounts, paste0(outdir, "/", fn, "_sct.csv"), row.names = T)

print(Sys.time())

##########################################################################################
### Restore output to console
##########################################################################################

sink() 
sink(type="message")
