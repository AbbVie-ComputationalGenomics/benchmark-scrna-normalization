#!/usr/bin/env bash

all_data=(kolod darmanis pollen zeisel aml appps1 hca-bm hca-cb hipsc pbmc)

# mkdir scone_output_full
# mkdir scone_output_1k

export normdir="/home/pvijay/projects/internal/sc_methods_eval/normalization"

indir_full="/home/pvijay/projects/internal/sc_methods_eval/datasets/scone_input_sce"
indir_1k="/home/pvijay/projects/internal/sc_methods_eval/datasets/scone_input_sce_top1K_HVG"

# scaling="seurat,scnorm,linnorm,zinbwave,census,sum,tmm,uq,fq,deseq"
# scaling="seurat,scnorm,linnorm,tpm,census,sum,tmm,uq,fq,deseq"
scaling="sum,tmm,fq" ##test

for f in ${all_data[@]}
do
  echo $f
  export data=${f}
  
#   ##all genes
  #export outdir="scone_output_full/${data}"
#   export outdir="test_output_full/${data}"
#   mkdir ${outdir}
#   export infile=${indir_full}/${data}_sconeinput.RData
# 
# 	sbatch --mem=200G --job-name=${data}_full --wrap="./scone.R ${infile} ${data} ${normdir} --scaling ${scaling} --outdir ${outdir} >> ${outdir}/${data}_log.txt 2>&1"

  ##1k subset of genes
  # export outdir="scone_output_1k/${data}"
  export outdir="test_output_1k/${data}"
  mkdir ${outdir}
  export infile=${indir_1k}/${data}_1Kgenes_sconeinput.RData

	sbatch --mem=200G --job-name=${data}_1k --wrap="./scone.R ${infile} ${data} ${normdir} --subset --scaling ${scaling} --outdir ${outdir} >> ${outdir}/${data}_log.txt 2>&1"


done
