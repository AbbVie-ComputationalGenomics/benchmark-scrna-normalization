#!/usr/bin/env bash

#save the invoking command as a record
date >>cmd.log
echo $0 $* >>cmd.log

outpath=/mnt/nvme/users/wachtax/
filen=$1
normmet=$2

for f in `cat $filen`
do
        export dataset=${f}
        echo ${f}
        sbatch --mem=250G -N 1 -n 5 --job-name=${dataset}_lowrank --mail-type=ALL -o ${outpath}/lowrank_run/Output_files/${dataset}_lowrank.out -e ${outpath}/lowrank_run/Error_files/${dataset}_lowrank.err --wrap="Rscript --vanilla /mnt/nvme/users/wachtax/lowrank_run/ALRA_imp.R $f $normmet"

done



