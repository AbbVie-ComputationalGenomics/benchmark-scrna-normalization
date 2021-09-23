#!/usr/bin/env bash

#save invoking command as a record
date >>cmd.log
echo $0 $* >>cmd.log

outpath=/mnt/nvme/users/wachtax/
filen=$1

for f in `cat $filen`
do
   	export dataset=${f}
	echo ${f}
	sbatch --mem=250G -N 1 -n 5 --job-name=${dataset}_scImpute --mail-type=ALL -o ${outpath}/scImpute_run/Output_files/${dataset}_scImpute.out -e ${outpath}/scImpute_run/Error_files/${dataset}_scImpute.err --wrap="Rscript --vanilla /mnt/nvme/users/wachtax/scImpute_run/scImpute.R $f"
done

