#!/usr/bin/env bash

#SBATCH -n 1
##SBATCH --mail-type=END
##SBATCH --mail-user=hongkityeeclarice@wustl.edu

##Change number of arrays to number of files you need to run
#SBATCH --array=1-19%19

file=$( sed -n ${SLURM_ARRAY_TASK_ID}p R2_lookup.txt )
output=$( sed -n ${SLURM_ARRAY_TASK_ID}p basename_lookup.txt )

python3 /home/hongkityeeclarice/cohen_lab/genome_wide_LPs/emerald_expression/count_SB_emerald_barcodes.py ${file} CTCGAG ${output}
