#!/usr/bin/env bash

#SBATCH -n 1
##SBATCH --mail-type=END
#sleep 1 # pause for 1 second between job submissions

##Change number of arrays to number of files you need to run
#SBATCH --array=1-12%4

R1=$( sed -n ${SLURM_ARRAY_TASK_ID}p R1_lookup.txt )
R2=$( sed -n ${SLURM_ARRAY_TASK_ID}p R2_lookup.txt )
basename=$( sed -n ${SLURM_ARRAY_TASK_ID}p basename_lookup.txt )


python3 /home/hongkityeeclarice/cohen_lab/genome_wide_LPs/mapping/iPCR_mapping/process_fastq.py CAACTG CTCGAG ${R1} ${R2} -o ${basename}
