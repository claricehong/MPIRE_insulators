#!/usr/bin/env bash

#SBATCH -n 1
#SBATCH --mail-type=END
#SBATCH --mail-user=hongkityeeclarice@wustl.edu
#SBATCH --mem=8196
#SBATCH --cpus-per-task=4

#SBATCH --array=1-24%8

R1=$( sed -n ${SLURM_ARRAY_TASK_ID}p R1_lookup.txt )
R2=$( sed -n ${SLURM_ARRAY_TASK_ID}p R2_lookup.txt )
basename=$( sed -n ${SLURM_ARRAY_TASK_ID}p basename_lookup.txt )

python3 /home/hongkityeeclarice/cohen_lab/genome_wide_LPs/recombined_expression/bc3_counts.py GCTAGC CTCGAG ${R1} ${R2} insulator_bcs.txt -i -o grouped_counts/${basename}