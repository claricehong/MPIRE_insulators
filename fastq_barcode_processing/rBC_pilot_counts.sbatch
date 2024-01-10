#!/usr/bin/env bash

#SBATCH -n 1
#SBATCH --mail-type=END
#SBATCH --mail-user=hongkityeeclarice@wustl.edu
#SBATCH --mem=8196
#SBATCH --cpus-per-task=4

#SBATCH --array=1-4%4

R1=$( sed -n ${SLURM_ARRAY_TASK_ID}p R1_lookup.txt )
R2=$( sed -n ${SLURM_ARRAY_TASK_ID}p R2_lookup.txt )
basename=$( sed -n ${SLURM_ARRAY_TASK_ID}p basename_lookup.txt )

python3 /home/hongkityeeclarice/cohen_lab/genome_wide_LPs/recombined_expression/rBC_pilot_error_correction.py CTAGC TCGAG ${R1} ${R2} -o ${basename}
