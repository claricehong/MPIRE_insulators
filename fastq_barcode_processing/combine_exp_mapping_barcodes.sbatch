#!/usr/bin/env bash

#SBATCH -n 1
#SBATCH --mail-type=END
#SBATCH --mail-user=hongkityeeclarice@wustl.edu
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1

#SBATCH --array=1-4%4

filepath='emerald_expression'

read exp map < <( sed -n ${SLURM_ARRAY_TASK_ID}p ${filepath}/correction_lookup.txt )

eval "$(conda shell.bash hook)"
conda activate trip_mapping

python3 /home/hongkityeeclarice/cohen_lab/genome_wide_LPs/correct_expression_barcodes/combine_exp_mapping_barcodes.py ${filepath}/${exp} ../${map} ${filepath}/${exp}
