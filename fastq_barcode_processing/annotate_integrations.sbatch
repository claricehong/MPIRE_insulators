#!/usr/bin/env bash

#SBATCH --mail-type=END
#SBATCH --mail-user=hongkityeeclarice@htcf.wustl.edu
#SBATCH --mem=8196
#SBATCH --cpus-per-task=4

#SBATCH --array=1-12%4

eval `spack load --sh bedtools2`

eval "$(conda shell.bash hook)"
conda activate trip_mapping

read1=$( sed -n ${SLURM_ARRAY_TASK_ID}p basename_lookup.txt )

python3 /home/hongkityeeclarice/cohen_lab/genome_wide_LPs/mapping/iPCR_mapping/annotate_integrations.py ${read1}.bam ${read1}

conda deactivate
