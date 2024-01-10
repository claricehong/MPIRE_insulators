#!/usr/bin/env bash

#SBATCH --mail-type=END
#SBATCH --mail-user=hongkityeeclarice@htcf.wustl.edu
#SBATCH --mem=8196
#SBATCH --cpus-per-task=4

#SBATCH --array=1-12%4

eval `spack load --sh bwa`
eval `spack load --sh samtools`

read1=$( sed -n ${SLURM_ARRAY_TASK_ID}p basename_lookup.txt )

bwa mem /scratch/ref/genomes/Homo_sapiens/Ensembl/hg38_ensembl/hg38 ${read1}.fastq > ${read1}.sam

samtools view -bS ${read1}.sam > ${read1}.bam

samtools sort ${read1}.bam | samtools index - ${read1}.bam.bai

rm ${read1}.sam
