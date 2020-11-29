#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=72:00:00
#SBATCH --partition=standard
#SBATCH --account=berglandlab

JAVAMEM=28G
WORKING_DIRECTORY=/scratch/yey2sn

module load qualimap

qualimap  multi-bamqc  -outdir /scratch/yey2sn/multi_bamQC --java-mem-size=$JAVAMEM -d $WORKING_DIRECTORY/bam_qc_guide_file.txt