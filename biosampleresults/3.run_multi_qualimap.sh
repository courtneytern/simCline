#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=72:00:00
#SBATCH --partition=standard
#SBATCH --account=berglandlab

JAVAMEM=28G
WORKING_DIRECTORY=/scratch/cat7ep/simCline/biosampleresults

module load qualimap

qualimap  multi-bamqc  -outdir /scratch/cat7ep/multi_bamQC --java-mem-size=$JAVAMEM -d $WORKING_DIRECTORY/bam_qc_guide_file.txt
