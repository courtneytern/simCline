#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=72:00:00
#SBATCH --partition=standard
#SBATCH --account=berglandlab
#SBATCH -o /scratch/cat7ep/slurmOut/multi_qualimap.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/multi_qualimap.%A_%a.err # Standard error

## sbatch /scratch/cat7ep/simCline/biosampleresults/3.run_multi_qualimap.sh

JAVAMEM=28G
WORKING_DIRECTORY=/scratch/cat7ep/simCline/individPipeline

module load qualimap

qualimap  multi-bamqc  -outdir /scratch/cat7ep/simCline/individPipeline/multi_bamQC --java-mem-size=$JAVAMEM -d $WORKING_DIRECTORY/bam_qc_guide_file.txt
