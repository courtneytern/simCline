#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=8:00:00
#SBATCH --partition=standard
#SBATCH --account=berglandlab
#SBATCH -o /scratch/cat7ep/slurmOut/deinterleave.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/deinterleave.%A_%a.err # Standard error
#SBATCH --array=591-773

## This script will reformat the interleaved file that was not split from fasterq-dump automatically
# The Signor files were not automatically split. the array will split them

module load gcc/9.2.0
module load bbmap

cd /scratch/cat7ep/fasterq

# get the appropriate accession number
sra=$( grep ^"${SLURM_ARRAY_TASK_ID}""," /scratch/cat7ep/simCline/biosampleresults/concatenated.csv | \
  awk -F"," '{
    split ($0,array,",")
    SRSnum= array[16]
    print SRSnum
  }' )

reformat.sh in=${sra}.fastq out1=${sra}_1.fq out2=${sra}_2.fq
