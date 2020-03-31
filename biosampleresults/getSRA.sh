#!/bin/sh
#
#SBATCH -J simCline # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1
#SBATCH -t 30:00 ### 30 min
#SBATCH --mem 1G
#SBATCH -o /scratch/cat7ep/slurmOut/simCline.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/simCline.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


##SLURM_ARRAY_TASK_ID=4
#sbatch --array=1-826 /scratch/cat7ep/simCline/biosampleresults/getSRA.sh
echo ${SLURM_ARRAY_TASK_ID}

##takes in parameter of the row number to access
SRA=$( grep ^${SLURM_ARRAY_TASK_ID}"," /scratch/cat7ep/simCline/biosampleresults/concatenated.csv | \
  awk -F"," '{
    split ($0,array,",")
    SRSnum= array[15]
    print SRSnum
  }' )

cp /scratch/cat7ep/fastq/${SRA}*.fastq /project/berglandlab/courtney/simCline/fastq
