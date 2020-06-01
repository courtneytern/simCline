#!/bin/sh
#
#SBATCH -J simCline # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/cat7ep/slurmOut/simCline.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/simCline.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

##############################
##In rivanna, un-comment below
module load sratoolkit/2.9.1

##SLURM_ARRAY_TASK_ID=4
#sbatch --array=1-811 /scratch/cat7ep/simCline/biosampleresults/getfastq.sh
###array=1-811 for total concatenated.csv
echo ${SLURM_ARRAY_TASK_ID}

##takes in parameter of the row number to access
SRS=$( grep ^${SLURM_ARRAY_TASK_ID}"," /scratch/cat7ep/simCline/biosampleresults/concatenated.csv | \
  awk -F"," '{
    split ($0,array,",")
    SRSnum= array[16]
    print SRSnum
  }' )

fasterq-dump ${SRS} -t /scratch/cat7ep/interDir --outdir /scratch/cat7ep/fasterq/
