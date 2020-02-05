#!/bin/sh
#
#SBATCH -J simCline # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/cat7ep/slurmOut/simCline.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/simCline.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#####srun ./getfastq.sh $SLURM_ARRAY_TASK_ID

##############################
##In rivanna, un-comment below
module load sratoolkit/2.9.1

##SLURM_ARRAY_TASK_ID=4
#sbatch --array=1-19 /scratch/cat7ep/simCline/biosampleresults/getfastq.sh
###array=1-1206 for total concatenated, 1-19 for bastide only
echo ${SLURM_ARRAY_TASK_ID}

##takes in parameter of the row number to access
SRS=$( grep ^${SLURM_ARRAY_TASK_ID}"," ~/simCline/biosampleresults/bastidecat.csv | \
  awk -F"," '{
    split ($0,array,",")
    SRSnum= array[15]
      #name= array[2]""array[6]""array[7]""array[11]
      #print name
    print SRSnum
  }' )

  ##only need fastq-dump command by itself in rivanna
  ##~/Downloads/sratoolkit.2.9.6-1-mac64/bin/fastq-dump -X 5 -Z

fastq-dump --split-files ${SRS} > /scratch/cat7ep/bastidefastq/${SRS}.fastq
