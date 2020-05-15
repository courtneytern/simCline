#!/bin/sh
#
#SBATCH -J bam2fastq # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/cat7ep/slurmOut/bam2fastq.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/bam2fastq.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

##############################
##In rivanna, un-comment below
module load samtools/1.10

##SLURM_ARRAY_TASK_ID=4
#sbatch /scratch/cat7ep/simCline/biosampleresults/bamtofastq.sh

#each line in destmapped.txt is a file name in /project/berglandlab/dest_mapped/
/scratch/cat7ep/simCline/destmapped.txt | \
while read filename; do
  samtools bam2fq  ./*/sim.bam > /scratch/cat7ep/fastq/*.fastq
done
