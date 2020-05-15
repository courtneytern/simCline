#!/bin/sh
#
#SBATCH -J bam2fastq # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/cat7ep/slurmOut/bam2fastq.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/bam2fastq.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

##############################
##In rivanna, un-comment below
module load samtools/1.10

#sbatch /scratch/cat7ep/simCline/biosampleresults/bamtofastq.sh

#each line in destmapped.txt is a file name in /project/berglandlab/dest_mapped/
/scratch/cat7ep/simCline/destmapped.txt | \
while read -r filename; do
  samtools bam2fq /project/berglandlab/dest_mapped/"$filename"/sim.bam > /scratch/cat7ep/fastq/"$filename".fastq
done
