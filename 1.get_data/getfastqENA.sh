#!/bin/sh
#
#SBATCH -J ENAfastq # A single job name for the array
#SBATCH --ntasks-per-node=20 # twenty cores
#SBATCH -N 1
#SBATCH -t 48:00:00 ### 24 hr
#SBATCH --mem 0
#SBATCH -o /scratch/cat7ep/slurmOut/ENAfastq2.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/ENAfastq2.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load gcc/9.2.0
module load mvapich2/2.3.3
module load openmpi/3.1.6
module load intel/20.0
module load intelmpi/20.0
module load python/3.7.7

#groupget command pulls all files from one project number from ENA
alias enaGroupGet="/scratch/cat7ep/simCline/software/enaBrowserTools/python3/enaGroupGet"

enaGroupGet -f fastq -d /scratch/cat7ep/fastqJackson PRJEB7673
