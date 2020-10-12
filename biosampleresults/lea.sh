#!/bin/sh
#
#SBATCH -J lea # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 30G
#SBATCH -o /scratch/cat7ep/slurmOut/lea.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/lea.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#load dependencies for R
module load gcc/7.1.0
module load openmpi/3.1.4
module load gcc/8.3.0
module load cuda/10.2.89
module load intel/18.0
module load intelmpi/18.0
module load R/4.0.0

#run lea.R
#right now it's just creating the lfmm file
Rscript lea.R
