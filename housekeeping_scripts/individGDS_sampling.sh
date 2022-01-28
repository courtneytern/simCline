#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=72:00:00
#SBATCH --partition=standard
#SBATCH --account=berglandlab
#SBATCH -o /scratch/cat7ep/slurmOut/individGDS.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/individGDS.%A_%a.err # Standard error

####### sbatch /scratch/cat7ep/simCline/biosampleresults/individGDS_sampling.sh

module load gcc/7.1.0  openmpi/3.1.4  intel/18.0  intelmpi/18.0
module load R
R

setwd("~/Downloads/GitHub/simCline/biosampleresults/")
