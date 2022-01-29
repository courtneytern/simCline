#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=120G
#SBATCH --time=2:00:00
#SBATCH --partition=largemem
#SBATCH --account=berglandlab_standard
#SBATCH -o /scratch/cat7ep/slurmOut/combinePooledIndivid.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/combinePooledIndivid.%A_%a.err # Standard error

####### sbatch /scratch/cat7ep/simCline/biosampleresults/combinePooledIndivid.sh

module load gcc/7.1.0  openmpi/3.1.4
module load R

Rscript --vanilla /scratch/cat7ep/simCline/biosampleresults/combinePooledIndividSNPs.R