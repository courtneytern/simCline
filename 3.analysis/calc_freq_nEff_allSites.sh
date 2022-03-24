#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --partition=largemem
#SBATCH --account=berglandlab_standard
#SBATCH -o /scratch/cat7ep/slurmOut/individNeff.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/individNeff.%A_%a.err # Standard error

####### sbatch /scratch/cat7ep/simCline/3.analysis/calc_freq_nEff_allSites.sh

module load gcc/7.1.0  openmpi/3.1.4
module load R

date
echo "Starting"

Rscript --vanilla /scratch/cat7ep/simCline/3.analysis/calc_freq_nEff_allSites.R

echo "done"
date
