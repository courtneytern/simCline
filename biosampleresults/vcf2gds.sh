#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --partition=largemem
#SBATCH --account=berglandlab
#SBATCH -o /scratch/cat7ep/slurmOut/vcf2gds100621.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/vcf2gds100621.%A_%a.err # Standard error

module load gcc/7.1.0  openmpi/3.1.4  intel/18.0  intelmpi/18.0
module load R

/scratch/cat7ep/simCline/biosampleresults/vcf2gds.R
