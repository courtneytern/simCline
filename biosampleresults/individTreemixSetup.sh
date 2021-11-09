#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=48:00:00
#SBATCH --partition=standard
#SBATCH --account=berglandlab_standard
#SBATCH -o /scratch/cat7ep/slurmOut/individTreemixSetup.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/individTreemixSetup.%A_%a.err # Standard error

### sbatch /scratch/cat7ep/simCline/biosampleresults/individTreemixSetup.sh

module load gcc/7.1.0  openmpi/3.1.4
module load R

Rscript --vanilla /scratch/cat7ep/simCline/biosampleresults/individTreemixSetup.R
