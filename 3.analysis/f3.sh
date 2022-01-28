#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=72:00:00
#SBATCH --partition=standard
#SBATCH --account=berglandlab_standard
#SBATCH -o /scratch/cat7ep/slurmOut/f3FINAL.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/f3FINAL.%A_%a.err # Standard error
#SBATCH --array=5-12,21,22,24-32

####### sbatch /scratch/cat7ep/simCline/biosampleresults/f3.sh
## This will take each of the smaller treemix files and run f3 stats on them
## The outputs from this file will need to be parsed through to get the A;B,C trees of interest

# setup and run treemix

# download tar.gz from Treemix source code
# cd into treemix folder

### load GSL package and boost to configure treemix
# module load gcc/9.2.0 intel/20.0 gsl/2.6
# module load intel/20.0  mvapich2/2.3.3 boost/1.68.0
# ./configure
# make
# export DESTDIR="$HOME/Software/LocalInstall" && make -j4 install

module load gcc/9.2.0 intel/20.0 gsl/2.6
module load intel/20.0  mvapich2/2.3.3 boost/1.68.0

# Set paths
THREEPOP_PATH="$HOME/Software/LocalInstall/usr/local/bin/threepop"

INPUT_FOLDER=/scratch/cat7ep/simCline/biosampleresults/treemixInputs
INPUT_FILE=${INPUT_FOLDER}/treemixInput_${SLURM_ARRAY_TASK_ID}.txt.gz

OUTPUT_FOLDER=/scratch/cat7ep/simCline/biosampleresults/f3_outputs
OUTPUT_FILE=${OUTPUT_FOLDER}/f3_output_${SLURM_ARRAY_TASK_ID}.txt

cd /scratch/cat7ep/simCline/biosampleresults

# run f3 statistics and save into out path
$THREEPOP_PATH -i $INPUT_FILE -k 500 > $OUTPUT_FILE
