#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=72:00:00
#SBATCH --partition=standard
#SBATCH --account=berglandlab_standard
#SBATCH -o /scratch/cat7ep/slurmOut/f3.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/f3.%A_%a.err # Standard error

####### sbatch /scratch/cat7ep/simCline/biosampleresults/f3.sh


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

INPUT_FILE="/scratch/cat7ep/simCline/biosampleresults/treemixCOMBINEDinput.txt.gz"
OUTPUT_FILE="/scratch/cat7ep/simCline/biosampleresults/f3_output_72H.txt"

cd /scratch/cat7ep/simCline/biosampleresults

# run f3 statistics and save into out path
$THREEPOP_PATH -i $INPUT_FILE -k 500 > $OUTPUT_FILE
