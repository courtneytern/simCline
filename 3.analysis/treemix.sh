#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=berglandlab_standard
#SBATCH -o /scratch/cat7ep/slurmOut/treemix3.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/treemix3.%A_%a.err # Standard error

####### sbatch /scratch/cat7ep/simCline/3.analysis/treemix.sh


# setup and run treemix

# download tar.gz from Treemix source code
# cd into treemix folder

### load GSL package and boost to configure treemix
# module load gcc/9.2.0 intel/20.0 gsl/2.6
# module load intel/20.0  mvapich2/2.3.3 boost/1.68.0
# ./configure
# make
# export DESTDIR="$HOME/Software/LocalInstall" && make -j4 install

## reformat treemix input
# cd /scratch/cat7ep/simCline/data/
# cut -d" " -f3- ./TREEMIX_032522.txt | \
# awk -F" " '{
#   row=$0
#   gsub("NA,NA","0,0",row)
#   print row
# }' > ./TREEMIX_FINAL_032522.txt

module load gcc/9.2.0 intel/20.0 gsl/2.6
module load intel/20.0  mvapich2/2.3.3 boost/1.68.0

# Set paths
TREEMIX_PATH="$HOME/Software/LocalInstall/usr/local/bin/treemix"

# input must be gzipped
INPUT_FILE="/scratch/cat7ep/simCline/data/TREEMIX_FINAL_032522.txt.gz"
OUTGROUP="Madagascar.Joffreville.2002.NA"
# run for m=0,1,2,3
OUT_STEM="simcline_treeOut_032722"

cd /scratch/cat7ep/simCline/data/treemixOutputs

# m is number of migratiion events
# k is SNP grouping
$TREEMIX_PATH -i $INPUT_FILE -o $OUT_STEM -root $OUTGROUP -m 0 -k 1000

# make plot in treemixPlot.R
