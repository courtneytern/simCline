#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --partition=largemem
#SBATCH --account=berglandlab

# This script is a pipeline which gather VCFs from all chromosomes.

#Load Modules
module load picard
module load tabix

#Name of pipeline
PIPELINE=OW_final_2018_2019

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/yey2sn

#Parameters
#Java
JAVAMEM=110G
CPU=2

###########################################################################
###########################################################################
# Begin Pipeline
###########################################################################
###########################################################################
#This part of the pipeline will generate log files to record warnings and completion status

# Move to working directory
cd $WORKING_FOLDER

###########################################################################
###########################################################################
# Gather VCFs to make a final VCF
###########################################################################
###########################################################################

# order 2L, 2R, 3L, 3R, 4, X, Y

java -Xmx$JAVAMEM \
 -jar $PICARD GatherVcfs \
  I=$WORKING_FOLDER/2L.recalibratedSNP.vcf.gz \
  I=$WORKING_FOLDER/2R.recalibratedSNP.vcf.gz \
  I=$WORKING_FOLDER/3L.recalibratedSNP.vcf.gz \
  I=$WORKING_FOLDER/3R.recalibratedSNP.vcf.gz \
  I=$WORKING_FOLDER/X.recalibratedSNP.vcf.gz \
  O=$WORKING_FOLDER/$PIPELINE.vcf
	
#bgzip and tabix
	bgzip $WORKING_FOLDER/$PIPELINE.vcf
	tabix $WORKING_FOLDER/$PIPELINE.vcf.gz
	

echo "done" $(date)