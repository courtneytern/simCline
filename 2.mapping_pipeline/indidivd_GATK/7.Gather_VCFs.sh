#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --partition=largemem
#SBATCH --account=berglandlab
#SBATCH -o /scratch/cat7ep/slurmOut/gatherVCFs.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/gatherVCFs.%A_%a.err # Standard error

# This script is a pipeline which gather VCFs from all chromosomes.
### Using genotyped.raw.vcf.gz because there's no recalibration step for simulans

####### sbatch /scratch/cat7ep/simCline/biosampleresults/8.Gather_VCFs.sh

#Load Modules
module load picard
module load tabix

#Name of pipeline
PIPELINE=Simcline_final_2021

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/cat7ep/individPipeline/MergeVCF

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
  I=$WORKING_FOLDER/Dsim_Scf_2L.genotyped.raw.vcf.gz \
  I=$WORKING_FOLDER/Dsim_Scf_2R.genotyped.raw.vcf.gz \
  I=$WORKING_FOLDER/Dsim_Scf_3L.genotyped.raw.vcf.gz \
  I=$WORKING_FOLDER/Dsim_Scf_3R.genotyped.raw.vcf.gz \
  I=$WORKING_FOLDER/Dsim_Scf_X.genotyped.raw.vcf.gz \
  O=$WORKING_FOLDER/$PIPELINE.vcf

#bgzip and tabix
	bgzip $WORKING_FOLDER/$PIPELINE.vcf
	tabix $WORKING_FOLDER/$PIPELINE.vcf.gz


echo "done" $(date)
