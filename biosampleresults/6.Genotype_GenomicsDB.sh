#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH --partition=largemem
#SBATCH --account=berglandlab
#SBATCH --array=1-5

# This script will conduct genotype calling on the GenomeDBI object

#Load Modules
module load gatk

#Name of pipeline
PIPELINE=GenotypeGVCFs

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/yey2sn

#Reference genome
REFERENCE=/project/berglandlab/Dmel_fasta_refences/holo_dmel_6.12.fa

#Intervals to analyze
intervals=/scratch/yey2sn/Intervals_Dmel.txt

#Parameters

#Java
JAVAMEM=190G
CPU=4

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
# Generate Folders and files
###########################################################################
###########################################################################

#Interval to analyze
i=`sed -n ${SLURM_ARRAY_TASK_ID}p $intervals`

echo ${i} "is being processed" $(date)

###########################################################################
###########################################################################
# Identify the Genome database to genotyoe
###########################################################################
###########################################################################

GenomeDB_path=`echo $WORKING_FOLDER/OVERWINTER_2018_2019_DBI_${i}`

###########################################################################
###########################################################################
# Genotype call the samples in the DBI merged set
###########################################################################
###########################################################################

 gatk --java-options "-Xmx${JAVAMEM}" GenotypeGVCFs \
   -R $REFERENCE \
   -V gendb://$GenomeDB_path \
   -O $WORKING_FOLDER/${i}.genotyped.raw.vcf.gz

echo ${i} "done" $(date)
