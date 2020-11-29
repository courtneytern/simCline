#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH --partition=largemem
#SBATCH --account=berglandlab
#SBATCH --array=1-5

# This script will merge gVCFs into a unified database for genotype calling.
# Because this and following steps are so memory intensive, this will be done using
# a per chromosome approach

#Load Modules
module load gatk

#Name of pipeline
PIPELINE=GenomicsDBImport

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/yey2sn

# User defined inputs -- this represents the name of the samples
OW_sample_map=/scratch/yey2sn/Samples_to_haplotype.txt
#This file looks like this
#  sample1      sample1.vcf.gz
#  sample2      sample2.vcf.gz
#  sample3      sample3.vcf.gz

#Where the bam files are located
BAMS_FOLDER=/scratch/yey2sn/joint_bams

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
# Generate Folders and files
###########################################################################
###########################################################################
# this part of the script will check and generate, if necessary, all of the output folders used in the script

if [[ -d "TEMP_GenomicsDBImport" ]]
then
	echo "Working TEMP_MERGEVCF folder exist"
	echo "lets move on"
	date
else 
	echo "folder doesnt exist. lets fix that"
	mkdir $WORKING_FOLDER/TEMP_GenomicsDBImport_${i}
	date
fi

###########################################################################
###########################################################################
# Merge VCFs using GenomicsDBImport
###########################################################################
###########################################################################

  gatk --java-options "-Xmx${JAVAMEM} -Xms${JAVAMEM}" \
       GenomicsDBImport \
       --genomicsdb-workspace-path $WORKING_FOLDER/OVERWINTER_2018_2019_DBI_${i} \
       --batch-size 50 \
       --sample-name-map $OW_sample_map \
       --tmp-dir=$WORKING_FOLDER/TEMP_MERGEVCF_${i} \
       --reader-threads $CPU \
       -L ${i}

echo ${i} "done" $(date)
