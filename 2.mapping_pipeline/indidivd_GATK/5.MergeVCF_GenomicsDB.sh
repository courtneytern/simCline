#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH --partition=largemem
#SBATCH --account=berglandlab
#SBATCH -o /scratch/cat7ep/slurmOut/mergeVCF_dbi.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/merrgeVCF_dbi.%A_%a.err # Standard error
#SBATCH --array=1-5

# This script will merge gVCFs into a unified database for genotype calling.
# Because this and following steps are so memory intensive, this will be done using
# a per chromosome approach

####### sbatch /scratch/cat7ep/simCline/biosampleresults/5.MergeVCF_GenomicsDB.sh


#Load Modules
module load gatk

#Name of pipeline
PIPELINE=GenomicsDBImport

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/cat7ep/individPipeline/MergeVCF

# User defined inputs -- this represents the name of the samples
sample_map=/scratch/cat7ep/individPipeline/Samples_to_haplotype.txt
#This file looks like this
#  sample1      sample1.vcf.gz
#  sample2      sample2.vcf.gz
#  sample3      sample3.vcf.gz
## second column is the path

#Intervals to analyze (chromosome arms)
intervals=/scratch/cat7ep/simCline/biosampleresults/intervals.txt

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

if [[ -d "TEMP_GenomicsDBImport"_$i ]]
then
	echo "Working TEMP_GenomicsDBImport_ folder exist"
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
       --genomicsdb-workspace-path $WORKING_FOLDER/SIMCLINE_DBI_${i} \
       --batch-size 50 \
       --sample-name-map $sample_map \
       --tmp-dir $WORKING_FOLDER/TEMP_GenomicsDBImport_${i} \
       --reader-threads $CPU \
       -L ${i}

echo ${i} "done" $(date)
