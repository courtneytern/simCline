#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH --partition=largemem
#SBATCH --account=berglandlab
#SBATCH -o /scratch/cat7ep/slurmOut/genotype_genomics.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/genotype_genomics.%A_%a.err # Standard error
#SBATCH --array=1-5

# This script will conduct genotype calling on the GenomeDBI object

####### sbatch /scratch/cat7ep/simCline/biosampleresults/6.Genotype_GenomicsDB.sh

#Load Modules
module load gatk

#Name of pipeline
PIPELINE=GenotypeGVCFs

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/cat7ep/GenotypeGenomics

#Reference genome
REFERENCE=/project/berglandlab/courtney/simCline/refgenomes/simulans/dsim-mod.fasta

#Intervals to analyze
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
# Identify the Genome database to genotyoe
###########################################################################
###########################################################################

GenomeDB_path=`echo $WORKING_FOLDER/SIMCLINE_DBI_${i}`

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
