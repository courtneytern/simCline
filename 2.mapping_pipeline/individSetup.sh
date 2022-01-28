#!/bin/bash

# This script makes all of the files required for the individual-samples GATK pipeline to be run smoothly.
# This includes reference txt files, sample lists, and other setup files

cd /scratch/cat7ep/simCline/biosampleresults/

##########################
### PREP #################
##########################

#first take the individual samples
sra=$( grep ",I," ./concatenated.csv | \
  awk -F"," '{
    split ($0,array,",")
    SRSnum= array[16]
    print SRSnum
  }' )

#make file with list of sample names
  for line in $sra; do
    echo $line >> ./individFileNames.txt
  done

# file with the file names for all individual samples
sed 's/ /\n/g' ./individFileNames.txt > ./individFileNames.txt

#########################
### Step 3 ##############
#########################

# file with bam qc file paths for step 3 run_multi_qualimap
WORKING_FOLDER=/scratch/cat7ep/individPipeline/TrimMap
WORKING_DIRECTORY=/scratch/cat7ep/individPipeline

  for line in $sra; do
    echo -e $line'\t'$WORKING_FOLDER/joint_bams_qualimap/Qualimap_JointBam_${line} >> $WORKING_DIRECTORY/bam_qc_guide_file.txt
  done

#########################
### Step 5 ##############
#########################

# file with vcf sample paths for step 5 MergeVCF_genomics
  WORKING_FOLDER=/scratch/cat7ep/individPipeline/HapCaller

    for line in $sra; do
      echo -e $line'\t'$WORKING_FOLDER/haplotype_calling/${line}.raw.g.vcf.gz >> $WORKING_DIRECTORY/Samples_to_haplotype.txt
    done

# # Generate interval file
# module load gatk
#
# REF_GENOME_PATH=/project/berglandlab/courtney/simCline/refgenomes/simulans
# gatk PreprocessIntervals \
#       -R $REF_GENOME_PATH/dsim-mod.fasta \
#       --bin-length 0 \
#       --padding 0 \
#       -O /scratch/cat7ep/individPipeline/preprocessed_intervals.interval_list
