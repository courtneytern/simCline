#!/bin/bash

# This script gets all of the fastq files for the individual samples.

cd /scratch/cat7ep/simCline/biosampleresults/

###make list of input bam files
# #first take the pooled samples that don't include Jackson or Kang
# sra=$( grep ",I," ./concatenated.csv | grep -vE "Jackson|Kang" | \
#   awk -F"," '{
#     split ($0,array,",")
#     SRSnum= array[16]
#     print SRSnum
#   }' )

#first take the individual samples
sra=$( grep ",I," ./concatenated.csv | \
  awk -F"," '{
    split ($0,array,",")
    SRSnum= array[16]
    print SRSnum
  }' )

#move the individual files to a different folder
  for line in $sra; do
    echo $line >> ./individFileNames.txt
  done

# file with the file names for all individual samples
sed 's/ /\n/g' ./individFileNames.txt > ./individFileNames.txt

# file with bam qc file paths for step 3 run_multi_qualimap
WORKING_FOLDER=/scratch/cat7ep/individPipeline/TrimMap
WORKING_DIRECTORY=/scratch/cat7ep/individPipeline

  for line in $sra; do
    echo -e $line'\t'$WORKING_FOLDER/joint_bams_qualimap/Qualimap_JointBam_${line} >> $WORKING_DIRECTORY/bam_qc_guide_file.txt
  done
