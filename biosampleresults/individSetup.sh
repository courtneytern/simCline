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

  for line in /scratch/cat7ep/simCline/biosampleresults/individFileNames.txt; do
    mv $line* /scratch/cat7ep/individFastq
  done

sed 's/ /\n/g' ./individFileNames.txt > ./individFileNames.txt
