#!/bin/bash

# This script gets all of the fastq files for the individual samples.

cd /scratch/cat7ep/simCline/biosampleresults/

###make list of input bam files
#first take the pooled samples that don't include Jackson or Kang
sra=$( grep ",I," ./concatenated.csv | grep -vE "Jackson|Kang" | \
  awk -F"," '{
    split ($0,array,",")
    SRSnum= array[16]
    print SRSnum
  }' )

#move the individual files to a different folder
  for line in $sra; do
    echo /scratch/cat7ep/fasterq/$line* >> ./individFileNames.txt
  done

sed 's/ /\n/g' ./individFileNames.txt > ./individFileNames.txt
