#!/bin/bash

#This script filters the pooled sequencing data from the final bam files
#set directories
bamDir="/project/berglandlab/courtney/simCline/bamfinal"

cd /scratch/cat7ep/simCline/biosampleresults/

###make list of input bam files
#first take the pooled samples that don't include european
sra=$( grep ",P," ./concatenated.csv | grep -v "DrosEU" | \
  awk -F"," '{
    split ($0,array,",")
    SRSnum= array[16]
    print SRSnum
  }' )
#now get european
idEuro=$( grep ",P," ./concatenatedEuro.csv| \
  awk -F"," '{
    split ($0,array,",")
    id= array[16]
    print id
  }' )

for line in $sra; do
  echo $line >> /scratch/cat7ep/simCline/biosampleresults/inputPooledSamps.txt
  echo $bamDir/$line.finalmap.bam >> /scratch/cat7ep/simCline/biosampleresults/inputPooledBam.txt
done
for line in $idEuro; do
  echo $line >> /scratch/cat7ep/simCline/biosampleresults/inputPooledSamps.txt
  echo $bamDir/$line.finalmap.bam >> /scratch/cat7ep/simCline/biosampleresults/inputPooledBam.txt
done

#go to makeVCF.sh
