#!/bin/bash

#This script filters the pooled sequencing data from the final bam files
#set directories
bamDir="/project/berglandlab/courtney/simCline/bamfinal"

cd /scratch/cat7ep/simCline

###make list of input bam files
#first take the pooled samples that don't include european
sra=$( grep ",P," ./metadata/concatenated.csv | grep -v "DrosEU" | \
  awk -F"," '{
    split ($0,array,",")
    SRSnum= array[16]
    print SRSnum
  }' )
#now get european
idEuro=$( grep ",P," ./metadata/concatenatedEuro.csv | \
  awk -F"," '{
    split ($0,array,",")
    id= array[16]
    print id
  }' )

for line in $sra; do
  echo $bamDir/$line.finalmap.bam >> ./2.mapping_pipeline/inputPooledBam.txt
done
for line in $idEuro; do
  echo $bamDir/$line.finalmap.bam >> ./2.mapping_pipeline/inputPooledBam.txt
done

##now create the sample name list
samps=$( grep ",P," ./metadata/concatenated.csv | \
  awk -F"," '{
    split ($0,array,",")
    name= array[17]
    print name
  }' )

for line in $samps; do
  echo $line >> ./2.mapping_pipeline/inputPooledSamps.txt
done

#go to mapping script 
