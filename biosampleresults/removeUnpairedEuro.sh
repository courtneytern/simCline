#!/bin/bash

#sort all fastq files, remove those without pairs (put in discard file?)
#then split by _1 and _2
fastqDir="/scratch/cat7ep/fastqEuro"

cd /scratch/cat7ep/simCline/biosampleresults
grep -v 'row' over5sim.csv | \
awk -F"," '{
  print $17 >> euroList.txt
}'

#sort fastq file
for line in euroList.txt; do
  cat ${fastqDir}/"$line"_1.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > ${fastqDir}/"$line"_1_sorted.fastq
  cat ${fastqDir}/"$line"_2.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > ${fastqDir}/"$line"_2_sorted.fastq
  awk 'NR % 4 == 0' "$line"_1_sorted.fastq> output
done
