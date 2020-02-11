#!/bin/sh

## set variable "file" to the parameter of findAdapter.sh
## file in format "/scratch/cat7ep/fastq/SRA###.fastq"
file=$1

## run fastq_detect.pl on given file
## grep the line with the version marked as the match
## split the line into just the name of the version
/scratch/cat7ep/simCline/biosampleresults/fastq_detect.pl ${file} > temp_out.txt

  grep 'x' temp_out.txt | \
  awk '{
    split($0,out," :")
    print out[1]
  }'
