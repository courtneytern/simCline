#!/bin/sh

## set variable "file" to the parameter of findAdapter.sh
## file in format "SRA###.fastq"
file=$1

## run fastq_detect.pl on given file
############## Uncomment below for testing on mac
# ~/Downloads/GitHub/simCline/biosampleresults/fastq_detect.pl ~/Downloads/GitHub/simCline/biosampleresults/${file} > temp_out.txt
############## Uncomment below for running on rivanna
/scratch/cat7ep/simCline/biosampleresults/fastq_detect.pl /scratch/cat7ep/fastq/${file} > temp_out.txt

## grep the line with the version marked as the match
## split the line into just the name of the version and print
## store version in var
version={  grep ' x ' temp_out.txt | \
  awk '{
    split($0,out," :")
    print out[1]
  }' }

## translate
if ${version}= "Sanger" | "Illumina 1.8+"{echo 33}
else {echo 64} ## is +64; change to +33
