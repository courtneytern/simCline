#!/bin/sh

#groupget command pulls all files from one project number from ENA
alias enaGroupGet="/scratch/cat7ep/simCline/biosampleresults/enaBrowserTools/python3/enaGroupGet"

enaGroupGet -f fastq -d /scratch/cat7ep/bastidefastq PRJEB3292
