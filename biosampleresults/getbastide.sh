#!/bin/sh
module load gcc/9.2.0
module load mvapich2/2.3.3
module load gcc/9.2.0
module load openmpi/3.1.4
module load python/3.7.7

#groupget command pulls all files from one project number from ENA
alias enaGroupGet="/scratch/cat7ep/simCline/biosampleresults/enaBrowserTools/python3/enaGroupGet"

enaGroupGet -f fastq -d /scratch/cat7ep/bastidefastq PRJEB3292
