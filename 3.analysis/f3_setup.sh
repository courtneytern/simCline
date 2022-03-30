#! /bin/bash

## Separate combined treemix input into separate files to be run through f3 more efficiently
## Will make 19 files; each file has one AMR pop + all eur and afr pops

WORKING_FOLDER=/scratch/cat7ep/simCline/data/treemixInputs
cd $WORKING_FOLDER

combinedInput=/scratch/cat7ep/simCline/data/TREEMIX_FINAL_032522.txt

#cols 1,13,14,17-32 are N.Am
# all files will keep 2-12,15,16
colArray=(1 13 14 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32)
for c in ${colArray[@]}; do
  awk -F" " -v col=$c '{ print $col,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$15,$16 }'  < $combinedInput > $WORKING_FOLDER/treemixInput_${c}.txt
done
# replace col and the first item in the print list with cols 5-12, 21, 22, 24-32 one at a time
col=32
awk -F" " '{ print $32,$1,$2,$3,$4,$13,$14,$15,$16,$17,$18,$19,$20,$23 }' < $combinedInput > $WORKING_FOLDER/treemixInput_${col}.txt

## after you've generated all files:
gzip ./treemixInput_*

####
## RUN f3.sh FIRST, THEN RUN LINES BELOW


######################
## Post running f3.sh:
## Grep just the A;B,C trees from all the files and cat them together
## the output will be run through f3_parse.R

cd /scratch/cat7ep/simCline/data/f3_outputs
cat ./f3_output_* | grep ';' > ./FINAL_f3_output.txt
