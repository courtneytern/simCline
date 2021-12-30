#! /bin/bash

## Separate combined treemix input into separate files to be run through f3 more efficiently
## Will make 19 files; each file has one AMR pop + all eur and afr pops

WORKING_FOLDER=/scratch/cat7ep/simCline/biosampleresults/treemixInputs
cd $WORKING_FOLDER

combinedInput=$WORKING_FOLDER/treemixCOMBINEDinput.txt
#f3_output="/scratch/cat7ep/simCline/biosampleresults/delete_later/f3_output.txt"

#cols 5-12 , 21, 22, 24-32 are N.Am
# all files will keep 1-4, 13-20, 23
col=32
awk -F" " '{ print $32,$1,$2,$3,$4,$13,$14,$15,$16,$17,$18,$19,$20,$23 }' < $combinedInput > $WORKING_FOLDER/treemixInput_${col}.txt
