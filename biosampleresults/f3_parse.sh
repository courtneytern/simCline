#! /bin/bash

## Parse through f3 output

f3_output="/scratch/cat7ep/simCline/biosampleresults/delete_later/f3_output.txt"

grep ';' $f3_output | \
awk -F" " '{
  
}'
