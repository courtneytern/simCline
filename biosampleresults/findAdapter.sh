## sort through just the fourth line of each read (fourth line is quality line)
awk '{if(NR%4==0){print $0}}' ERR1980726.fastq > /scratch/cat7ep/simCline/biosampleresults/ERR1980726out.txt

## Write an array of the ASCII scale of the quality markers
## Search through each row of quality, get the range
## match the range to certain read (with confidence interval?)
