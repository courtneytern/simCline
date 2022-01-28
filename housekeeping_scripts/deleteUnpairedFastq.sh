## delete unpaired Machado fastq files
cd /scratch/cat7ep/simCline/biosampleresults
fastqDir=/scratch/cat7ep/fasterq

#get the SRS numbers from the Machado files only
machado=$( grep "Machado" ./concatenated.csv | \
  awk -F"," '{
    split ($0,array,",")
    SRSnum= array[16]
    print SRSnum
  }' )

for line in $machado; do
  if [[ -e "${fastqDir}/${line}.fastq" ]] # if the unpaired file exists
  then
    rm ${fastqDir}/${line}.fastq # remove it
  fi # otherwise do nothing
done

##QC: SRR2397083 doesn't exist as unpaired. SRR2859127 does exist as unpaired

## delete unpaired Signor fastq files
cd /scratch/cat7ep/simCline/biosampleresults
fastqDir=/scratch/cat7ep/fasterq

#get the SRS numbers from the Machado files only
signor=$( grep "Signor" ./concatenated.csv | \
  awk -F"," '{
    split ($0,array,",")
    SRSnum= array[16]
    print SRSnum
  }' )

for line in $signor; do
  if [[ -e "${fastqDir}/${line}.fastq" ]] # if the unpaired file exists
  then
    rm ${fastqDir}/${line}.fastq # remove it
  fi # otherwise do nothing
done
