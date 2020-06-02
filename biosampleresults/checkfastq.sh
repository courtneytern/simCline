#initial check
id=$( cat /scratch/cat7ep/simCline2err.txt | awk -F" " '{
   split($1,array,"_")
   print array[2]
}' )

for n in $id; do
SRS=$( grep ^$n"," /scratch/cat7ep/simCline/biosampleresults/concatenated.csv | \
  awk -F"," '{
    split ($0,array,",")
    SRSnum= array[16]
    print SRSnum
  }' )

  if [ -f /scratch/cat7ep/fasterq/${SRS}.fastq ] || [ -f /scratch/cat7ep/fasterq/${SRS}_1.fastq ]
  then
    echo "$n ok"
  else
    echo "$n NOT OK"
  fi
done

#only the ones that are not ok
notOK=$( cat /scratch/cat7ep/simCline/biosampleresults/notok.txt | awk -F" " '{
  print $1
}' )

for n in $notOK; do
SRS=$( grep ^$n"," /scratch/cat7ep/simCline/biosampleresults/concatenated.csv | \
  awk -F"," '{
    split ($0,array,",")
    SRSnum= array[16]
    print SRSnum
  }' )

  fasterq-dump ${SRS} -t /scratch/cat7ep/interDir --outdir /scratch/cat7ep/fasterq/
done
