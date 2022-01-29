#!/bin/sh
## TO BE RUN ON RIVANNA

module load samtools
#testing with UA_Cho_15_65

###This script will run through each of the dest_mapped European mel samples,
### then identify which ones have a higher contamination of sim data

cd /project/berglandlab/dest_mapped/
mkdir /scratch/cat7ep/interDirEuro

for dir in *_*; do
  #get only the 2L, 2R, 3L, 3R columns from mel.bam and sim.bam
  samtools idxstats ./$dir/mel.bam | \
    grep -E "2L|2R|3L|3R" | grep -vE "mapped|sim" \
    > /scratch/cat7ep/interDirEuro/"$dir".mel.txt

  samtools idxstats ./$dir/sim.bam | \
    grep -E "sim_2L|sim_2R|sim_3L|sim_3R" | grep -v "het" \
    > /scratch/cat7ep/interDirEuro/"$dir".sim.txt

  #sum up the mel and sim
  mel=$( cat /scratch/cat7ep/interDirEuro/"$dir".mel.txt | awk -F '\t' '{sum+=$3} END {print sum}' )
  sim=$( cat /scratch/cat7ep/interDirEuro/"$dir".sim.txt | awk -F '\t' '{sum+=$3} END {print sum}' )
  total="$(($mel+$sim))"
  #find proportion
  prop="$(echo $sim $total | awk '{ print $1/$2 }')"

  echo "$dir $prop" >> /scratch/cat7ep/simCline/metadata/propsOut.txt
done

################
cd /scratch/cat7ep/simCline/metadata/
#only get samps with proportions >10%
cat ./propsOut.txt | awk -F " " '{
 if ($2*100>10)
  print $1 " " $2*100
}' > propsOutEdited.txt

#find the num simulans
sampName=$( cat ./propsOut.txt | awk -F " " '{print $1}' )
prop=$( cat ./propsOut.txt | awk -F " " '{print $2}' )

echo "Samp,Prop,nFlies,nSimulans" > propsOutEdited.csv
cat ./propsOut.txt | \
while read line; do
 echo $line > tempLine.txt
 sampName=$( awk -F " " '{print $1}' tempLine.txt )
 prop=$( awk -F " " '{print $2}' tempLine.txt )
 nFlies=$( grep -a $sampName ./euro_meta_pre/euroMetadata.csv | awk -F "," '{print $12}' )

 nSim=$( echo $prop $nFlies | awk '{ print $1*$2 }' )
 echo "$sampName,$prop,$nFlies,$nSim" >> propsOutEdited.csv
done

#only get those with over 5 count
cat ./propsOutEdited.csv | awk -F "," '{
 if ($4>5)
  print $0
}' > over5sim.csv

###############
#now move the fastq files from scratch to fastqEuro
grep -v "Samp" ./over5sim.csv | \
while read line; do
  echo $line > tempLine.txt
  samp=$( awk -F "," '{ print $1 }' tempLine.txt)
  mv /scratch/cat7ep/fastq/$samp.fastq /scratch/cat7ep/fastqEuro/$samp.fastq
done

#split fastq into separate files per read
cd /scratch/cat7ep/fastqEuro/
grep -v "Samp" /scratch/cat7ep/simCline/metadata/over5sim.csv | \
while read line; do
  echo $line > tempLine.txt
  samp=$( awk -F "," '{ print $1 }' tempLine.txt)
  cat $samp.fastq | grep '^@.*/1$' -A 3 --no-group-separator > /scratch/cat7ep/fastqEuro/"$samp"_1.fastq
  cat $samp.fastq | grep '^@.*/2$' -A 3 --no-group-separator > /scratch/cat7ep/fastqEuro/"$samp"_2.fastq
done

#create metadata table with same formatting at concatenated.csv
#****EDIT LOCALLY*****
cd /scratch/cat7ep/simCline/metadata
grep -v "Samp" ./over5sim.csv | \
while read line; do
  echo $line > tempLine.txt
  samp=$( awk -F "," '{ print $1 }' tempLine.txt )
  grep "$samp" ./euro_meta_pre/euroMetadata.csv >> ./euro_meta_pre/euroMetadataFinal.csv
done

#parse metadata for European samples
cd /scratch/cat7ep/simCline/metadata/euro_meta_pre
{
  echo "row,author,species,numInd,p/i,continent,country,state,city,lat,long,year,month,day,biosamp,sra,identifier"
  cat euroMetadataFinal.csv | \

    awk -F "," '{
     split($4,date,"/")
     month=date[1]
       if(date[1]==8) {month="Aug"}
       else if (date[1]==9) {month="Sep"}
       else if (date[1]==10) {month="Oct"}
     day=date[2]
      if(index($0,"Trentino")!=0) {state="Trentino"; city="San Michele"}
      else if(index($0,"Purullena")!=0) {state="Purullena"; city="El Bejarin"}
      else {state="NA"; city=$3}

     print NR","$11",Drosophila_simulans,"$12",P,"$10","$2","state","city","$5","$6","$15","month","day","$14","$13","$1
   }'
 } > ../concatenatedEuro.csv
