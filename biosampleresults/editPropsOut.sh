#!/bin/bash
cd ~/Downloads/GitHub/simCline/biosampleresults

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
 nFlies=$( grep -a $sampName ./euroMetadata.csv | awk -F "," '{print $12}' )

 nSim=$( echo $prop $nFlies | awk '{ print $1*$2 }' )
 echo "$sampName,$prop,$nFlies,$nSim" >> propsOutEdited.csv
done
