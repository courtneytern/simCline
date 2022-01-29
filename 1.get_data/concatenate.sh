#!/bin/bash/
###########################################
### Script for concatenating .txt files ###
###########################################

cd ./metadata

{
 echo "row,author,species,numInd,p/i,continent,country,state,city,lat,long,year,month,day,biosamp,sra,identifier"
 cat ./parseTables/parse*.csv concatenatedEuro.csv | \
 #grep -v "author" | \

 awk -F"," '{
 row=$0
 gsub(" ","_",row)
 print NR"," ,row}'
} >concatenated.csv
