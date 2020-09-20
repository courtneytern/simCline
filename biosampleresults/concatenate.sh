###########################################
### Script for concatenating .txt files ###
###########################################
#!/bin/bash/

{
 echo "row,author,species,numInd,p/i,continent,country,state,city,lat,long,year,month,day,biosamp,sra,identifier"
 cat parse*.csv concatenatedEuro.csv | \
 #grep -v "author" | \

 awk -F"," '{
 row=$0
 gsub(" ","_",row)
 print NR"," ,row}'
} >concatenated.csv
