### ************** biosamp 7 ************** ### Bastide
grep -E 'BioSample|Organism|description|Title' biosample_result-7.txt  | \
awk -F";" '

BEGIN{
print "author,species,numInd,p/i,country,state,city,lat,long,year,month,day,biosamp,sra,identifier"
}
{
 if(NR%4==1) {
  split($1, biosamp, ": ")
  split($2, sra, ": ")
  }
 else if(NR%4==2){
  split($1,organism,": ")
 }
 else if(NR%4==3){
  split($1,num,": ")
  sub(/\)"/,"",num[2])
 }
 else{
   split($1,identifier,"=\"")
   sub(/"/,"",identifier[2])
   if(index($0,"Italy")!=0) {
      country="Italy"; city="Bolzano";
      lat="46.4983N"; long="11.3548E";
      year=2011; month="September";
   }
   else{
      country="Austria"; city="Vienna";
      lat="48.2082N"; long="16.3738E";
      year=2010; month="October";
   }


    print "Bastide,"organism[2]","num[2]",P,"country",NA,"city \
       ","lat","long","year","month",NA,"biosamp[3]","sra[2]",Bastide:"city":NA:"country":"month":NA:"year":"identifier[2]
 }

} ' >  bastideParse.csv

{
 echo "row,author,species,numInd,p/i,country,state,city,lat,long,year,month,day,biosamp,sra,identifier"
 cat bastideParse.csv | \
 grep -v "author" | \

 awk -F"," '{
 row=$0
 gsub(" ","_",row)
 print NR"," ,row
  }'
} > bastidecat.csv