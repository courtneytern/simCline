#####################################
### Script for parsing .txt files ###
#####################################


### ************** biosamp 1 ************** ### Barghi SraRunTable5
##ignore D. mel lines with -v

grep -v -E 'Run|Dmel' SraRunTable5.txt.csv | \
awk -F"," '

BEGIN{
print "author,species,numInd,p/i,country,state,city,lat,long,year,month,day,biosamp,sra,identifier"
}
{
  organism=$43
  individ=1000
  pi="P"
  country="USA"
  state="Florida"
  city="Tallahassee"
  lat=30.4383
  long=84.2807
  year="2010"
  month= "Nov"
  day="NA"
  biosamp=$7
  sra=$1
  identifier=$2

  print "Barghi,"organism","individ","pi","country","state","city","\
    lat","long","year","month","day","biosamp","sra",Barghi:"city":"state":"country":"month":"day":"year":"identifier

} ' > parse1.csv

### ************** biosamp 2 ************** ### Machado
grep -v 'Run' SraRunTable.txt.csv | \
awk -F"," '

BEGIN{
print "author,species,numInd,p/i,country,state,city,lat,long,year,month,day,biosamp,sra,identifier"
}
{
  organism=$33 #col AC
  city="NA"
  country=$22
  split($24,loc,": ")
    state=loc[2]
  split($10,date,"-")
    year=date[3]
    month=date[2]
    day=date[1]
  biosamp=$6
  sra=$1
  identifier=$26

  ##number of individuals varies by the date
    if( index($0,"1-May-2011")!=0 ) {individ=44}
    else if( index($0,"1-Sep-2010")!=0 ) {individ=50}
    else if( index($0,"1-Jul-2010")!=0 ) {individ=65}
    else if( index($0,"1-Aug-2011")!=0 ) {individ=44}
    else if( index($0,"1-Sep-2011")!=0 ) {individ=50}
    else if( index($0,"1-Nov-2011")!=0 ) {individ=41}
    else if( index($0,"1-Oct-2011")!=0 ) {individ=38}
    else {individ="NA"}

   ##latitude and p/i varies by location
   pi="I"
   if( index($0,"Florida")!=0 ) {
      lat=25.5
      long=-80.477
      city="Homestead"
    }
   else if( index($0,"Virginia")!=0 ) {
     lat=38
     long=-78.47
     city="Charlottesville"
   }
   else if( index($0,"Penn")!=0 ) {
       lat=39.88
       long=-75.41
       city="Linvilla"
       if( index($0,"2010 pool")!=0 ) {pi="P"}
     }
   else if( index($0,"Maine")!=0 ) {
     lat=44
     long=-69.96
     city="Bowdoin"
   }
   else {lat="NA"}

   print "Machado,"organism","individ","pi","country","state","city","\
      lat","long","year","month","day","biosamp","sra",Machado:"city":"state":"country":"month":"day":"year":"identifier

} ' > parse2.csv

# ### ************** biosamp 3 ************** ### Kang
# grep -E 'BioSample|Organism|ecotype' biosample_result-3.txt | \
# awk -F";" '
#
# BEGIN{
# print "author,species,p/i,numInd,country,state,city,lat,long,year,month,day,biosamp,sra,identifier"
# }
# {
#  if(NR%3==1) {
#   split($1, biosamp, ": ")
#   split($2, identifier, ": ")
#   split($3, sra, ": ")
#   }
#  else if (NR%3==2){
#   split($1,organism,": ")
#  }
#  else{
#    split($1, ecotype, "=\"")
#    split(ecotype[2],country,", ")
#    print "Kang,"organism[2]",20,P,Israel,"country[1]",Mount Carmel,32.7427 N,35.0484 E,2014,Oct,26,"biosamp[3]","sra[2]","identifier[2]
#  }
#
# } ' > parse3.csv

### ************** biosamp 4 ************** ### Palmieri
## Don't want the RNA-seq data
grep -v -E 'Run|RNA-Seq' SraRunTable6.txt.csv | \
awk -F"," '

BEGIN{
print "author,species,numInd,p/i,country,state,city,lat,long,year,month,day,biosamp,sra,identifier"
}
{
  organism=$29
  individ=43
  pi="P"
  country=$21
  state="NA"
  city="NA"
  lat="NA"
  long="NA"
  year="1998"
  month= "NA"
  day="NA"
  biosamp=$8
  sra=$1
  identifier=$25

print "Palmieri,"organism","individ","pi","country","state","city","\
  lat","long","year","month","day","biosamp","sra",Palmieri:"city":"state":"country":"month":"day":"year":"identifier

} ' > parse4.csv

### ************** biosamp 5 ************** ### Sedghifar
##don't need info from accession, but need to pull the line to allow uniform printing
grep -E 'BioSample|Organism|source|Accession' biosample_result-5.txt | \
#grep -v 'Australia' | \
awk -F";" '

BEGIN{
print "author,species,numInd,p/i,country,state,city,lat,long,year,month,day,biosamp,sra,identifier"
} #end of begin
{
  day="NA"

 if( index($0,"BioSample")!= 0 ) {
  split($1, biosamp, ": ")
  split($2, identifier, ": ")
     gsub(" ","",identifier[2])
  split($3, sra, ": ")

   ##treat Australia vs US country and dates
   country="USA"
   year="2011"
  month="Sep"
  if ( index($0,"Aus")!=0 ){
   country="Australia"
   year="2004"
   month="Feb"
  }

  ##treat locations with missing "location" lines
  lat="NA"
  if ( index(identifier[2],"RhodeIsland")!=0) {city[2]="Rhode_Island"; lat=41.84}
  else if ( index(identifier[2],"Florida")!=0) {city[2]="Florida"; lat=25.47}
  else if ( index(identifier[2],"South_Aus")!=0) {city[2]="Queensland"; lat=-42.77}
  else if ( index(identifier[2],"North_Aus")!=0) {city[2]="Tasmania"; lat=-5.54}
 }

 else if ( index($0,"Organism")!= 0 ){
  split($1,organism,": ")
 }

 else if ( index($0,"source")!=0 ){
  split($1, location, "=\"")
  split(location[2],city,", ")
  sub(/"/,"",city[2])

  ##add latitudes and longitudes, write out full states
  lat="NA";long="NA"
  if ( index($0,"Panama")!=0 ) {lat=8.9824 N; long=-79.5199}
  else if ( index($0,"Savannah")!=0 ) {lat=32.0809; long=-81.0912; city[2]="Georgia"}
  else if ( index($0,"Conway")!=0 ) {lat=33.8360; long=-79.0478; city[2]="South Carolina"}
  else if ( index($0,"Raleigh")!=0 ) {lat=35.7796; long=-78.6382; city[2]="North Carolina"}
  else if ( index($0,"Richmond")!=0 ) {lat=37.5407; long=-77.4360; city[2]="Virginia"}
  else if ( index($0,"Princeton")!=0 ) {lat=40.3573; long=-74.6672; city[2]="New Jersey"}
  else if ( index($0,"Fairfield")!=0 ) {lat=44.5884; long=-69.5986; city[2]="Maine"}
 }

 else {
   print "Sedghifar,"organism[2]",NA,P,"country","city[2]","city[1]","lat","long","year","month","day","biosamp[3]","sra[2]",Sedghifar:"city[1]":"city[2]":"country":"month":"day":"year":"identifier[2]
   city[1]="NA"; city[2]="NA"; sra[2]="NA"; lat="NA";long="NA"
 }
}
  ' | \
grep -v -E ",,|Australia" > parse5.csv

### ************** biosamp 6 ************** ### Signor
grep -E 'BioSample|date|Organism|location' biosample_result-6.txt  | \
awk -F";" '

BEGIN{
print "author,species,numInd,p/i,country,state,city,lat,long,year,month,day,biosamp,sra,identifier"
}
{
  month="Feb"
  day="11"
  year="2012"

 if(NR%4==1) {
  split($1, biosamp, ": ")
  split($2, identifier, ": ")
  split($3, sra, ": ")
  }
 else if(NR%4==2){
  split($1,organism,": ")
 }
 else if(NR%4==3){
  split($1,location,"=\"")
  split(location[2],country,": ")
  split(country[2],state,", ")
  sub(/"/,"",state[2])

   print "Signor,"organism[2]",20,I,"country[1]","state[1]","state[2] \
      ",34.0259,118.7798,"year","month","day","biosamp[3]","sra[2]",Palmieri:"state[2]":"state[1]":"country[1]":"month":"day":"year":"identifier[2]
 }

} ' > parse6.csv

### ************** sra run table 2 ************** ### Barghi...Kofler
grep -v 'Run' SraRunTable2.txt.csv | \
awk -F"," '

BEGIN{
print "author,species,numInd,p/i,country,state,city,lat,long,year,month,day,biosamp,sra,identifier"
}
{
  organism=$22
  individ=1000
  pi="P"
  country="USA"
  state="Florida"
  city="Tallahassee"
  lat=30.4383
  long=-84.2807
  year="2010" #previous barghi paper cited in current, copying date from other barghi data
  month= "Nov"
  day="NA"
  biosamp=$7
  sra=$1
  identifier=$2

   print "BarghiKofler,"organism","individ","pi","country","state","city","\
      lat","long","year","month","day","biosamp","sra",BarghiKofler:"city":"state":"country":"month":"day":"year":"identifier

} ' | \
grep "Base" > parse7.csv #only include base populations

# ### ************** sra run table 3 ************** ### Nouhaud
# grep -v 'Run' SraRunTable3.txt.csv | \
# awk -F"," '
#
# BEGIN{
# print "author,species,numInd,p/i,country,state,city,lat,long,year,month,day,biosamp,sra,identifier"
# }
# {
#   organism=$38
#   individ=50
#   pi="P"
#   country="Madagascar"
#   state="NA"
#   city="NA"
#   lat="NA"
#   long="NA"
#   year="1998"
#   month= "NA"
#   day="NA"
#   biosamp=$7
#   sra=$1
#   identifier=$2
#
#    print "Nouhaud,"organism","individ","pi","country","state","city","\
#       lat","long","year","month","day","biosamp","sra","identifier
#
# } ' > parse8.csv
