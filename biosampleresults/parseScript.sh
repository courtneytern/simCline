#####################################
### Script for parsing .txt files ###
#####################################


### Barghi ###
##ignore D. mel lines
grep -v -E 'Run|Dmel' SraRunTable_Barghi.txt.csv | \
awk -F"," '

# BEGIN{
# print "author,species,numInd,p/i,country,state,city,lat,long,year,month,day,biosamp,sra,identifier"
# }
{
  organism=$43
  individ=1000
  pi="P"
  continent="North America"
  country="USA"
  state="Florida"
  city="Tallahassee"
  lat=30.4383
  long=-84.2807
  year="2010"
  month= "Nov"
  day="NA"
  biosamp=$7
  sra=$1
  identifier=$2

  print "Barghi,"organism","individ","pi","continent","country","state","city","\
    lat","long","year","month","day","biosamp","sra",Barghi:"city":"state":"continent":"month":"day":"year":"identifier

} ' | \
grep "Base" | grep "ERR1980726" > parseBarghi.csv #only want one sample of Base population

### Machado ###
grep -v 'Run' SraRunTable_Machado.txt.csv | \
awk -F"," '

# BEGIN{
# print "author,species,numInd,p/i,country,state,city,lat,long,year,month,day,biosamp,sra,identifier"
# }
{
  organism=$33 #col AC
  pi="I"
  city="NA"
  continent="North America"
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

   print "Machado,"organism","individ","pi","continent","country","state","city","\
      lat","long","year","month","day","biosamp","sra",Machado:"city":"state":"continent":"month":"day":"year":"identifier

} ' > parseMachado.csv

### Palmieri ###
## Don't want the RNA-seq data
grep -v -E 'Run|RNA-Seq' SraRunTable_Palmieri.txt.csv | \
awk -F"," '

# BEGIN{
# print "author,species,numInd,p/i,country,state,city,lat,long,year,month,day,biosamp,sra,identifier"
# }
{
  organism=$29
  individ=43
  pi="P"
  continent="Africa"
  country=$21
  state="NA"
  city="Antsirabe"
  lat=-19.8730
  long=47.0291
  year="1998"
  month= "Mar"
  day="19"
  biosamp=$8
  sra=$1
  identifier=$25

print "Palmieri,"organism","individ","pi","continent","country","state","city","\
  lat","long","year","month","day","biosamp","sra",Palmieri:"city":"state":"continent":"month":"day":"year":"identifier

} ' > parsePalmieri.csv

### Sedghifar ###
#don't need Australia info
grep -v -E 'Run|Aus' SraRunTable_Sedghifar.txt.csv | \
awk -F"," '

# BEGIN{
# print "author,species,numInd,p/i,country,state,city,lat,long,year,month,day,biosamp,sra,identifier"
# }
{
  organism=$20
  pi="P"
  continent="North America"
  country="USA"
  year=2011
  month="Sep"
  day="NA"
  biosamp=$6
  sra=$1
  identifier=$28

  if ( index($0,"Panama")!=0 ) {lat=8.60 N; long=-79.3;
    state="Panama"; city="Panama City";
    individ=18; year=2012; month="Jan"}
  else if ( index($0,"Savannah")!=0 ) {lat=32.05; long=-81.11;
    state="Georgia"; city="Savannah";
    individ=45}
  else if ( index($0,"Conway")!=0 ) {lat=33.89; long=-78.95;
    state="South Carolina"; city="Conway";
    individ=50}
  else if ( index($0,"Raleigh")!=0 ) {lat=35.50; long=-78.38;
    state="North Carolina"; city= "Raleigh";
    individ=50}
  else if ( index($0,"Richmond")!=0 ) {lat=37.57; long=-77.40;
    state="Virginia"; city="Richmod";
    individ=25}
  else if ( index($0,"Princeton")!=0 ) {lat=40.33; long=-74.73;
    state="New Jersey"; city="Princeton";
    individ=34}
  else if ( index($0,"Fairfield")!=0 ) {lat=44.68; long=-69.69;
    state="Maine"; city="Fairfield";
    individ=50}
  else if ( index($0,"Rhode Island")!=0 ) {lat=41.84; long=-71.40;
    state="Rhode_Island"; city="Providence";
    individ=33}
  else if ( index($0,"Florida")!=0 ) {lat=25.74; long=-80.42;
    state="Florida"; city="Miami";
    individ=33}

 print "Sedghifar,"organism","individ","pi","continent","country","state","city","\
   lat","long","year","month","day","biosamp","sra",Sedghifar:"city":"state":"continent":"month":"day":"year":"identifier

} ' > parseSedghifar.csv

### Signor ###
grep -v 'Run' SraRunTable_Signor.txt.csv  | \
awk -F"," '

# BEGIN{
# print "author,species,numInd,p/i,country,state,city,lat,long,year,month,day,biosamp,sra,identifier"
# }
{
  organism=$36
  individ=20
  pi="I"
  continent="North America"
  country=$26
  state="California"
  city="Zuma"
  lat=34.0259
  long=-118.7798
  year=2012
  month="Feb"
  day="11"
  biosamp=$7
  sra=$1
  identifier=$39


 print "Signor,"organism","individ","pi","continent","country","state","city","\
   lat","long","year","month","day","biosamp","sra",Signor:"city":"state":"continent":"month":"day":"year":"identifier

} ' > parseSignor.csv

### Kang ###
grep -v 'Run' SraRunTable_Kang.txt.csv  | \
awk -F"," '

# BEGIN{
# print "author,species,numInd,p/i,country,state,city,lat,long,year,month,day,biosamp,sra,identifier"
# }
{
  organism=$28
  individ=20
  pi="P"
  continent="Asia"
  country="Israel"
  state="Mount Carmel"
  city="Nahal Oren"
  lat=32.43
  long=34.58
  year=2014
  month="Oct"
  day="26"
  biosamp=$6
  sra=$1
  identifier=$31


 print "Kang,"organism","individ","pi","continent","country","state","city","\
   lat","long","year","month","day","biosamp","sra",Kang:"city":"state":"continent":"month":"day":"year":"identifier

} ' > parseKang.csv

### Jackson ###
grep -v 'Run' SraRunTable_Jackson.txt.csv  | \
awk -F"," '

# BEGIN{
# print "author,species,numInd,p/i,country,state,city,lat,long,year,month,day,biosamp,sra,identifier"
# }
{
  organism="Drosophila simulans"
  individ=25
  pi="I"
  continent="Africa"

  if(index($0,"Kenyan")!=0){
    country="Kenya"
    state="NA"
    city="Nairobi"
    lat=-1.2921
    long=36.8219
    year=2006
  }
  else {
    country="Madagascar"
    state="NA"
    city="Joffreville"
    lat=-12.4938
    long=49.2056
    year=2002
  }
  month="NA"
  day="NA"
  biosamp=$7
  sra=$1
  identifier=$43


 print "Jackson,"organism","individ","pi","continent","country","state","city","\
   lat","long","year","month","day","biosamp","sra",Jackson:"city":"state":"continent":"month":"day":"year":"identifier

} ' > parseJackson.csv
