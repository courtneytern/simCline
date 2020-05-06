SRS=$( grep ^[1-826]"," /scratch/cat7ep/simCline/biosampleresults/concatenated.csv | \
  awk -F"," '{
    split ($0,array,",")
    SRSnum= array[15]
    print SRSnum
  }' )

if [ -f /project/berglandlab/courtney/simCline/bamfiles/*${SRS}* ]
then
{
   echo ""
}

else
{
   echo ${SRS}
}
fi
