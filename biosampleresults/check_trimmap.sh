# this script checks which directories were not properly made after step 2

SAMPLE_FILE=/scratch/cat7ep/simCline/biosampleresults/individFileNames.txt
WORKING_DIRECTORY=/scratch/cat7ep/simCline/biosampleresults
WORKING_FOLDER=/scratch/cat7ep/individPipeline/TrimMap

# run through all individ accession nos
  # if directory exists, do nothing
  # otherwise, get the line number to rerun
while read line
do
  if [[ -e $WORKING_FOLDER/joint_bams_qualimap/Qualimap_JointBam_${line} ]]
  then
    echo ""
  else
    grep -n $line /scratch/cat7ep/simCline/biosampleresults/individFileNames.txt | cut -d: -f1
  fi
done < $SAMPLE_FILE > $WORKING_DIRECTORY/checkTrimmap.txt

# remove empty lines
sed -i '/^$/d' $WORKING_DIRECTORY/checkTrimmap.txt
