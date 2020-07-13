#unzip all files from the Jackson fastq folder and move them to one location

while read line
do
{
   cd /scratch/cat7ep/fastqJackson/PRJEB7673/$line
   gunzip ./*.gz
   mv ./*.fastq /scratch/cat7ep/fastqJackson
}
done <  /scratch/cat7ep/simCline/biosampleresults/jacksonAccessionNos.txt
