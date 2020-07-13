#unzip all files from the Jackson fastq folder and move them to one location

while read line
do
{
   cd /scratch/cat7ep/fastqJackson/PRJEB7673/$line
   gunzip ./$line_1 > /scratch/cat7ep/fastqJackson/$line_1.fastq
   gunzip ./$line_2 > /scratch/cat7ep/fastqJackson/$line_2.fastq
}
done <  /scratch/cat7ep/simCline/biosampleresults/jacksonAccessionNos.txt
