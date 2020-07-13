#unzip all files from the Jackson fastq folder and move them to one location

for line in /scratch/cat7ep/simCline/biosampleresults/jacksonAccessionNos.txt
do
{
   cd /scratch/cat7ep/fastqJackson/PRJEB7673/{$line}
   gunzip ./{$line}_1 > /scratch/cat7ep/fastqJackson/{$line}_1.fastq
   gunzip ./{$line}_2 > /scratch/cat7ep/fastqJackson/{$line}_2.fastq
}
done
