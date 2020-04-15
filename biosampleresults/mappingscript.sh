module load gcc/7.1.0
module load bwa/0.7.17
module load samtools/1.10
module load picard/2.20.6
module load gatk/4.0.0.0

# define some parameters. Take in SRA accession number and unique identifier
sra=${1}
identifier=${2}

#Create directories
inputDir="/project/berglandlab/courtney/simCline/fastq"
interDir="/scratch/cat7ep/interDir"
outputDir="/project/berglandlab/courtney/simCline/bamfiles"

#AdapterRemoval path
AdapterRemoval="/project/berglandlab/courtney/adapterremoval-2.3.1/build/AdapterRemoval"
#set current working directory to the intermediate directory
cd ${interDir}

#if paired, trim and merge
if [ -f ${inputDir}/${sra}_1.fastq ] && [ -f ${inputDir}/${sra}_2.fastq ]; then
  {
    ${AdapterRemoval} --file1 ${inputDir}/${sra}_1.fastq --file2 ${inputDir}/${sra}_2.fastq \
                      --basename ${sra}_paired --trimns --trimqualities --collapse
  }
#if unpaired, trim
else if [ -f ${inputDir}/${sra}_1.fastq ]; then
  {
    ${AdapterRemoval} --file1 ${inputDir}/${sra}_1.fastq --basename ${sra}_unpaired --trimns --trimqualities
  }

#map to reference genome
#paired reads
  bwa mem -R "@RG\tID:${sra}\tSM:${identifier}\PL:illumina" \
          /project/berglandlab/courtney/simCline/refgenomes/combinedref.fasta \
          ${interDir}/${sra}_paired.collapsed.truncated | \
  samtools view -uh -q 20 -F 0x100 | \
  samtools sort -o ${interDir}/${sra}.paired.sort.bam
  samtools index ${interDir}/${sra}.paired.sort.bam

#unpaired reads
#need to include more files?
  bwa mem -R "@RG\tID:${sra}\tSM:${identifier}\PL:illumina" \
          /project/berglandlab/courtney/simCline/refgenomes/combinedref.fasta \
          ${interDir}/${sra}_unpaired.truncated | \
  samtools view -uh -q 20 -F 0x100 | \
  samtools sort -o ${interDir}/${sra}.unpaired.sort.bam
  samtools index ${interDir}/${sra}.unpaired.sort.bam

#merge bam files
  samtools merge ${interDir}/${sra}.mergedbam.bam \
           ${interDir}/${sra}.paired.sort.bam \
           ${interDir}/${sra}.unpaired.sort.bam
  samtools index ${interDir}/${sra}.mergedbam.bam

#remove duplicates
java -jar $EBROOTPICARD/picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=5000 \
     INPUT=${interDir}/${sra}.mergedbam.bam \
     OUTPUT=${outputDir}/${sra}.finalmap.bam \
     METRICS_FILE=${interDir}/${sra}.metrics \
     COMMENT= ${identifier} \
     REMOVE_DUPLICATES=true \
     CREATE_INDEX=true

# #indel realignment
# gatk HaplotypeCaller \
#     -i ${outputDir}/${sra}.finalmap.bam \
#     -o ${outputDir}/${sra}.removeindels.vcf \
#     -r /project/berglandlab/courtney/simCline/refgenomes/combinedref.fasta
