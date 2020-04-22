module load gcc/7.1.0
module load bwa/0.7.17
module load samtools/1.10
module load picard/2.20.6
module load gatk/4.0.0.0

# define some parameters. Take in SRA accession number and unique identifier
#test SRR2396839_1.fastq SRR2396839_2.fastq paired
#### ./mappingscript.sh SRR2396839 TESTPAIR
#test SRS3924975_1.fastq unpaired
sra=${1}
identifier=${2}

echo "SRA= "${sra}
echo "Identifier= "${identifier}

#Create directories
inputDir="/project/berglandlab/courtney/simCline/fastq"
interDir="/scratch/cat7ep/interDir"
outputDir="/project/berglandlab/courtney/simCline/bamfiles"

#AdapterRemoval path
AdapterRemoval="/project/berglandlab/courtney/adapterremoval-2.3.1/build/AdapterRemoval"
#set current working directory to the intermediate directory
cd ${interDir}
paired="null"
echo "changed directory"
#if paired, trim and merge
if [ -f ${inputDir}/${sra}_1.fastq ] && [ -f ${inputDir}/${sra}_2.fastq ]
then
  {
    paired="true"
    echo ${paired}

    ${AdapterRemoval} --file1 ${inputDir}/${sra}_1.fastq --file2 ${inputDir}/${sra}_2.fastq \
                      --basename ${sra}_paired --trimns --trimqualities --collapse

    #paired reads
    bwa mem -R "@RG\tID:${sra}\tSM:${identifier}\PL:illumina" \
            /project/berglandlab/courtney/simCline/refgenomes/combinedref.fasta \
            ${interDir}/${sra}_paired.pair1.truncated | \
    samtools view -uh -q 20 -F 0x100 | \
    samtools sort -o ${interDir}/${sra}.pair1.sort.bam
    samtools index ${interDir}/${sra}.pair1.sort.bam

    bwa mem -R "@RG\tID:${sra}\tSM:${identifier}\PL:illumina" \
            /project/berglandlab/courtney/simCline/refgenomes/combinedref.fasta \
            ${interDir}/${sra}_paired.pair2.truncated | \
    samtools view -uh -q 20 -F 0x100 | \
    samtools sort -o ${interDir}/${sra}.pair2.sort.bam
    samtools index ${interDir}/${sra}.pair2.sort.bam

    bwa mem -R "@RG\tID:${sra}\tSM:${identifier}\PL:illumina" \
            /project/berglandlab/courtney/simCline/refgenomes/combinedref.fasta \
            ${interDir}/${sra}_paired.collapsed.truncated | \
    samtools view -uh -q 20 -F 0x100 | \
    samtools sort -o ${interDir}/${sra}.collapsed.sort.bam
    samtools index ${interDir}/${sra}.collapsed.sort.bam

    #merge collapsed and uncollapsed bam files
      samtools merge ${interDir}/${sra}.mergedbam.bam \
               ${interDir}/${sra}.collapsed.sort.bam \
               ${interDir}/${sra}.pair1.sort.bam \
               ${interDir}/${sra}.pair2.sort.bam
      samtools index ${interDir}/${sra}.mergedbam.bam

      inputBam=${interDir}/${sra}.mergedbam.bam
  }
#if unpaired
elif [ -f ${inputDir}/${sra}_1.fastq ]
then
  {
    paired="false"
    echo ${paired}
    ${AdapterRemoval} --file1 ${inputDir}/${sra}_1.fastq --basename ${sra}_unpaired --trimns --trimqualities

      bwa mem -R "@RG\tID:${sra}\tSM:${identifier}\PL:illumina" \
              /project/berglandlab/courtney/simCline/refgenomes/combinedref.fasta \
              ${interDir}/${sra}_unpaired.truncated | \
      samtools view -uh -q 20 -F 0x100 | \
      samtools sort -o ${interDir}/${sra}.unpaired.sort.bam
      samtools index ${interDir}/${sra}.unpaired.sort.bam

      inputBam=${interDir}/${sra}.unpaired.sort.bam
  }
else{}
fi

echo "Outside of if/else"

#remove duplicates
java -jar $EBROOTPICARD/picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=5000 \
     I=${inputBam} \
     O=${outputDir}/${sra}.finalmap.bam \
     METRICS_FILE=${interDir}/${sra}.metrics \
     COMMENT= ${identifier} \
     REMOVE_DUPLICATES=true \
     CREATE_INDEX=true

# go to indel realignment script
# if ${paired}=="true"; then
#   #open paired indel realignment script
# else if ${paired}=="false"; then
#   #open unpaired indel realignment script
# else
#   echo "ERROR in paired variable"

# gatk HaplotypeCaller \
#     -i ${outputDir}/${sra}.finalmap.bam \
#     -o ${outputDir}/${sra}.removeindels.vcf \
#     -r /project/berglandlab/courtney/simCline/refgenomes/combinedref.fasta
