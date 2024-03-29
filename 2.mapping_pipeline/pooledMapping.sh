#!/bin/sh
#
#SBATCH -J pooledMapping # A single job name for the array
#SBATCH --ntasks-per-node=20 # twenty cores
#SBATCH -N 1
#SBATCH -t 24:00:00 ### 24 hr
#SBATCH --mem 0
#SBATCH -o /scratch/cat7ep/slurmOut/Jacksonmapping.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/Jacksonmapping.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


###SLURM_ARRAY_TASK_ID=4
#sbatch --array=1-811 /scratch/cat7ep/simCline/2.mapping_pipeline/pooledMapping.sh
module load gcc/7.1.0
module load bwa/0.7.17
module load samtools/1.10
module load picard/2.20.6
module load gatk/4.0.0.0

echo "${SLURM_ARRAY_TASK_ID}"

# define some parameters. Take in SRA accession number and unique identifier
#test SRR2396839_1.fastq SRR2396839_2.fastq paired
#### ./pooledMapping.sh SRR2396839 TESTPAIR
#test SRS3924975_1.fastq unpaired
sra=$( grep ^"${SLURM_ARRAY_TASK_ID}""," /scratch/cat7ep/simCline/metadata/concatenated.csv | \
  awk -F"," '{
    split ($0,array,",")
    SRSnum= array[16]
    print SRSnum
  }' )
identifier=$( grep ^"${SLURM_ARRAY_TASK_ID}""," /scratch/cat7ep/simCline/metadata/concatenated.csv | \
  awk -F"," '{
    split ($0,array,",")
    id= array[17]
    print id
  }' )

echo "SRA= ""${sra}"
echo "Identifier= ""${identifier}"

#Create directories
inputDir="/scratch/cat7ep/simCline/fasterq"
interDir="/scratch/cat7ep/simCline/interDir"
outputDir="/project/berglandlab/courtney/simCline/bamfinal"

#AdapterRemoval path
alias AdapterRemoval="/project/berglandlab/courtney/adapterremoval-2.3.1/build/AdapterRemoval"
#set current working directory to the intermediate directory
cd ${interDir}
paired="null"
echo "changed directory"
#if paired, trim and merge
if [ -f ${inputDir}/"${sra}"_1.fastq ] && [ -f ${inputDir}/"${sra}"_2.fastq ]
then
  {
    paired="true"
    echo ${paired}

    AdapterRemoval --file1 ${inputDir}/"${sra}"_1.fastq --file2 ${inputDir}/"${sra}"_2.fastq \
                      --basename "${sra}"_paired --trimns --trimqualities --collapse

    #paired reads
    bwa mem -R "@RG\tID:${sra}\tSM:${identifier}\PL:illumina" \
            /project/berglandlab/courtney/simCline/refgenomes/combinedref.fasta \
            ${interDir}/"${sra}"_paired.pair1.truncated \
            ${interDir}/"${sra}"_paired.pair2.truncated | \
    samtools view -uh -q 20 -F 0x100 | \
    samtools sort -o ${interDir}/"${sra}".pairs.sort.bam
    samtools index ${interDir}/"${sra}".pairs.sort.bam

    bwa mem -R "@RG\tID:${sra}\tSM:${identifier}\PL:illumina" \
            /project/berglandlab/courtney/simCline/refgenomes/combinedref.fasta \
            ${interDir}/"${sra}"_paired.collapsed.truncated | \
    samtools view -uh -q 20 -F 0x100 | \
    samtools sort -o ${interDir}/"${sra}".collapsed.sort.bam
    samtools index ${interDir}/"${sra}".collapsed.sort.bam

    #merge collapsed and uncollapsed bam files
      samtools merge ${interDir}/"${sra}".mergedbam.bam \
               ${interDir}/"${sra}".collapsed.sort.bam \
               ${interDir}/"${sra}".pairs.sort.bam
      samtools index ${interDir}/"${sra}".mergedbam.bam

      inputBam=${interDir}/${sra}.mergedbam.bam
  }
#if unpaired
else
  {
    paired="false"
    echo ${paired}
    AdapterRemoval --file1 ${inputDir}/"${sra}".fastq --basename "${sra}"_unpaired --trimns --trimqualities

      bwa mem -R "@RG\tID:${sra}\tSM:${identifier}\PL:illumina" \
              /project/berglandlab/courtney/simCline/refgenomes/combinedref.fasta \
              ${interDir}/"${sra}"_unpaired.truncated | \
      samtools view -uh -q 20 -F 0x100 | \
      samtools sort -o ${interDir}/"${sra}".unpaired.sort.bam
      samtools index ${interDir}/"${sra}".unpaired.sort.bam

      inputBam=${interDir}/${sra}.unpaired.sort.bam
  }
fi

echo "Outside of if/else"

#remove duplicates
java -jar "$EBROOTPICARD"/picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=5000 \
     I="${inputBam}" \
     O=${outputDir}/"${sra}".finalmap.bam \
     METRICS_FILE=${interDir}/"${sra}".metrics \
     COMMENT= "${identifier}" \
     REMOVE_DUPLICATES=true \
     CREATE_INDEX=true

## go to makeVCF.sh
