#!/bin/sh
#
#SBATCH -J Euromapping # A single job name for the array
#SBATCH --ntasks-per-node=20 # twenty cores
#SBATCH -N 1
#SBATCH -t 24:00:00 ### 24 hr
#SBATCH --mem 30G
#SBATCH -o /scratch/cat7ep/slurmOut/Euromapping.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/Euromapping.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


###SLURM_ARRAY_TASK_ID=4
#sbatch --array=1-9 /scratch/cat7ep/simCline/biosampleresults/mappingscriptEuro.sh
module load gcc/7.1.0
module load bwa/0.7.17
module load samtools/1.10
module load picard/2.20.6
module load gatk/4.0.0.0

echo "${SLURM_ARRAY_TASK_ID}"

# define some parameters. Take in SRA accession number and unique identifier
sra=$( grep ^"${SLURM_ARRAY_TASK_ID}""," /scratch/cat7ep/simCline/biosampleresults/concatenatedEuro.csv | \
  awk -F"," '{
    split ($0,array,",")
    SRSnum= array[16]
    print SRSnum
  }' )
identifier=$( grep ^"${SLURM_ARRAY_TASK_ID}""," /scratch/cat7ep/simCline/biosampleresults/concatenatedEuro.csv | \
  awk -F"," '{
    split ($0,array,",")
    id= array[17]
    print id
  }' )
echo "SRA= ""${sra}"
echo "Identifier= ""${identifier}"

#Create directories
inputDir="/scratch/cat7ep/fastqEuro"
interDir="/scratch/cat7ep/interDirEuro"
outputDir="/project/berglandlab/courtney/simCline/bamfinal"

#AdapterRemoval path
alias AdapterRemoval="/project/berglandlab/courtney/adapterremoval-2.3.1/build/AdapterRemoval"
#set current working directory to the intermediate directory
cd ${interDir}
paired="null"
echo "changed directory"

#all paired for the European samples
#fastq files come from the fastqCombinePairedEnd output
paired="true"
echo ${paired}

    AdapterRemoval --file1 ${inputDir}/"${identifier}"_1.fastq_pairs_R1.fastq --file2 ${inputDir}/"${identifier}"_2.fastq_pairs_R2.fastq \
                      --basename "${identifier}"_paired --trimns --trimqualities --collapse

    #paired reads
    bwa mem -R "@RG\tID:${sra}\tSM:${identifier}\PL:illumina" \
            /project/berglandlab/courtney/simCline/refgenomes/combinedref.fasta \
            ${interDir}/"${identifier}"_paired.pair1.truncated \
            ${interDir}/"${identifier}"_paired.pair2.truncated | \
    samtools view -uh -q 20 -F 0x100 | \
    samtools sort -o ${interDir}/"${identifier}".pairs.sort.bam
    samtools index ${interDir}/"${identifier}".pairs.sort.bam

    bwa mem -R "@RG\tID:${sra}\tSM:${identifier}\PL:illumina" \
            /project/berglandlab/courtney/simCline/refgenomes/combinedref.fasta \
            ${interDir}/"${identifier}"_paired.collapsed.truncated | \
    samtools view -uh -q 20 -F 0x100 | \
    samtools sort -o ${interDir}/"${identifier}".collapsed.sort.bam
    samtools index ${interDir}/"${identifier}".collapsed.sort.bam

    #merge collapsed and uncollapsed bam files
      samtools merge ${interDir}/"${identifier}".mergedbam.bam \
               ${interDir}/"${identifier}".collapsed.sort.bam \
               ${interDir}/"${identifier}".pairs.sort.bam
      samtools index ${interDir}/"${identifier}".mergedbam.bam

      inputBam=${interDir}/${identifier}.mergedbam.bam

#remove duplicates
java -jar "$EBROOTPICARD"/picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=5000 \
     I="${inputBam}" \
     O=${outputDir}/"${identifier}".finalmap.bam \
     METRICS_FILE=${interDir}/"${identifier}".metrics \
     COMMENT= "${identifier}" \
     REMOVE_DUPLICATES=true \
     CREATE_INDEX=true
