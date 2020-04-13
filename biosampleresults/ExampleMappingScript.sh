module load gcc/7.1.0
module load bwa/0.7.17
module load samtools/1.10
module load picard/2.20.6

# define some parameters. Take in SRA accession number and unique identifier
sra=${1}
identifier=$2
# sampleB=($(echo ${sample} | cut -d"_" -f1-8))
# sampID=($(echo ${sample} | cut -d"_" -f1-7))
# pond=($(echo ${sample} | cut -d"_" -f9-12))

# threads=10
echo $sra $identifier
#echo $sampleB
inputDir="/scratch/cat7ep/fastq"
interDir="/scratch/cat7ep/interDir"
outputDir="/scratch/cat7ep/bamfiles"

### trim out nextera and index seq
java -jar /scratch/cat7ep/simCline/biosampleresults/trimmomatic-0.39 PE -phred33 \
        ${inputDir}/${sra}_1.fastq \
        ${inputDir}/${sra}_2.fastq \
        ${interDir}/${sra}_1.P_trimm.fastq \
        ${interDir}/${sra}_1.U_trimm.fastq \
        ${interDir}/${sra}_2.P_trimm.fastq \
        ${interDir}/${sra}_2.U_trimm.fastq \
        ILLUMINACLIP:/scratch/cat7ep/simCline/biosampleresults/trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:8:true


### first, merge overlapping reads from the paired reads
### f= forward reads, r= reverse reads, o=base name for output files, j=threads
~/pear-0.9.11/bin/pear \
        -f ${interDir}/${sra}_1.P_trimm.fastq \
        -r ${interDir}/${sra}_2.P_trimm.fastq \
        -o ${interDir}/${sra} \
      #  -j ${threads}

### next, map to reference genome
bwa mem -R "@RG\tID:${sra}\tSM:${identifier}\PL:illumina" \
        /project/berglandlab/courtney/simCline/refgenomes/combinedref.fasta \
        ${interDir}/${sra}.assembled.fastq | \
samtools view -L /scratch/kbb7sh/genomefiles/D84Agoodscaffstouse.bed -Suh -q 20 -F 0x100 | \
samtools sort -@ ${threads} -o ${interDir}/${sra}.sort.bam
samtools index ${interDir}/${sra}.sort.bam

# ## unassembled reads
bwa mem -t ${threads} -K 100000000 -Y \
        -R "@RG\tID:${sampID};${cell};${lane}\tSM:${pond}\tPL:illumina\tPU:${sampID};${cell};${lane}" \
        /scratch/kbb7sh/genomefiles/totalHiCwithallbestgapclosed.fa \
        ${interDir}/${cell}_${lane}_${sampID}.unassembled.forward.fastq \
        ${interDir}/${cell}_${lane}_${sampID}.unassembled.reverse.fastq | \
samtools view -L /scratch/kbb7sh/genomefiles/D84Agoodscaffstouse.bed -Suh -q 20 -F 0x100 | \
samtools sort -@ ${threads} -o ${interDir}/${cell}_${lane}_${pond}.filt.unassembled.sort.bam
samtools index ${interDir}/${cell}_${lane}_${pond}.filt.unassembled.sort.bam
#
# ## Next, merge assembled and unassembled bam files and mark duplicates
samtools merge ${interDir}/${cell}_${lane}_${pond}.filt.merged.bam \
        ${interDir}/${cell}_${lane}_${pond}.sort.bam \
        ${interDir}/${cell}_${lane}_${pond}.filt.unassembled.sort.bam
samtools index ${interDir}/${cell}_${lane}_${pond}.filt.merged.bam

### next, merge bam files to single bam file
# samtools merge ${interDir}/${pond}_finalmap.bam ${interDir}/*${pond}.filt.merged.bam
# samtools index ${interDir}/${pond}_finalmap.bam

### Mark duplicates
java -jar $EBROOTPICARD/picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=5000 \
REMOVE_DUPLICATES=true \
INPUT=${interDir}/${pond}_finalmap.bam \
OUTPUT=${outputDir}/${pond}_finalmap_mdup.bam METRICS_FILE=${interDir}/${pond}_finalmap_mdups.metrics CREATE_INDEX=true
