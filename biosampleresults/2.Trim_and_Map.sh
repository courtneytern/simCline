#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=12:00:00
#SBATCH --partition=standard
#SBATCH --account=berglandlab
#SBATCH -o /scratch/cat7ep/slurmOut/trimmap-fix2.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/trimmap-fix2.%A_%a.err # Standard error
#SBATCH --array=13,15-18,22-24,28-29,34-39,41-50,54,57,59-63,75,84,92-99,101,102,104,105,108,110,111,113,114,129,138,158,159,164,165,182,185-190,198,199

####### sbatch /scratch/cat7ep/simCline/biosampleresults/2.Trim_and_Map.sh
## 1-577 non-Signor
## 577-759 is just Signor. See 2.x to run Signor (interleaved)

# This script will initiate a pipeline which will do some quality QC on the reads and then will proceed to map the reads to a reference genome.
# Prepared by Joaquin C. B. Nunez, PhD -- Sep 24, 2020
# yey2sn@virginia.edu

# NOTES ON NOMENCLATURE: This script uses nomenclature which can be confusing. The first part of the script split raw reads into insert-"merged"-reads (hereby called merged) and unmerged reads (those which were not merged). As such, all operations done using ether of these reads will have the term "merged" or "unmerged" attached to them. At a later point in the script, I combine bam files using "samtools merge" the output of this combination is a joint-bam file (hereby called "joint"). Thus, the joint suffix referes to this step. Other suffix used here are: "srt" which mean picard sorted, and "rmdp" which mean picard-removed duplicated reads.

#Load Rivanna modules
module load gcc/9.2.0 # for bbmap
module load bwa
module load bbmap
module load fastqc
module load samtools
module load qualimap
module load picard

#Define important file locations
#RAW READS indicates the folder where the raw reads are stored. (As fastq)
RAW_READS=/scratch/cat7ep/fasterq

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/cat7ep/individPipeline/TrimMap

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/project/berglandlab/courtney/simCline/refgenomes/simulans

# This is a file with the name all the samples to be processed. one sample name per line
SAMPLE_FILE=/scratch/cat7ep/simCline/biosampleresults/individFileNames.txt

#This is a unique number id which identifies this run
unique_run_id=`date +%N`

#Name of pipeline
PIPELINE=TrimMap

#Define parameters
CPU=1 # number of cores
QUAL=40 # Quality threshold for samtools
JAVAMEM=18G # Java memory

###########################################################################
###########################################################################
# Determine sample to process, "i"
###########################################################################
###########################################################################

i=`sed -n ${SLURM_ARRAY_TASK_ID}p $SAMPLE_FILE`

###########################################################################
###########################################################################
# Begin Pipeline
###########################################################################
###########################################################################
#This part of the pipeline will generate log files to record warnings and completion status

# Welcome message
echo "your unique run id is" $unique_run_id
date

if [[ -e "${PIPELINE}.warnings.log" ]]
then
	echo "Warning log exists"
else
	echo "Log doesnt exist. lets fix that"
	touch $WORKING_FOLDER/${PIPELINE}.warnings.log
fi

if [[ -e "${PIPELINE}.completion.log" ]]
then
	echo "Warning log exists"
else
	echo "Log doesnt exist. lets fix that"
	touch $WORKING_FOLDER/${PIPELINE}.completion.log
fi


# Move to working directory
cd $WORKING_FOLDER

###########################################################################
###########################################################################
# Generate Folders and files
###########################################################################
###########################################################################
# this part of the script will check and generate, if necessary, all of the output folders used in the script

#Generating new folders
echo "Check if folders have been made"
if [[ -d "merged_reads" ]]
then
	echo "Working merged_reads folder exists"
else
	echo "folder doesnt exist. lets fix that"
	mkdir $WORKING_FOLDER/merged_reads
fi

if [ -d "unmerged_reads" ]
then
	echo "Working unmerged_reads folder exist"
else
	echo "folder doesnt exist. lets fix that"
	mkdir $WORKING_FOLDER/unmerged_reads
fi

if [ -d "mapping_stats" ]
then
	echo "Working mapping_stats folder exist"
else
	echo "folder doesnt exist. lets fix that"
	mkdir $WORKING_FOLDER/mapping_stats
fi

if [ -d "read_stats" ]
then
	echo "Working read_stats folder exist"
else
	echo "folder doesnt exist. lets fix that"
	mkdir $WORKING_FOLDER/read_stats
fi

if [ -d "joint_bams" ]
then
	echo "Working joint_bams folder exist"
else
	echo "folder doesnt exist. lets fix that"
	mkdir $WORKING_FOLDER/joint_bams
fi

if [ -d "joint_bams_qualimap" ]
then
	echo "Working joint_bams_qualimap folder exist"
else
	echo "folder doesnt exist. lets fix that"
	mkdir $WORKING_FOLDER/joint_bams_qualimap
fi

###########################################################################
###########################################################################
# Trim and merge reads
###########################################################################
###########################################################################
# This part of the pipeline will trim and merge the reads. It is very likely that the reads will be split into merged and unmerged.
# Both reads will be mapped. This loop operates using a while-read-do-done structure. the while loop is feed a file "SAMPLE_FILE"
# where  all sample names are stored, one name per line. This can be leveraged for parallelization.

# Setting sample name to user input

#while read $i #${files}
	#do #---- Open Do------ <----

	echo "now merging reads for" ${i}

	mkdir $WORKING_FOLDER/merged_reads/${i}
	mkdir $WORKING_FOLDER/unmerged_reads/${i}

	read1=`echo $RAW_READS/${i}_1.fastq`
	read2=`echo $RAW_READS/${i}_2.fastq`

	echo "read 1: "$read1
	echo "read 2: "$read2

  bbmerge.sh \
	in1=$read1 in2=$read2 \
	out=$WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.fq \
	outu1=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.1.fq \
	outu2=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.2.fq \
	-strict

	#Sanity checks
	if [ -s $WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.fq ]; then
	echo ${i} "all good"
	else
	echo "File is empty -- WARNING ISSUED!"
	echo ${i} "Merged reads is empty! check the script, debug, and rerun" >> $WORKING_FOLDER/${PIPELINE}.warnings.log
	fi

	if [ -s $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.1.fq ]; then
	echo ${i} "all good"
	else
	echo "File is empty -- WARNING ISSUED!"
	echo ${i} "Pair 1 reads is empty! check the script, debug, and rerun" >> $WORKING_FOLDER/${PIPELINE}.warnings.log
	fi

	if [ -s $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.2.fq ]; then
	echo ${i} "all good"
	else
	echo "File is empty -- WARNING ISSUED!"
	echo ${i} "Pair 2 reads is empty! check the script, debug, and rerun" >> $WORKING_FOLDER/${PIPELINE}.warnings.log
	fi

	########################################
	#Now do some light trimming on the reads
	echo ${i} "Trimming merge reads"

	bbduk.sh \
	in=$WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.fq \
	out=$WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.trim.fq \
	ftl=15 ftr=285 qtrim=w trimq=20

	rm  $WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.fq

	fastqc $WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.trim.fq \
	-t 20 \
	-o $WORKING_FOLDER/read_stats

	########################################
	#Now do some light trimming on the reads
	echo ${i} "Trimming unmerged reads"

	bbduk.sh \
	in=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.1.fq \
	in2=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.2.fq \
	out=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.trim.1.fq \
	out2=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.trim.2.fq \
	ftl=15 qtrim=w trimq=20

	rm $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.1.fq
	rm $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.2.fq

	fastqc $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.trim.1.fq \
	-t 20 \
	-o $WORKING_FOLDER/read_stats

	fastqc $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.trim.2.fq \
	-t 20 \
	-o $WORKING_FOLDER/read_stats

###########################################################################
###########################################################################
# Map reads to a reference
###########################################################################
###########################################################################
# this part will map reads to the reference genome. Because the reads are likely split into two groups, this script will loop over both types of reads. After
# reads have been mapped, they will be compressed into bam files, sorted, and duplicates will be removed. I will also conduct an intermediary QC step with Qualimap.
# Because there are inherent QC steps here, I have avoided adding extra "warnings" in the log. Remember to take a look at the qualimap and the flagstat outputs to check for inconsistencies.

	for j in merged unmerged
	do # Begin loop of j

		########################################
#J loop#	# Starting mapping
	echo "I will first map ${j} reads of" ${i}

#J loop#	# I will conduct the mapping with BWA-MEM

	if [[ ${j} == "merged" ]]; then
		echo "seems this is merged data, lets map it"
		bwa mem \
		-M \
		-t $CPU \
		$REFERENCE/dsim-mod.fasta \
		$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.reads.strict.trim.fq \
		> $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.sam

	elif [[ ${j} == "unmerged" ]]; then
		echo "seems this is unmerged data, lets map it using a 1-2 approach"
		bwa mem \
		-M \
		-t $CPU \
		$REFERENCE/dsim-mod.fasta \
		$WORKING_FOLDER/unmerged_reads/${i}/${i}.${j}.reads.trim.1.fq \
		$WORKING_FOLDER/unmerged_reads/${i}/${i}.${j}.reads.trim.2.fq \
		> $WORKING_FOLDER/unmerged_reads/${i}/${i}.${j}.sam

	else
		echo "I cant tell what type of data is this -- WARNING!"
		echo ${i} "Something is wrong at the mapping stage" $(date) \
		  $Project_name.warnings.$unique_run_id.log
	fi

#J loop#	#I will now extract some summary stats
	samtools flagstat \
	--threads $CPU \
	$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.sam \
	> $WORKING_FOLDER/${j}_reads/${i}/${i}.flagstats_raw_${j}.sam.txt

#J loop#	#build bam files
	samtools view \
	-b \
	-q $QUAL \
	--threads $CPU  \
	$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.sam \
	> $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.bam

#J loop#	# Sort with picard
	# Notice that once a file has been sorted it is added the "srt" suffix
	java -Xmx$JAVAMEM \
	-jar $PICARD SortSam \
	I=$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.bam \
	O=$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.srt.bam \
	SO=coordinate \
	VALIDATION_STRINGENCY=SILENT

#J loop# Remove duplicates with picard
	# Notice that once a file has been sorted it is added the "rmdp" suffix
	java -Xmx$JAVAMEM \
	-jar $PICARD MarkDuplicates \
	I=$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.srt.bam \
	O=$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.srt.rmdp.bam \
	M=$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.dupstat.txt \
	VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

#J loop# Lets do QC on the bam file
	qualimap bamqc \
	-bam $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.srt.rmdp.bam  \
	-outdir $WORKING_FOLDER/mapping_stats/Qualimap_${i} \
	--java-mem-size=$JAVAMEM

#J loop#	# Clean intermediate files
	rm $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.sam
	rm $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.bam
	rm $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.srt.bam

#J loop#	# Housekeeping
	mv $WORKING_FOLDER/${j}_reads/${i}/${i}.flagstats_raw_${j}.sam.txt \
	$WORKING_FOLDER/mapping_stats
	mv $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.dupstat.txt \
	$WORKING_FOLDER/mapping_stats

	mv $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.reads.strict.trim_fastqc.html \
	$WORKING_FOLDER/read_stats
	mv $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.reads.strict.trim_fastqc.zip \
	$WORKING_FOLDER/read_stats

#J loop#
	done # End loop of j

###########################################################################
###########################################################################
# Merge and asses the final file
###########################################################################
###########################################################################
# Here I will merge the bam outputs from the merge and unmerged portions of the pipeline. Subsequently, I will once again sort and remove duplicated, before performing the final QC on the aligment.

# Merge bams
java -Xmx$JAVAMEM \
 -jar $PICARD MergeSamFiles \
 I=$WORKING_FOLDER/merged_reads/${i}/${i}.merged.srt.rmdp.bam  \
 I=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.srt.rmdp.bam  \
 O=$WORKING_FOLDER/joint_bams/${i}.joint.bam

# Sort merge bams
java -Xmx$JAVAMEM \
 -jar $PICARD SortSam \
 I=$WORKING_FOLDER/joint_bams/${i}.joint.bam \
 O=$WORKING_FOLDER/joint_bams/${i}.joint.srt.bam \
 SO=coordinate \
 VALIDATION_STRINGENCY=SILENT

# Remove duplicates of final file
java -Xmx$JAVAMEM \
 -jar $PICARD MarkDuplicates \
 I=$WORKING_FOLDER/joint_bams/${i}.joint.srt.bam \
 O=$WORKING_FOLDER/joint_bams/${i}.joint.srt.rmdp.bam  \
 M=$WORKING_FOLDER/mapping_stats/${i}.joint.dupstat.txt \
 VALIDATION_STRINGENCY=SILENT \
 REMOVE_DUPLICATES=true

# Assess quality of final file
qualimap bamqc \
 -bam $WORKING_FOLDER/joint_bams/${i}.joint.srt.rmdp.bam   \
 -outdir $WORKING_FOLDER/joint_bams_qualimap/Qualimap_JointBam_${i} \
 --java-mem-size=$JAVAMEM

# Remove intermediary files
rm $WORKING_FOLDER/joint_bams/${i}.joint.bam
rm $WORKING_FOLDER/joint_bams/${i}.joint.srt.bam

###########################################################################
###########################################################################
# Inform that sample is done
###########################################################################
###########################################################################
# This part of the pipeline will produce a notification of completion.

echo ${i} " completed" >> $WORKING_FOLDER/${PIPELINE}.completion.log

echo "pipeline completed" $(date)
