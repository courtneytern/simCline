#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=8:00:00
#SBATCH --partition=standard
#SBATCH --account=berglandlab
#SBATCH -o /scratch/cat7ep/slurmOut/trimmap-Signor.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/trimmap-Signor.%A_%a.err # Standard error
#SBATCH --array=577-759

#### sbatch /scratch/cat7ep/simCline/biosampleresults/2.x.Trim_and_Map_interleaved.sh

## This script is adapted from the 2.Trim_and_Map.sh file written by Joaquin C.B. Nunez, adapted to work for interleaved files
## This assumes that you have already successfully run the original 2.Trim_and_Map.sh, and therefore all relevant
## initial directories are already created

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

# Move to working directory
cd $WORKING_FOLDER

###########################################################################
###########################################################################
# Trim reads
###########################################################################
###########################################################################
########################################
#Now do some light trimming on the reads
echo ${i} "Trimming reads"

bbduk.sh \
in=$RAW_READS/${i}.fastq \
out=$WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.trim.fq \
ftl=15 ftr=285 qtrim=w trimq=20

fastqc $WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.trim.fq \
-t 20 \
-o $WORKING_FOLDER/read_stats

###########################################################################
###########################################################################
# Map reads to a reference
###########################################################################
###########################################################################
j=merged

echo "I will now map merged reads of " ${i}

bwa mem \
-M \
-t $CPU \
$REFERENCE/dsim-mod.fasta \
$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.reads.strict.trim.fq \
> $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.sam

#J loop#	#I will now extract some summary stats
	samtools flagstat \
	--threads $CPU \
	$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.sam \
	> $WORKING_FOLDER/${j}_reads/${i}/${i}.flagstats_raw_${j}.sam.txt

#J loop#	#build bam files
	samtools view \
	-b \
	-q $QUAL \
	--threads $CPU \
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
	# qualimap bamqc \
	# -bam $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.srt.rmdp.bam\
	# -outdir $WORKING_FOLDER/mapping_stats/Qualimap_${i} \
	# --java-mem-size=$JAVAMEM

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

  ###########################################################################
  ###########################################################################
  # Copy over to final file, Bam QC
  ###########################################################################
  ###########################################################################
  # Housekeeping step to keep file names consistent with the other files
  # No merging needed for interleaved between merged/unmerged, so just copy over the merged files to the "joint_bams" folder
  # Makes subsequent steps easier with everything in the same location/nomenclature
  cp $WORKING_FOLDER/merged_reads/${i}/${i}.merged.srt.rmdp.bam $WORKING_FOLDER/joint_bams/${i}.joint.srt.rmdp.bam

  # now run QC on the final version of the file
  qualimap bamqc \
   -bam $WORKING_FOLDER/joint_bams/${i}.joint.srt.rmdp.bam   \
   -outdir $WORKING_FOLDER/joint_bams_qualimap/Qualimap_JointBam_${i} \
   --java-mem-size=$JAVAMEM

  ###########################################################################
  ###########################################################################
  # Inform that sample is done
  ###########################################################################
  ###########################################################################
  # This part of the pipeline will produce a notification of completion.

  echo ${i} " completed" >> $WORKING_FOLDER/${PIPELINE}.completion.log

  echo "pipeline completed" $(date)
