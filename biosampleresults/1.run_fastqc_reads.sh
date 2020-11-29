#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=20G
#SBATCH --time=72:00:00
#SBATCH --partition=standard
#SBATCH --account=berglandlab


#Load fastqc
module load fastqc

#change to folder containing reads
cd /project/berglandlab/Overwintering_Experiment/Year_2018_2019/usftp21.novogene.com/raw_data

#Make an array of all folder containing reads
files=OW*
#echo ${files}

#Run loop running fastqc
for i in  ${files}
	do
		echo "now processing" ${i}
		fastqc ${i}/*.fq.gz
	done
	
# Now consolidate all files into a single folder
mkdir ../QC_raw_reads
find * -name "*fastqc.zip" -exec cp {} ../QC_raw_reads \;

cd ../QC_raw_reads

=
#Now run multi QC
# I downloaded a multiQC singularity thanks to Cory.
module load singularity
#singularity pull --name multiqc.sif shub://cory-weller/singularity:multiqc
singularity run /home/yey2sn/software/multiqc.sif ./

echo "completed"