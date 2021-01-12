#! /bin/bash
#
#SBATCH -J fastqc # A single job name for the array
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=20G
#SBATCH --time=72:00:00
#SBATCH --partition=standard
#SBATCH -o /scratch/cat7ep/slurmOut/fastqc.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/fastqc.%A_%a.err # Standard error
#SBATCH --account=berglandlab


######sbatch /scratch/cat7ep/simCline/biosampleresults/1.run_fastqc_reads.sh
#Load fastqc
module load fastqc

#change to folder containing reads
cd /scratch/cat7ep/fasterq

#Make an array of all folder containing reads
files=*
#echo ${files}
#set output directory
outdir="/scratch/cat7ep/fastqc"

#Run loop running fastqc
for i in  ${files}
	do
		echo "now processing" ${i}
		fastqc -o ${outdir} ${i}
	done

# Now consolidate all files into a single folder
# mkdir ../QC_raw_reads
# find * -name "*fastqc.zip" -exec cp {} ../QC_raw_reads \;
#
# cd ../QC_raw_reads
#
# =
# #Now run multi QC
# # I downloaded a multiQC singularity thanks to Cory.
# module load singularity
# #singularity pull --name multiqc.sif shub://cory-weller/singularity:multiqc
# singularity run /home/yey2sn/software/multiqc.sif ./
#
# echo "completed"
