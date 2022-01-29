#!/bin/sh
#
#SBATCH -J makeVCF # A single job name for the array
#SBATCH --ntasks-per-node=20 #twenty cores
#SBATCH -N 1
#SBATCH -t 96:00:00 ### 96 hrs= 4 days
#SBATCH --mem 0 #take as much as it needs
#SBATCH -o /scratch/cat7ep/slurmOut/makeVCF8.%A_%a.out # Standard output
#SBATCH -e /scratch/cat7ep/slurmOut/makeVCF8.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load samtools/1.10
module load varscan/2.4.4

cd /scratch/cat7ep/simCline/2.mapping_pipeline/
####sbatch /scratch/cat7ep/simCline/2.mapping_pipeline/makeVCF.sh

###Make VCF file with samtools to VarScan pipe
samtools mpileup -f /project/berglandlab/courtney/simCline/refgenomes/combinedref.fasta \
   -b /scratch/cat7ep/simCline/2.mapping_pipeline/inputPooledBam.txt \
   -q 10 -Q 15 | \
java -jar $EBROOTVARSCAN/VarScan.v2.4.4.jar mpileup2snp \
   --output-vcf 1 \
   --vcf-sample-list /scratch/cat7ep/simCline/2.mapping_pipeline/inputPooledSamps.txt \
   --variants - > /project/berglandlab/courtney/simCline/data_files/pooledData.vcf
