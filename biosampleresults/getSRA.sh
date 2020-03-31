## run through concatenated file with awk
## for every remaining SRA number (save as variable from awk'{}' > ), copy all fastq files from fastq to project
## $15 is row O sra $16 is row P identifier


cat /scratch/cat7ep/simCline/biosampleresults/concatenated.csv | \
#SRA={ awk -F"," '{$15}' }
SRA=ERR2864366
cp /scratch/cat7ep/fastq/${SRA}*.fastq /project/berglandlab/courtney/simCline/fastq
