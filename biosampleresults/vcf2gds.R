#Convert VCF to GDF and filter 
## Run this script in Rivanna with vcf2gds.sh 
## Comment out one code block at a time, or allocate enough time to run both 

library(SeqArray)
library(data.table)

#########################
### Pooled ##############
#########################

#load vcf and gds files 
vcf.file <- "/scratch/cat7ep/simCline/biosampleresults/pooledData.vcf"
# vcf.file <- "~/Downloads/GitHub/simCline/biosampleresults/pooledData.vcf"
gds.output <- "/scratch/cat7ep/simCline/biosampleresults/pooled.gds"
# gds.output <- "~/Downloads/GitHub/simCline/biosampleresults/pooled.gds"
#create GDS file
seqVCF2GDS(vcf.fn = vcf.file, out.fn = gds.output)

#########################################

#go to pooledStatsPlots.R
# or combinePooledIndividSNPs.R to run both pooled and individ


#########################
### Individ #############
#########################

#load vcf and gds files 
vcf.file <- "/scratch/cat7ep/individPipeline/MergeVCF/Simcline_final_2021.vcf.gz"
gds.output <- "/scratch/cat7ep/simCline/biosampleresults/individ.gds"
#create GDS file 
seqVCF2GDS(vcf.fn = vcf.file, out.fn = gds.output)

#########################################

#go to individTreemixSetup.R
# or combinePooledIndividSNPs.R to run both pooled and individ
