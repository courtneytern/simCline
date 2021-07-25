#Convert VCF to GDF and filter 

library(SeqArray)
library(data.table)

#########################
### Pooled ##############
#########################

#load vcf and gds files 
vcf.file <- "~/Downloads/GitHub/simCline/biosampleresults/pooledData.vcf"
gds.output <- "~/Downloads/GitHub/simCline/biosampleresults/pooled.gds"
#create GDS file 
seqVCF2GDS(vcf.fn = vcf.file, out.fn = gds.output)

#########################################

#go to pooledStatsPlots.R


#########################
### Individ #############
#########################

#load vcf and gds files 
vcf.file <- "~/Downloads/GitHub/simCline/biosampleresults/Simcline_final_2021.vcf.gz"
gds.output <- "~/Downloads/GitHub/simCline/biosampleresults/individ.gds"
#create GDS file 
seqVCF2GDS(vcf.fn = vcf.file, out.fn = gds.output)

#########################################

#go to treemix.R