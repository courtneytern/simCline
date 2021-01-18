#Convert VCF to GDF and filter 

library(SeqArray)
library(ggplot2)
library(data.table)

#load vcf and gds files 
vcf.file <- "~/Downloads/GitHub/simCline/biosampleresults/pooledData.vcf"
gds.output <- "~/Downloads/GitHub/simCline/biosampleresults/pooled.gds"
#create GDS file 
seqVCF2GDS(vcf.fn = vcf.file, out.fn = gds.output)

#########################################

#go to pooledStatsPlots.R
