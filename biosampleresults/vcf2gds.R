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

#set up data table
gds.file <- seqOpen(gds.output)
snp.dt <- data.table(chr=seqGetData(gds.file, "chromosome"),
                     pos=seqGetData(gds.file, "position"),
                     nAlleles=seqGetData(gds.file, "$num_allele"),
                     id=seqGetData(gds.file, "variant.id"),
                     seqMissing(gds.file, per.variant=T))
# #filtering where chromosome is not Dmel 
snp.dt <- snp.dt[grepl("Dsim_Scf_2L|Dsim_Scf_2R|Dsim_Scf_3L|Dsim_Scf_3R", chr)]

#go to filtering.R
