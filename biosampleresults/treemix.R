## This script will execute TreeMix for pooled and individual samples 
## To be run on Rivanna

# module load gcc/7.1.0  openmpi/3.1.4  intel/18.0  intelmpi/18.0
# module load R 
# R 

library(SeqArray)
library(data.table)

setwd("/scratch/cat7ep/simCline/biosampleresults/")
# setwd("~/Downloads/GitHub/simCline/biosampleresults")

source("./treemix-1.13/src/plotting_funcs.R")
plot_tree("pooled_tree")


# 
# ##########################
# ### Setup/Filtering ######
# ##########################
# 
# # get gds files 
# output.pooled <- "./pooled.gds"
# pooled.gds <- seqOpen(output.pooled)
# output.individ <- "./individ.gds"
# individ.gds <- seqOpen(output.individ)
# 
# ## Filtering for individ gds 
# pooled.dt <- data.table(chr=seqGetData(pooled.gds, "chromosome"),
#                      pos=seqGetData(pooled.gds, "position"),
#                      nAlleles=seqGetData(pooled.gds, "$num_allele"),
#                      varID=seqGetData(pooled.gds, "variant.id"),
#                      seqMissing(pooled.gds, per.variant=T))
# 
# # keep only Dsim chromosomes
# dsim.pooled <- pooled.dt[grepl("Dsim_Scf_2L|Dsim_Scf_2R|Dsim_Scf_3L|Dsim_Scf_3R", chr)]
# # keep only where 2 alleles and d. sim
# nAlleles <- seqGetData(pooled.gds, "$num_allele")
# biallelic.dsim <- dsim.pooled[nAlleles==2]$varID
# 
# ##########################
# ### Function #############
# ##########################
# 
# # get populations
# pooled.pops <- seqGetData(pooled.gds, "sample.id")
# individ.pops <- seqGetData(individ.gds, "sample.id")
# 
# # Run through each population one at a time 
# # get genotypes for the individuals in the population 
# # count ref/alt alleles 
# # output column to be pasted into treemix input 
# #  genome-wide matrix of genotypes for population (?)
# 
# # input one population at a time 
# getTreeMixCol <- function(gds.file, pop){
#   # Set filter to biallelic D. sim 
#   seqSetFilter(pooled.gds, variant.id=biallelic.dsim)
#   # add filter for specific population 
#   seqSetFilter(pooled.gds, sample.id="Barghi:Tallahassee:Florida:North_America:Nov:NA:2010:Dsim_Fl_Base_1") # seqSetFilter(gds.file,sample.id=pop,action="intersect")
#   
#   seqResetFilter(pooled.gds) # seqSetFilter(gds.file)
# } #getTreeMixCol(pop)
# 
# ##seqGetData(pooled.gds, "genotype")
