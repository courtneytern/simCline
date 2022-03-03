# Functions that will take in variant IDs of interest
# and calculate the allele frequency and effective read depth

# module load gcc/7.1.0  openmpi/3.1.4  R
# R

library(SeqArray)
library(data.table)

# read in the filtered SNP ouput from filter_sites.R
## work in /project for easier collaboration
allSitesInfo<- fread("/project/berglandlab/courtney/simCline/data_files/all_SNPs_stats.txt")
pooledVariants<- allSitesInfo$id.x
individVariants<- allSitesInfo$id.y

# calculate freq and n eff
calcPooled<- function(pooledVariants){
  
}# calcPooled
