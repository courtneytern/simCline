# Functions that will take in variant IDs of interest
# and calculate the allele frequency and effective read depth

# module load gcc/7.1.0  openmpi/3.1.4  R
# R

library(SeqArray)
library(data.table)
library(tidyr)
library(dplyr)

setwd("/project/berglandlab/courtney/simCline/data_files")

# read in the filtered SNP ouput from filter_sites.R
# individ meta: individGDS_SNP_meta.txt // pooled meta: pooledGDS_SNP_meta.txt
# combined: all_SNPs_stats.txt
allSitesInfo<- fread("./all_SNPs_stats.txt")
pooledVariants<- allSitesInfo$id.x
individVariants<- allSitesInfo$id.y
## read in GDS files
pooled.gds<- seqOpen("./pooled.gds")
individ.gds<- seqOpen("./individ.gds")
## metadata including nflies
metadata<- fread("./concatenated.csv")

# calculate freq and n eff
# look at DEST_freeze1 if necessary for help
calcPooled<- function(pooledVariants){
  # keep only the relevant var ids 
  print("Starting function")
  seqSetFilter(pooled.gds,variant.id= pooledVariants)
  
  # make table with ad,rd,nflies,dp
  adList<- seqGetData(pooled.gds, "annotation/format/AD"); print("AD done")
  rdList<- seqGetData(pooled.gds, "annotation/format/RD"); print("RD done")
  dpList<- seqGetData(pooled.gds, "annotation/format/DP"); print("DP done")
  print("Making dat...")
  dat <- data.table(population=rep(seqGetData(pooled.gds, "sample.id"), dim(adList)[2]),
                    variant.id=rep(seqGetData(pooled.gds, "variant.id"), each=dim(adList)[1]),
                    ad=expand.grid(adList)$Var1,
                    rd=expand.grid(rdList)$Var1,
                    dp=expand.grid(dpList)$Var1
  ); print("done!")
  # get nflies for each population from metadata
  nflies<- metadata[metadata$"p/i"=="P",c("identifier","numInd")]
  nflies<- distinct(nflies) #remove machado duplicates
  #merge in nflies info
  print("Merging nflies...")
  dat2<- merge(dat,nflies,by.x="population",by.y="identifier",all.x=T)
  # calculate neff per pool per site
  dat2[,nEff:=round((2*numInd*dp)/(2*numInd + dp))]
  dat2[,af_nEff:=round((ad/dp) * nEff)/nEff]
  # then corrected allele counts for each site across all samples
  dat2[,AD_nEff:=af_nEff*nEff]
  dat2[,RD_nEff:=(1-af_nEff)*nEff]
  print("All done. Returning.")
  dat2
}# calcPooled
pooledFreqNeff<- calcPooled(pooledVariants)
