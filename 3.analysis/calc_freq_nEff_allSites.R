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
  # aggregate machado rows
  agg<- dat2 %>% group_by(population,variant.id) %>% summarise_all(sum)
  dat2<- as.data.table(agg)
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
fwrite(pooledFreqNeff,"/scratch/cat7ep/simCline/data/pooled_dosage_table.txt")

#######
# function for individual data. Takes in population and variants of interest.

# SETUP individual populations
ind_metadata<- metadata[metadata$"p/i"=="I",]
#separate evo canyon NFS and SFS
## reset city to NFS or SFS appropriately
ind_metadata$city<- as.character(ind_metadata$city)
ind_metadata$city[grepl(":NFS",ind_metadata$identifier)]<- "NFS"
ind_metadata$city[grepl(":SFS",ind_metadata$identifier)]<- "SFS"
# make column to identify population
ind_metadata$population<- paste(ind_metadata$country,ind_metadata$city,
                                ind_metadata$year,ind_metadata$month,sep=".")
# Get sample IDs per population
popSampsList<- list()
for(pop in unique(ind_metadata$population)){
  currentPop<- ind_metadata[ind_metadata$population==pop,]
  popSampsList[[length(popSampsList)+1]] <- as.character(currentPop$sra)
}# end for
names(popSampsList)<- unique(ind_metadata$population)

#for testing purposes
# popName<- names(popSampsList)[1]

## Will need to run this for each set of popSamps in popSampsList
calcIndivid<- function(popName,individVariants){
  print("***New pop***")
  popTable<- data.table()
  # run in 21 chunks of 100,000 variants
  for(i in 0:20){
    print(paste("i=",i))
    seqResetFilter(individ.gds)
    seqSetFilter(individ.gds, sample.id=popSampsList[[popName]],
                 variant.id=individVariants[(i*100000+1):(i*100000+100000)])

    dosage.mat<- seqGetData(individ.gds,"$dosage")
    print(paste("Starting dat for i=",i))
    dat <- data.table(population=popName,
                      variant.id=rep(seqGetData(individ.gds, "variant.id"),each=length(seqGetData(individ.gds,"sample.id"))),
                      ind.id=rep(seqGetData(individ.gds,"sample.id"),length(seqGetData(individ.gds, "variant.id"))),
                      dosage=expand.grid(dosage.mat)[,1]
    )
    print(paste("Making dat.haflo for i=",i))
    dat.haflo <- dat[,list(haflo=rbinom(1, 1, dosage/2), dosage), list(population, ind.id, variant.id)]
    #table(dat.haflo$haflo, dat.haflo$dosage)
    print(paste("Making dat.haflo.ag for i=",i))
    dat.haflo.ag <- dat.haflo[,list(nRef=sum(haflo==1, na.rm=T), nAlt=sum(haflo==0, na.rm=T)), list(population,variant.id)]
    popTable<- rbind(popTable,dat.haflo.ag)
  }

  popTable
}# calcIndivid

# Run calcIndivid for each population and rbind all together
individTable<- data.table()
for(n in names(popSampsList)){
  print(paste("n=",n))
  individTable<- rbind(individTable, calcIndivid(n,individVariants))
}
fwrite(individTable,"/scratch/cat7ep/simCline/data/individ_dosage_table.txt",row.names = F,
            quote=F,sep=" ")

###################
## MAKE TREEMIX ###
###################
# read in the individ table made from previous step
# just remake pooledFreqNeff from above fxn. should take just a few minutes
ind_dosage<- fread("/scratch/cat7ep/simCline/data/individ_dosage_table.txt")
pooledFreqNeff<- fread("/scratch/cat7ep/simCline/data/pooled_dosage_table.txt")

# format for treemix
ind_dosage[,X:=paste(nRef, nAlt, sep=",")]
pooledFreqNeff[,X:=paste(RD_nEff, RD_nEff+AD_nEff, sep=",")]

# merge in chr and pos info
seqSetFilter(pooled.gds,variant.id= pooledVariants)
dictP <- data.table(variant.id=seqGetData(pooled.gds, "variant.id"),
                  chr=seqGetData(pooled.gds,"chromosome"),
                  pos=seqGetData(pooled.gds,"position")
                  )
seqSetFilter(individ.gds,variant.id=individVariants)
dictI <- data.table(variant.id=seqGetData(individ.gds, "variant.id"),
                   chr=seqGetData(individ.gds,"chromosome"),
                   pos=seqGetData(individ.gds,"position")
                  )
mergedP<- merge(pooledFreqNeff,dictP,by="variant.id",all.x=T)
mergedI<- merge(ind_dosage,dictI,by="variant.id",all.x=T)

#get both to chr,pos,x,pop
m <- rbind(mergedP[,c("population","chr","pos","X")], mergedI[,c("population","chr","pos","X")])
m.wide <- dcast(m, chr+pos ~ population, value.var="X")
fwrite(m.wide,"/scratch/cat7ep/simCline/data/TREEMIX_032522.txt",row.names = F,
       quote=F,sep=" ")

## prep for neis_nj.R
mergedP<- mergedP[,c("population","chr","pos","AD_nEff","RD_nEff")]
mergedI<- mergedI[,c("population","chr","pos","nAlt","nRef")]
colnames(mergedP)<- c("population","chr","pos","ad","rd")
colnames(mergedI)<- c("population","chr","pos","ad","rd")
allSites<- rbind(mergedP,mergedI)
fwrite(allSites,"/scratch/cat7ep/simCline/data/allPops_adrd.txt")

