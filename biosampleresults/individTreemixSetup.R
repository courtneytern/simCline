# Run R analyses on individual 
## To be run on Rivanna

# module load gcc/7.1.0  openmpi/3.1.4  R
# R 

library(SeqArray)
library(ggplot2)
library(data.table)
library(dplyr)

setwd("/scratch/cat7ep/simCline/biosampleresults/")
# setwd("~/Downloads/GitHub/simCline/biosampleresults")

##########################
### Setup/filtering ######
##########################
individ.gds <- seqOpen("./individ.gds")
metadata<- read.csv("./concatenated.csv")
# keep just individ samps
ind_metadata<- metadata[metadata$p.i=="I",]

#separate evo canyon NFS and SFS 
## reset city to NFS or SFS appropriately
ind_metadata$city<- as.character(ind_metadata$city)
ind_metadata$city[grepl(":NFS",ind_metadata$identifier)]<- "NFS"
ind_metadata$city[grepl(":SFS",ind_metadata$identifier)]<- "SFS"

# make column to identify population
ind_metadata$population<- paste(ind_metadata$country,ind_metadata$city,
                                ind_metadata$year,ind_metadata$month,sep=".")

# make list of each group of population sample ids
popSampsList<- list()
for(pop in unique(ind_metadata$population)){
  currentPop<- ind_metadata[ind_metadata$population==pop,]
  popSampsList[[length(popSampsList)+1]] <- as.character(currentPop$sra)
}# end for
names(popSampsList)<- unique(ind_metadata$population)
popSampsList

#####################
## Write function ###
#####################

####sample 10,000 ids for testing
### REMOVE LATER WHEN TESTED PROPERLY
#samp.ids <- as.numeric(sample(x=seqGetData(individ.gds,"variant.id"), size=10000))

# takes in parameters 'popSamps,' which includes all sample ids from one population
# outputs a list of "AD,RD" for treemix
popTreemix<- function(popSamps){
  seqSetFilter(individ.gds, sample.id=popSamps)
  #TEST: seqSetFilter(individ.gds,variant.id=samp.ids)
  
  dosage.mat<- seqGetData(individ.gds,"$dosage_alt")
  numAlt<- apply(dosage.mat,2, function(x) {sum(x==2,na.rm=T)*2 + sum(x==1,na.rm=T)} )
  numRef<- apply(dosage.mat,2, function(x) {sum(x==0,na.rm=T)*2 + sum(x==1,na.rm=T)} )
  paste(numAlt,numRef,sep=",")
}#end popTreemix
## TEST popTreemix
#popTreemix(popSampsList$USA.Zuma.2012.Feb)

## This needs to be pasted into columns. each column with be headed with the population name
## and have as many rows as SNPs.

###################
## Use function ###
###################
# init matrix with as many rows as snps
treemixTable<- matrix(nrow=length(seqGetData(individ.gds,"variant.id")))
for(samps in popSampsList){
  treemixTable<- cbind(treemixTable,popTreemix(samps))
}# end for
# remove the first column of NAs
treemixTable<- treemixTable[,-1]
colnames(treemixTable)<- unique(ind_metadata$population)
head(treemixTable)

write.csv(treemixTable, file="individTreemixInput.csv", quote=F, row.names=F)
write.table(treemixTable, file="individTreemixInput.txt", sep=" ", quote=F, row.names=F)


#################
## Old drafts ###
#################

# ids<- seqGetData(individ.gds, "sample.id")
# snp.dt <- data.table(chr=seqGetData(individ.gds, "chromosome"),
#                      pos=seqGetData(individ.gds, "position"),
#                      id=rep(seqGetData(individ.gds,"sample.id"),each=length(ids)),
#                      dosage=expand.grid(seqGetData(individ.gds, "$dosage_alt"))
#                      )
# colnames(snp.dt)<- c("chr","pos","id","dosage_alt")

# ## calculate allele freq 
# # get population name by merging with concatenated.csv on sample.id 
# # sum ad divide by num individuals 
# #### presently this only calculates rbinom once and uses that value for all ==1
# snp.dt$dosageFreq<- ifelse(snp.dt$dosage_alt==1,mean(rbinom(2,1,.5)),snp.dt$dosage_alt/2)
# snp.dt
# 
# 
# ##########################
# ### Merge pop names ######
# ##########################
# concatenated<- read.csv("./concatenated.csv")



