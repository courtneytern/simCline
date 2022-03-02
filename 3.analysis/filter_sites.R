# Make summary table of all SNPs and alt/ref + read depth
## for filtering 

# module load gcc/7.1.0  openmpi/3.1.4  R
# R 

library(SeqArray)
library(data.table)
library(dplyr)
library(ggplot2)

setwd("/scratch/cat7ep/simCline/data")
individ.gds<- seqOpen("./individ.gds")
pooled.gds<- seqOpen("./pooled.gds")

# get only the SNPs in common
individ.dt <- data.table(chr_pos=seqGetData(individ.gds,"$chrom_pos"),
                         chr=seqGetData(individ.gds, "chromosome"),
                         pos=seqGetData(individ.gds, "position"))
pooled.dt <- data.table(chr_pos=seqGetData(pooled.gds,"$chrom_pos"),
                        chr=seqGetData(pooled.gds, "chromosome"),
                        pos=seqGetData(pooled.gds, "position"))
joined_snps<- intersect(individ.dt,pooled.dt)

# filter to keep only the SNPs in common
seqSetFilterPos(pooled.gds,chr=joined_snps$chr,pos=joined_snps$pos)
seqSetFilterPos(individ.gds,chr=joined_snps$chr,pos=joined_snps$pos)

#############
## Pooled ###
#############
snp.dt <- data.table(chr=seqGetData(pooled.gds, "chromosome"),
                     pos=seqGetData(pooled.gds, "position"),
                     nAlleles=seqGetData(pooled.gds, "$num_allele"),
                     id=seqGetData(pooled.gds, "variant.id"),
                     alleles=seqGetData(pooled.gds,"allele"))
# keep only where 2 alleles
nAlleles=seqGetData(pooled.gds, "$num_allele")
# get ids for previous filters
ids <- snp.dt[nAlleles==2]$id
seqSetFilter(pooled.gds, variant.id=ids)

freqList<- seqGetData(pooled.gds, "annotation/format/FREQ")
# compile variables of interest
dat <- data.table(population=rep(seqGetData(pooled.gds,"sample.id"),dim(freqList)[2]),
                  chromosome=rep(seqGetData(pooled.gds,"chromosome"),each=dim(freqList)[1]),
                  position=rep(seqGetData(pooled.gds,"position"),each=dim(freqList)[1]),
                  freq=expand.grid(freqList)$Var1
)
# FREQ is given as x%. convert to just x and then make numeric
dat$freq<-as.numeric(gsub("%","",as.character(dat$freq)))
# Calculate average frequency per SNP
groupedDat<- dat %>% group_by(chromosome,position)
aveFreqPerSNP<- groupedDat %>% summarize(mean = mean(freq,na.rm=T))
aveFreqPerSNP$mean<- aveFreqPerSNP$mean/100 # divide by 100 to convert pct to decimal
aveFreqPerSNP # table of avg frequency by chr and pos

missingRatePerSNP<- groupedDat %>% summarize(length(freq[is.na(freq)])/length(freq))
colnames(missingRatePerSNP)<- c("chr","pos","missingRate")
missingRatePerSNP # table of missing rate by chr and pos

pooledInfo<- merge(aveFreqPerSNP,missingRatePerSNP,by.x=c("chromosome","position"),by.y=c("chr","pos"))
# merge in variant.id, alleles, and read depth. keep chr pos and read depth cols from dat
## multiply read depth by 23 (# pooled pops) because ADP gives avg DP per sample
pooledRD<- data.table(chromosome=seqGetData(pooled.gds,"chromosome"),
                      position=seqGetData(pooled.gds,"position"),
                      id=seqGetData(pooled.gds,"variant.id"),
                      alleles=seqGetData(pooled.gds,"allele"),
                      readDepth=seqGetData(pooled.gds,"annotation/info/ADP")*dim(freqList)[1]
                      ) 
pooledInfo<- merge(pooledRD,pooledInfo,by=c("chromosome","position"))
colnames(pooledInfo)<- c("chromosome","position","id","pooled.alleles","pooled.readDepth","pooled.AF","pooled.missing")

## SAVE pooledInfo for quick loading later
fwrite(pooledInfo,"/scratch/cat7ep/simCline/data/pooledGDS_SNP_meta.txt",
       quote=F,sep=" ")
# pooledInfo<- fread("/scratch/cat7ep/simCline/data/pooledGDS_SNP_meta.txt")

##############
## individ ###
##############
individInfo <- data.table(chromosome=seqGetData(individ.gds, "chromosome"),
                  position=seqGetData(individ.gds, "position"),
                  id=seqGetData(individ.gds,"variant.id"),
                  alleles=seqGetData(individ.gds,"allele"),
                  freq=seqAlleleFreq(individ.gds),
                  missing=seqMissing(individ.gds),
                  readDepth=seqGetData(individ.gds,"annotation/info/DP")
)
colnames(individInfo)<- c("chromosome","position","id","individ.alleles","individ.AF","individ.missing","individ.readDepth")

## SAVE individInfo for quick loading later 
fwrite(individInfo,"/scratch/cat7ep/simCline/data/individGDS_SNP_meta.txt",
       quote=F,sep=" ")
# individInfo<- fread("/scratch/cat7ep/simCline/data/individGDS_SNP_meta.txt")


## MERGE 
allCombined<- merge(pooledInfo,individInfo,by=c("chromosome","position"))

#############
## FILTER ###
#############
# remove X chromosome, filter for AF and missing thresholds
filtered<- allCombined[ (allCombined$chromosome!="Dsim_Scf_X") &
                        (allCombined$pooled.AF>0.05) & (allCombined$pooled.AF<0.95) &
                        (allCombined$individ.AF>0.05) & (allCombined$individ.AF<0.95) &
                        (allCombined$pooled.missing<0.5) &
                        (allCombined$individ.missing<0.75)
                      ,] # keeps 2075919 sites
filtered$totalReadDepth<- filtered$pooled.readDepth + filtered$individ.readDepth

filtered.dt<- data.table(filtered)
filtered.dt<- filtered.dt[ totalReadDepth<5600 ] # keeps 2,072,013 sites

## MAKE filtered_sites.txt for easy loading later
fwrite(filtered.dt,"/scratch/cat7ep/simCline/data/all_SNPs_stats.txt",
       quote=F,sep=" ")
fwrite(filtered.dt,"/project/berglandlab/courtney/simCline/data_files/all_SNPs_stats.txt",
       quote=F,sep=" ")

############
## PLOTS ###
############
## Summary plots pre-filtering
pdf(file="/scratch/cat7ep/simCline/figures/SNP_summary_plots.pdf")
x<- ggplot(data=aveFreqPerSNP, aes(x=mean)) + 
  labs(title='Average Frequency Per SNP (pooled)') + geom_histogram()
x
y<- ggplot(data=missingRatePerSNP, aes(x=missingRate)) + 
  labs(title='Missing Rate Per SNP (pooled)') + geom_histogram()
y
w<- ggplot(data=dat2, aes(x=freq)) + 
  labs(title='Average Frequency Per SNP (individual)') + geom_histogram()
w
z<- ggplot(data=dat2, aes(x=missing)) + 
  labs(title='Missing Rate Per SNP (individual)') + geom_histogram()
z
dev.off()

## Total read depth post-filtering
pdf(file="/scratch/cat7ep/simCline/figures/totalReadDepth_plotVLINE.pdf")
r<- ggplot(data=filtered,aes(x=totalReadDepth)) +
  labs(title='Total Read Depth per SNP') + geom_density() +
  geom_vline(xintercept = 5600, color="red")
r
dev.off()
  quantile(filtered$totalReadDepth)
  # 0%   25%   50%   75%  100% 
  # 733  3327  3704  4056 79136 