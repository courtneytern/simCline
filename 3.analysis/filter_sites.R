# Make summary table of all SNPs and alt/ref + read depth
## for filtering 

#########
## MAYBE YOU DON'T NEED TO DO ALL THIS
#### Just merge with the treemix output. check to make sure this logic checks out though

# module load gcc/7.1.0  openmpi/3.1.4  R
# R 

library(SeqArray)
library(data.table)
library(dplyr)

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

snp.dt <- data.table(chr=seqGetData(pooled.gds, "chromosome"),
                     pos=seqGetData(pooled.gds, "position"),
                     nAlleles=seqGetData(pooled.gds, "$num_allele"),
                     id=seqGetData(pooled.gds, "variant.id"))
# keep only where 2 alleles
nAlleles=seqGetData(pooled.gds, "$num_allele")
# get ids for previous filters
ids <- snp.dt[nAlleles==2]$id
seqSetFilter(pooled.gds, variant.id=ids)

## Start making table
adList<- seqGetData(pooled.gds, "annotation/format/AD")
rdList<- seqGetData(pooled.gds, "annotation/format/RD")
# compile variables of interest
dat <- data.table(population=rep(seqGetData(pooled.gds, "sample.id"), dim(adList)[2]),
                  chromosome=rep(seqGetData(pooled.gds, "chromosome"), each=dim(adList)[1]),
                  position=rep(seqGetData(pooled.gds, "position"), each=dim(adList)[1]),
                  ref=rep(seqGetData(pooled.gds, "$ref"), each=dim(adList)[1]),
                  alt=rep(seqGetData(pooled.gds, "$alt"), each=dim(adList)[1]),
                  variant.id=rep(seqGetData(pooled.gds, "variant.id"), each=dim(adList)[1]),
                  ad=expand.grid(adList)$Var1,
                  rd=expand.grid(rdList)$Var1
)
# aggregate duplicate Machado rows
agg<- dat %>% group_by(population,variant.id,position,chromosome,ref,alt) %>% summarise_all(sum)
agg<- as.data.table(agg)
# now calc freq alt
agg[,freqAlt:=ad/(ad+rd)]

dat.ag <- dat[,list(nmissing=mean(is.na(ad)), aveAD=mean(ad, na.rm=T), aveRD=mean(rd, na.rm=T),
                    freqAlt=sum(ad, na.rm=T)/sum(ad+rd, na.rm=T),
                    chrom= chromosome[1], pos= position[1]),
              list(variant.id)]
# calc ad/rd per variant id -- merge dat.ag aveAD and aveRD with dat
setkey(agg,variant.id)
setkey(dat.ag,variant.id)
# left outer join. necessitate keeping all of dat
## subset dat.ag to be just the variant id, aveAD, aveRD, chrom, pos
merged <- merge(agg, dat.ag[,c(1,3,4,6,7)], all.x=TRUE)
returnAD <- function(a,r,aa,ar){
  if( is.na(a)|is.na(r) ) {
    return(aa)
  } # if
  else
    return(a)
}
returnRD <- function(a,r,aa,ar){
  if( is.na(a)|is.na(r) ) {
    return(ar)
  } # if
  else
    return(r)
}
merged<- merged[,list( ad=mapply(returnAD,merged$ad,merged$rd,merged$aveAD,merged$aveRD),
       rd=mapply(returnRD,merged$ad,merged$rd,merged$aveAD,merged$aveRD) )
       ]
head(merged)
merged<- merged[,c("population","chrom","pos","ref","alt","ad","rd")]
head(merged)
fwrite(merged,"summary_table_allSites.txt",quote=F,sep="\t")

# col 14 is newCol = AD,RD. need to split that 
# output of this is in /scratch/cat7ep/slurmOut/filterSites.32429916_4294967294.err

# from each GDS, need to get population, chr, pos, ref, alt, 
## seqGetData(gds.file,"$ref") and $alt 
