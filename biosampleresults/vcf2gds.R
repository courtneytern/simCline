#Convert VCF to GDF and filter 

library(SeqArray)
library(ggplot2)
library(data.table)

#load vcf and gds files 
vcf.file <- "~/Downloads/GitHub/simCline/biosampleresults/pooledData.vcf"
gds.output <- "~/Downloads/GitHub/simCline/biosampleresults/pooledData.gds"
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

table(snp.dt$chr)
prop.table(table(snp.dt$nAlleles))
##keep only where 2 alleles
nAlleles=seqGetData(gds.file, "$num_allele")
seqSetFilter(gds.file, variant.id=snp.dt[nAlleles==2]$id[1:1000])

adList<- seqGetData(gds.file, "annotation/format/AD")
rdList<- seqGetData(gds.file, "annotation/format/RD")


dat <- data.table(ad=expand.grid(adList$data)$Var1,
                  rd=expand.grid(rdList$data)$Var1,
                  population=rep(seqGetData(gds.file, "sample.id"), dim(adList$data)[2]),
                  variant.id=rep(seqGetData(gds.file, "variant.id"), each=dim(adList$data)[1]))
dat[,freqAlt:=ad/(ad+rd)]

dat.ag <- dat[,list(population=population,nmissing=sum(is.na(ad)), aveAD=mean(ad, na.rm=T), freqAlt=sum(ad, na.rm=T)/sum(ad+rd, na.rm=T)), list(variant.id)]

#calculate the unweighted frequency of alt alleles
unweightFreqAlt<-mean(dat.ag$freqAlt)
unweightFreqAlt

ggplot(data=dat, aes(x=freqAlt)) + geom_histogram() + facet_wrap(~population)

#calculate missing rate
freqMissing<-dat.ag$nmissing/sum(dat$ad,dat$rd,dat.ag$nmissing,na.rm=T)
mean(freqMissing)
var(freqMissing)

ggplot(data=dat.ag, aes(x=freqMissing)) + geom_histogram()
