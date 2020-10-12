# this generates the freqAlt plots per pooled population

library(SeqArray)
library(ggplot2)
library(data.table)

setwd("~/Downloads/GitHub/simCline/biosampleresults/")

gds.output <- "~/Downloads/GitHub/simCline/biosampleresults/pooled.gds"

gds.file <- seqOpen(gds.output)
snp.dt <- data.table(chr=seqGetData(gds.file, "chromosome"),
                     pos=seqGetData(gds.file, "position"),
                     nAlleles=seqGetData(gds.file, "$num_allele"),
                     id=seqGetData(gds.file, "variant.id"),
                     seqMissing(gds.file, per.variant=T))
# filter where not dmels
snp.dt <- snp.dt[grepl("Dsim_Scf_2L|Dsim_Scf_2R|Dsim_Scf_3L|Dsim_Scf_3R", chr)]

# keep only where 2 alleles
nAlleles=seqGetData(gds.file, "$num_allele")
#sample 100,000 ids
samp.ids <- as.numeric(sample(x=as.character(snp.dt[nAlleles==2]$id), size=10000))
seqSetFilter(gds.file, variant.id=samp.ids)

################
adList<- seqGetData(gds.file, "annotation/format/AD")
rdList<- seqGetData(gds.file, "annotation/format/RD")

# compile variables of interest
dat <- data.table(ad=expand.grid(adList$data)$Var1,
                  rd=expand.grid(rdList$data)$Var1,
                  population=rep(seqGetData(gds.file, "sample.id"), dim(adList$data)[2]),
                  variant.id=rep(seqGetData(gds.file, "variant.id"), each=dim(adList$data)[1]))
dat[,freqAlt:=ad/(ad+rd)]

#functions of dat[] variables
dat.ag <- dat[,list(population=population,nmissing=sum(is.na(ad)), aveAD=mean(ad, na.rm=T), freqAlt=sum(ad, na.rm=T)/sum(ad+rd, na.rm=T)), list(variant.id)]

#don't remember why we wanted to merge 
concatenated<- fread("concatenated.csv")
setnames(dat, "population", "identifier")
setkey(concatenated, identifier)
setkey(dat, identifier)
dat.merged <- merge(dat, concatenated)


#plot frequency of alternate alleles for each of the 42 populations
ggplot(data=dat, aes(x=freqAlt)) + geom_histogram() + facet_wrap(~identifier)

####export gds as vcf
vcf.fn<- "pooledData2.vcf"
seqGDS2VCF(gds.file, vcf.fn, info.var=NULL, fmt.var=NULL, use_Rsamtools=TRUE,
           verbose=TRUE)


