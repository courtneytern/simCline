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

####
adList<- seqGetData(gds.file, "annotation/format/AD")
rdList<- seqGetData(gds.file, "annotation/format/RD")

# compile variables of interest
dat <- data.table(ad=expand.grid(adList$data)$Var1,
                  rd=expand.grid(rdList$data)$Var1,
                  population=rep(seqGetData(gds.file, "sample.id"), dim(adList$data)[2]),
                  variant.id=rep(seqGetData(gds.file, "variant.id"), each=dim(adList$data)[1]))
dat[,freqAlt:=ad/(ad+rd)]

#######functions of dat[] variables
# dat.ag <- dat[,list(population=population,nmissing=sum(is.na(ad)), aveAD=mean(ad, na.rm=T), freqAlt=sum(ad, na.rm=T)/sum(ad+rd, na.rm=T)), list(variant.id)]
# 
# #don't remember why we wanted to merge 
# concatenated<- fread("concatenated.csv")
# setnames(dat, "population", "identifier")
# setkey(concatenated, identifier)
# setkey(dat, identifier)
# dat.merged <- merge(dat, concatenated)

pdf(file="altAlleles.pdf")
#plot frequency of alternate alleles for each of the 42 populations
ggplot(data=dat, aes(x=freqAlt)) + geom_histogram() + facet_wrap(~population)
dev.off()

####export gds as reduced vcf
vcf.fn<- "pooledData2.vcf"
seqGDS2VCF(gds.file, vcf.fn, info.var=NULL, fmt.var=NULL, use_Rsamtools=TRUE,
           verbose=TRUE)

###lfmm
library(LEA)
dat3<- adList$data/(adList$data+rdList$data)
colnames(dat3) <- paste("snp", seqGetData(gds.file, "variant.id"), sep="")
rownames(dat3) <- seqGetData(gds.file, "sample.id")
  
write.lfmm(dat3,"pooled3.lfmm")

##run pca
 pc<- pca("pooled3.lfmm",K=10,center = TRUE, scale = FALSE)
 tw = tracy.widom(pc)
 pc.dt <- as.data.table(pc$projections)
 
 setnames(pc.dt, names(pc.dt), gsub("V", "PC", names(pc.dt)))
 pc.dt[,sampleId:=seqGetData(gds.file, "sample.id")]
 
 # pc.dt <- merge(pc.dt, samps, by="sampleId")
 pc.dt
##save
 write.csv(pc.dt, file="pc_pooled.csv")

 
 library(ggplot2)
 library(data.table)
 
 pc.dt <- fread(file="pc_pooled.csv")
 
 pc.k <- kmeans(pc.dt[,c("PC1", "PC2", "PC3", "PC4", "PC5"), with=F], centers=3)
 pc.dt[,cluster:=pc.k$cluster]
 
pdf(file="pca.pdf")
 ggplot() +
   geom_point(data=pc.dt, aes(x=PC1, y=PC2, color=as.factor(cluster))) 
dev.off()
