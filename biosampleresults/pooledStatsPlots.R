# this generates the freqAlt plots per pooled population

library(SeqArray)
library(ggplot2)
library(data.table)
library(dplyr)
library(plyr) # for groupby and summarise_all 
library(R.utils) # for gzip

setwd("~/Downloads/GitHub/simCline/biosampleresults/")

##########################
### Setup/filtering ######
##########################

gds.output <- "./pooled.gds"
gds.file <- seqOpen(gds.output, allow.duplicate=T)
individ.gds <- seqOpen("./individ.gds")
snp.dt <- data.table(chr=seqGetData(gds.file, "chromosome"),
                     pos=seqGetData(gds.file, "position"),
                     nAlleles=seqGetData(gds.file, "$num_allele"),
                     id=seqGetData(gds.file, "variant.id"),
                     seqMissing(gds.file, per.variant=T))

# filter where not dmels
snp.dt <- snp.dt[grepl("Dsim_Scf_2L|Dsim_Scf_2R|Dsim_Scf_3L|Dsim_Scf_3R", chr)]
# keep only where 2 alleles
nAlleles=seqGetData(gds.file, "$num_allele")
# get ids for previous filters
ids <- snp.dt[nAlleles==2]$id
seqSetFilter(gds.file, variant.id=ids)
####sample 100,000 ids
# samp.ids <- as.numeric(sample(x=as.character(snp.dt[nAlleles==2]$id), size=10000))
# seqSetFilter(gds.file, variant.id=samp.ids)

####
adList<- seqGetData(gds.file, "annotation/format/AD")
rdList<- seqGetData(gds.file, "annotation/format/RD")
# compile variables of interest
dat <- data.table(population=rep(seqGetData(gds.file, "sample.id"), dim(adList$data)[2]),
                  variant.id=rep(seqGetData(gds.file, "variant.id"), each=dim(adList$data)[1]),
                  ad=expand.grid(adList$data)$Var1,
                  rd=expand.grid(rdList$data)$Var1,
                  position=rep(seqGetData(gds.file, "position"), each=dim(adList$data)[1]),
                  chromosome=rep(seqGetData(gds.file, "chromosome"), each=dim(adList$data)[1])
                  )
write.csv(dat, file="dat_pooledpops.csv")
## dat <- read.csv("dat_pooledpops.csv")
dat.ag2 <- dat[,list(nmissing=mean(is.na(ad)), aveAD=mean(ad, na.rm=T), aveRD=mean(rd, na.rm=T), 
                     freqAlt=sum(ad, na.rm=T)/sum(ad+rd, na.rm=T),
                     chrom= chromosome[1], pos= position[1]), 
               list(variant.id)]

##########################
### Make data tables #####
##########################

# aggregate duplicate Machado rows 
agg<- dat %>% group_by(population,variant.id,position,chromosome) %>% summarise_all(sum)
# save agg in new csv
write.csv(agg, file="aggregated_ad_rd.csv")
## agg <- read.csv("aggregated_ad_rd.csv")

agg<- as.data.table(agg)
# now calc freq alt 
agg[,freqAlt:=ad/(ad+rd)]

# calc ad/rd per variant id -- merge dat.ag2 aveAD and aveRD with dat 
setkey(agg,variant.id)
setkey(dat.ag2,variant.id)
# left outer join. necessitate keeping all of dat
## subset dat.ag2 to be just the variant id, aveAD, and aveRD
merged <- merge(agg, dat.ag2[,c(1,3,4)], all.x=TRUE)

# if NA in either ad or rd, construct treemix "x,y" with both avgs
## if not NA, use the actual val for "x,y"
# take as parameters: ad col, rd col, aveAD col, aveRD col 
checkNA <- function(a,r,aa,ar){
  if( is.na(a)|is.na(r) ) {
    return(paste(aa,ar,sep=","))
  } # if
  else
    return(paste(a,r,sep=","))
}
merged[,newCol:=mapply(checkNA,merged$ad,merged$rd,merged$aveAD,merged$aveRD)]
write.csv(merged,"merged.csv")

# #replace NA with 0
# ### faulty-- fix this 
# agg$ad[is.na(agg$ad)] <- 0
# agg$rd[is.na(agg$rd)] <- 0
# # paste together for treemix format
# agg[,newCol:=paste(ad,rd,sep=",")]
## for individ, sum ad and rd first per pop, then paste and datw

# try this with ad and rd first, then do newCol 
datw <- dcast(merged, variant.id ~ population, value.var="newCol")
# drop the var.id column
datw <- datw[,-1]
# rename columns to something shorter
colnames(datw)
newColNames <- c("Barghi:FL:Tallahassee","ES_Gim_14_34","ES_Gim_14_35","ES_Gim_16_33",
                 "ES_Pur_16_35","FR_Got_15_48","IT_Mez_15_43","IT_Mez_15_44","IT_Tre_16_15",
                 "Machado:Linvilla:65","Machado:Linvilla:50","PT_Rec_15_16","Sedghifar:SC:Conway",
                 "Sedghifar:ME:Fairfield","Sedghifar:FL:Miami","Sedghifar:Panama",
                 "Sedghifar:NJ:Princeton","Sedghifar:RI:Providence","Sedghifar:NC:Raleigh",
                 "Sedghifar:VA:Richmond","Sedghifar:GA:Savannah"
                 )

colnames(datw) <- newColNames
# Write to txt file to be zipped for treemix 
write.table(datw,file="treemix_input.txt",quote=FALSE,row.names=FALSE)
gzip(filename="treemix_input.txt",destname="treemix_input.gz")

#### vvv This one works if you need to aggregate with duplicates
# datw <- dcast(dat, variant.id ~ population, value.var=c("ad","rd"),fun.aggregate=sum)
# datw ## output for TreeMix 


#get unique chromosome and snps per chromosome
chroms<- unique(dat$chromosome)
length(chroms)

numSnpsList<- list()
for(i in 0:length(chroms) ){
  numSnpsList[i]<- length(which(dat$chromosome==chroms[i]))
}
length(unique(dat$variant.id))
#snps per chromosome as a table
numSnps<- data.table(Chromosome=chroms,numSnps=numSnpsList)
numSnps

#missing rate, mean RD, median RD, lower5th/upper 95th quantile, per population
dat.ag <- dat[,list(propMissing=mean(is.na(ad)), aveRD=mean(rd, na.rm=T), medRD=as.double(median(rd,na.rm=T)), 
                    lower5RD=quantile(rd, 0.05, na.rm=T), upper95RD=quantile(rd, 0.95, na.rm=T) ), 
              list(population)]
write.csv(dat.ag, file="pooled_pops_summary.csv")

########################
### Plots/tables #######
########################

# prop missing vs average read depth per population
p<- ggplot(data=dat.ag, aes(x=propMissing,y=aveRD)) + geom_point()
p

# nmissing histogram per snp 
q<- ggplot(data=dat.ag2, aes(x=nmissing)) + geom_histogram() 
q

# alt frequency histogram per snp
v<- ggplot(data=dat.ag2, aes(x=freqAlt)) + geom_histogram()
v

### Make summary table 
# check some quantiles for nmissing
quantile(dat.ag2$nmissing, c(.75,.8,.9,.95))

threshold<- quantile(dat.ag2$nmissing, .9) # set threshold
# sets list of variant id, then whether it passed the missing threshold or not 
nmissingPassed<- data.frame(chrom=dat.ag2$chrom,
                         pos=dat.ag2$pos,
                         var.id=dat.ag2$variant.id,
                         nmissing= dat.ag2$nmissing, 
                         thresholdPass= dat.ag2$nmissing<=threshold)
#summary of nmissing by variant id, remove pos and chrom info 
nmissingPassed2<- distinct(data.frame(var.id=dat.ag2$variant.id,
                            nmissing= dat.ag2$nmissing, 
                            thresholdPass= dat.ag2$nmissing<=threshold))

####### Make PDF of summary graphs #########
pdf(file="pooled_data_summary_plots.pdf")
# prop missing vs average read depth per population
p<- ggplot(data=dat.ag, aes(x=propMissing,y=aveRD)) + 
  labs(title='Prop Missing vs Average Read Depth Per Population') + geom_point()
p

x<- ggplot(data=dat.ag, aes(x=aveRD)) + 
  labs(title='Average Read Depth Per Population') + geom_histogram()
x

# nmissing histogram per snp 
q<- ggplot(data=dat.ag2, aes(x=nmissing)) + 
  labs(title='N Missing Per SNP') + geom_histogram() 
q

# alt frequency histogram per snp
v<- ggplot(data=dat.ag2, aes(x=freqAlt)) + 
  labs(title='Alternate Frequency Per SNP') + geom_histogram()
v

#alt alleles freq per population
w<- ggplot(data=dat, aes(x=freqAlt)) + 
  labs(title='Alternate Frequency Per Population') + geom_histogram() + facet_wrap(~population)
w
dev.off()

####export gds as reduced vcf
vcf.fn<- "pooledData_reduced.vcf"
seqGDS2VCF(gds.file, vcf.fn, info.var=NULL, fmt.var=NULL, use_Rsamtools=TRUE,
           verbose=TRUE)


#####################
##### LEA ###########
#####################

library(LEA)

###lfmm
dat3<- adList$data/(adList$data+rdList$data)
colnames(dat3) <- paste("snp", seqGetData(gds.file, "variant.id"), sep="")
rownames(dat3) <- seqGetData(gds.file, "sample.id")
  
write.lfmm(dat3,"pooled.lfmm")

### run pca
 pc<- pca("pooled.lfmm",K=10,center = TRUE, scale = FALSE)
 tw <- tracy.widom(pc)
 pc.dt <- as.data.table(pc$projections)
 
 setnames(pc.dt, names(pc.dt), gsub("V", "PC", names(pc.dt)))
 pc.dt[,sampleId:=seqGetData(gds.file, "sample.id")]
 # pc.dt <- merge(pc.dt, samps, by="sampleId")
 
##save datatable of pc components
 pc.k <- kmeans(pc.dt[,c("PC1", "PC2", "PC3", "PC4", "PC5"), with=F], centers=2)
 pc.dt[,cluster:=pc.k$cluster]
 write.csv(pc.dt, file="pc_pooled.csv")
 pc.dt <- fread(file="pc_pooled.csv")
 
pdf(file="pca_pooled.pdf")
 ggplot() +
   geom_point(data=pc.dt, aes(x=PC1, y=PC2, color=as.factor(cluster))) 
dev.off()
