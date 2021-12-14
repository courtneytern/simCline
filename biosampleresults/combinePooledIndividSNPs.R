# Get only the SNPs that overlap in both pooled and individual populations
# Make a venn diagram showing this overlap 
## To be run on Rivanna

# module load gcc/7.1.0  openmpi/3.1.4  R
# R 

library(SeqArray)
library(ggplot2)
library(data.table)
library(dplyr)

setwd("/scratch/cat7ep/simCline/biosampleresults")
#setwd("~/Downloads/GitHub/simCline/biosampleresults")
individ.gds<- seqOpen("./individ.gds")
pooled.gds<- seqOpen("./pooled.gds")

# get the intersection of chromosome and position between pooled and individ samples
# individ.snps<- seqGetData(individ.gds,"$chrom_pos")
# pooled.snps<- seqGetData(pooled.gds,"$chrom_pos")
# intersectPos<- intersect(individ.snps,pooled.snps)
#chroms<- c("Dsim_Scf_2L", "Dsim_Scf_2R", "Dsim_Scf_3L","Dsim_Scf_3R","Dsim_Scf_X" )

# make table with the chromosome and position assoc with variant id
individ.dt <- data.table(chr_pos=seqGetData(individ.gds,"$chrom_pos"),
                     chr=seqGetData(individ.gds, "chromosome"),
                     pos=seqGetData(individ.gds, "position"))
pooled.dt <- data.table(chr_pos=seqGetData(pooled.gds,"$chrom_pos"),
                         chr=seqGetData(pooled.gds, "chromosome"),
                         pos=seqGetData(pooled.gds, "position"))
joined_snps<- intersect(individ.dt,pooled.dt)
### or inner_join()

#set filter to the snps that are common between both
#indvidIntersect<- seqSetFilterPos(individ.gds,chr=chroms, pos=intersectPos)
# individTreemix<- read.table("./individTreemixInput.txt")
# pooledTreemix<- read.table("./pooledTreemixInput.txt")

# filter by only the chr/pos in common between ind and pooled
# seqSetFilterPos(individ.gds,chr=joined_snps$chr,pos=joined_snps$pos)
# seqSetFilterPos(pooled.gds,chr=joined_snps$chr,pos=joined_snps$pos)

########### 
## rerun the treemix formatting for each 

## Pooled  #### 
pooledPath<- "/scratch/cat7ep/simCline/biosampleresults/pooled.gds"
makePooledTreemix<- function(pooledPath)  {
  gds.file<- seqOpen(pooledPath,allow.duplicate = T)
  # keep only joined snps
  seqSetFilterPos(gds.file,chr=joined_snps$chr,pos=joined_snps$pos)

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

  ## Get AD/RD list
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
  dat.ag2 <- dat[,list(nmissing=mean(is.na(ad)), aveAD=mean(ad, na.rm=T), aveRD=mean(rd, na.rm=T),
                      freqAlt=sum(ad, na.rm=T)/sum(ad+rd, na.rm=T),
                       chrom= chromosome[1], pos= position[1]),
                 list(variant.id)]
  # aggregate duplicate Machado rows
  agg<- dat %>% group_by(population,variant.id,position,chromosome) %>% summarise_all(sum)
  agg<- as.data.table(agg)
  # now calc freq alt
  agg[,freqAlt:=ad/(ad+rd)]

  # calc ad/rd per variant id -- merge dat.ag2 aveAD and aveRD with dat
  setkey(agg,variant.id)
  setkey(dat.ag2,variant.id)
  # left outer join. necessitate keeping all of dat
  ## subset dat.ag2 to be just the variant id, aveAD, aveRD, chrom, pos
  merged <- merge(agg, dat.ag2[,c(1,3,4,6,7)], all.x=TRUE)

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
  # long to wide
  datw <- dcast(merged, chrom+pos ~ population, value.var="newCol")
  # drop the var.id column
  #datw <- datw[,-1]
  # rename columns to something shorter
  colnames(datw)
  newColNames <- c("chrom","pos",
                   "Barghi:FL:Tallahassee","ES_Gim_14_34","ES_Gim_14_35","ES_Gim_16_33",
                   "ES_Pur_16_35","FR_Got_15_48","IT_Mez_15_43","IT_Mez_15_44","IT_Tre_16_15",
                   "Machado:Linvilla:65","Machado:Linvilla:50","PT_Rec_15_16","Sedghifar:SC:Conway",
                   "Sedghifar:ME:Fairfield","Sedghifar:FL:Miami","Sedghifar:Panama",
                   "Sedghifar:NJ:Princeton","Sedghifar:RI:Providence","Sedghifar:NC:Raleigh",
                   "Sedghifar:VA:Richmond","Sedghifar:GA:Savannah" )

  colnames(datw) <- newColNames
  datw
}
pooledOutput<- makePooledTreemix(pooledPath)
write.table(pooledOutput,"/scratch/cat7ep/simCline/biosampleresults/treemixPooledOut.txt",row.names = F)


## Individual  #### 
individPath<- "/scratch/cat7ep/simCline/biosampleresults/individ.gds"
makeIndividTreemix<- function(individPath) {
  individ.gds<- seqOpen(individPath,allow.duplicate = T)
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
  
  popTreemix<- function(popSamps){    
    seqResetFilter(individ.gds)
    seqSetFilter(individ.gds, sample.id=popSamps)
    #TEST: seqSetFilter(individ.gds,variant.id=samp.ids)
    
    dosage.mat<- seqGetData(individ.gds,"$dosage_alt")
    numAlt<- apply(dosage.mat,2, function(x) {sum(x==2,na.rm=T)*2 + sum(x==1,na.rm=T)} )
    numRef<- apply(dosage.mat,2, function(x) {sum(x==0,na.rm=T)*2 + sum(x==1,na.rm=T)} )
    paste(numAlt,numRef,sep=",")
  }#end popTreemix
  
  # init matrix with chromosome and position; as many rows as snps
  treemixTable<- matrix(c(seqGetData(individ.gds,"chromosome"),seqGetData(individ.gds,"position")),
                        nrow=length(seqGetData(individ.gds,"variant.id")))
  for(samps in popSampsList){
    treemixTable<- cbind(treemixTable,popTreemix(samps))
  }# end for
  colnames(treemixTable)<- c("chr","pos",unique(ind_metadata$population))
  treemixTable
}
individOutput<- makeIndividTreemix(individPath)
write.table(individOutput,"/scratch/cat7ep/simCline/biosampleresults/treemixIndividOut.txt",row.names = F)



