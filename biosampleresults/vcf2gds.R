#Convert VCF to GDF and filter 

library(SeqArray)

vcf.file <- "~/Downloads/GitHub/simCline/biosampleresults/pooledData.vcf"
gds.output <- "~/Downloads/GitHub/simCline/biosampleresults/pooledData.gds"
seqVCF2GDS(vcf.fn = vcf.file, out.fn = gds.output)

gds.file <- seqOpen(gds.output)
snp.dt <- data.frame(chr=seqGetData(gds.file, "chromosome"),
                     pos=seqGetData(gds.file, "position"),
                     nAlleles=seqGetData(gds.file, "$num_allele"),
                     id=seqGetData(gds.file, "variant.id"),
                     seqMissing(gds.file, per.variant=T))
table(snp.dt$chr)
prop.table(table(snp.dt$nAlleles))

#now filter
nAlleles=seqGetData(gds.file, "$num_allele")
seqSetFilter(gds.file, variant.id=snp.dt[which(snp.dt$nAlleles==2), ]$id)
adList<- seqGetData(gds.file, "annotation/format/AD")
rdList<- seqGetData(gds.file, "annotation/format/RD" )
