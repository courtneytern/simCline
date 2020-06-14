#Convert VCF to GDF and filter 

library(SeqArray)

vcf.file <- "/project/berglandlab/courtney/simCline/pooledData.vcf"
gds.output <- "/project/berglandlab/courtney/simCline/pooledData.gds"
seqVCF2GDS(vcf.fn = vcf.file, out.fn = gds.output)