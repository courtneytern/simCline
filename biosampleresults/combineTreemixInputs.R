# merge pooled and individ treemix input files
## To be run on rivanna

# module load gcc/7.1.0  openmpi/3.1.4  R
# R 

library(data.table)

setwd("/scratch/cat7ep/simCline/biosampleresults")
individInput<- fread("./treemixIndividInput_FINAL.txt",header = TRUE ,sep=" ")
pooledInput<- fread("./treemixPooledInput_FINAL.txt",header = TRUE ,sep=" ")

setkeyv(individInput,c("chr","pos"))
setkeyv(pooledInput,c("chrom","pos"))
combined<- merge(individInput,pooledInput,by.x=c("chr","pos"),by.y=c("chrom","pos"))

# drop the chr and pos cols now 
finalInput<- combined[,c(-1,-2)]
fwrite(finalInput,file="./treemixCOMBINEDinput.txt",quote=F,sep = " ")

# ready to run treemix.sh now