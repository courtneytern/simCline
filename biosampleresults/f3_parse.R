## Parse through f3 output
## To be run on Rivanna

# module load gcc/7.1.0  openmpi/3.1.4  R
# R 

############
## Setup ###
############
library(tidyr)

f3_output<- read.table("/scratch/cat7ep/simCline/biosampleresults/f3_outputs/FINAL_f3_output.txt",
                       sep=" ", fill=TRUE)
## f3_output<- read.table("~/Downloads/FINAL_f3_output.txt", sep=" ",fill=TRUE)

# sort populations (both pooled and individ) by continent
amr_pops<- c("USA.Linvilla.2011.Aug","USA.Linvilla.2011.Sep","USA.Charlottesville.2010.Sep","USA.Bowdoin.2011.Oct",
             "USA.Homestead.2011.May","USA.Linvilla.2011.Nov","USA.Zuma.2012.Feb","Barghi:FL:Tallahassee",
             "Machado:Linvilla:65","Machado:Linvilla:50","Sedghifar:ME:Fairfield","Sedghifar:FL:Miami","Sedghifar:Panama",
             "Sedghifar:NJ:Princeton","Sedghifar:RI:Providence","Sedghifar:NC:Raleigh","Sedghifar:VA:Richmond",
             "Sedghifar:GA:Savannah","Sedghifar:SC:Conway")
eur_pops<- c("Israel.NFS.2014.Oct","Israel.SFS.2014.Oct","ES_Gim_14_34","ES_Gim_14_35","ES_Gim_16_33",
             "ES_Pur_16_35","FR_Got_15_48","IT_Mez_15_43","IT_Mez_15_44","IT_Tre_16_15","PT_Rec_15_16")
afr_pops<- c("Kenya.Nairobi.2006.NA","Madagascar.Joffreville.2002.NA")

############
## Parse ###
############
library(ggplot2)

colnames(f3_output)<- c("tree","f3_stat","std_err","z_score")
# separate the tree elements into A, B, and C cols 
f3_sep<- separate(f3_output,1,sep="(;|,)",into=c("A","B","C"))

# keep only AMR;AFR,EUR trees
## AMR;EUR,AFR == AMR;AFR,EUR
keep_trees<- f3_sep[(is.element(f3_sep$A,amr_pops))&(is.element(f3_sep$B,afr_pops))&
               (is.element(f3_sep$C,eur_pops)),]
keep_trees<- keep_trees[order(keep_trees$z_score),]

ggplot(data=keep_trees) + geom_tile(aes(x=A,y=C,fill=z_score)) + facet_grid(cols=vars(B)) +
  theme(axis.text.x=element_text(angle=90)) + xlab("North American Populations") + ylab("European Populations")

# split by African population
madagascar<- keep_trees[keep_trees$B=="Madagascar.Joffreville.2002.NA",]
kenya<- keep_trees[keep_trees$B=="Kenya.Nairobi.2006.NA",]

madagascar_z<- madagascar[,c("A","C","z_score")]
m<- ggplot(data=madagascar_z) + geom_tile(aes(x=A, y=C, fill=z_score)) + 
  theme(axis.text.x=element_text(angle=90)) + 
  labs(title="Z-score Heatmap",subtitle="Ancestral Population: Madagascar") +
  xlab("North American Populations") + ylab("European Populations")

kenya_z<- kenya[,c("A","C","z_score")]
k<- ggplot(data=kenya_z) + geom_tile(aes(x=A, y=C, fill=z_score)) + 
  theme(axis.text.x=element_text(angle=90)) + 
  labs(title="Z-score Heatmap",subtitle="Ancestral Population: Kenya") +
  xlab("North American Populations") + ylab("European Populations")


###testing
# f3_sep[(is.element(f3_sep$A,amr_pops)),]
# f3_sep[(is.element(f3_sep$A,amr_pops))&(is.element(f3_sep$B,afr_pops)),]
