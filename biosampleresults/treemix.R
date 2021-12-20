## This script will execute TreeMix for pooled and individual samples 
## To be run on Rivanna

# module load gcc/7.1.0  openmpi/3.1.4  R
# R 


setwd("/scratch/cat7ep/simCline/biosampleresults/treemixOutputs")
# setwd("~/Downloads/GitHub/simCline/biosampleresults")

source("/scratch/cat7ep/simCline/biosampleresults/treemix-1.13/src/plotting_funcs.R")
# run the following for m0,1,2,3
pdf("./m3.pdf")
plot_tree("simcline_treeOut_m3")
dev.off()
