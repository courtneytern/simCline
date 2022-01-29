## This script will execute TreeMix for pooled and individual samples 
## To be run on Rivanna

# module load gcc/7.1.0  openmpi/3.1.4  R
# R 
### Run locally and save as 10x20 pdf


setwd("/scratch/cat7ep/simCline/biosampleresults/treemixOutputs")
# setwd("~/Downloads/GitHub/simCline/biosampleresults/treemixOutputs")

source("/scratch/cat7ep/simCline/biosampleresults/treemix-1.13/src/plotting_funcs.R")
# source("~/Downloads/GitHub/simCline/biosampleresults/treemix-1.13/src/plotting_funcs.R")
# run the following for m0,1,2,3
pdf("./treemixPlots.pdf")
par(adj=0,cex=0.5,mar=c(2,4,4,4),mfrow=c(2,2))
plot_tree("simcline_treeOut_m0") 
  mtext("(A)",side=3,adj=0,cex=3)
plot_tree("simcline_treeOut_m1")
  mtext("(B)",side=3,adj=0,cex=3)
plot_tree("simcline_treeOut_m2")
  mtext("(C)",side=3,adj=0,cex=3)
plot_tree("simcline_treeOut_m3")
  mtext("(D)",side=3,adj=0,cex=3)
dev.off()