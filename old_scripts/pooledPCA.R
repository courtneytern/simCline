#################LEA###############
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("LEA")

library(LEA)

# setwd("~/Downloads/GitHub/simCline/biosampleresults/")
setwd("/scratch/cat7ep/simCline/biosampleresults/")

# vcf.file<- "/scratch/cat7ep/simCline/biosampleresults/pooledData2.vcf"
lfmm.file<- "/scratch/cat7ep/simCline/biosampleresults/pooled.lfmm"
# vcf2lfmm(vcf.file, lfmm.file, force = TRUE)
#vcf.file<- "~/Downloads/GitHub/simCline/biosampleresults/pooledData2.vcf"
#lfmm.file<- "~/Downloads/GitHub/simCline/biosampleresults/pooledData.lfmm"

pdf(file="/scratch/cat7ep/simCline/biosampleresults/pca.pdf")

pc<- pca("/scratch/cat7ep/simCline/biosampleresults/pooled.lfmm",K=10)

# Display information on analysis.
show(pc)
# Summarize analysis.
summary(pc)
par(mfrow=c(2,2))
# Plot eigenvalues.
plot(pc, lwd=5, col="blue", cex = .7, xlab=("Factors"), ylab="Eigenvalues")
# PC1-PC2 plot.
plot(pc$projections)
# PC3-PC4 plot.
plot(pc$projections[,3:4])
# Plot standard deviations.
plot(pc$sdev)

dev.off()