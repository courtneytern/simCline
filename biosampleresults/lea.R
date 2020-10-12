#################LEA###############
library(LEA)

# setwd("~/Downloads/GitHub/simCline/biosampleresults/")
setwd("/scratch/cat7ep/simCline/biosampleresults/")

# vcf.file<- "pooledData.vcf"
# lfmm.file<- "pooledData.lfmm"
# vcf2lfmm(vcf.file, lfmm.file, force = TRUE)

pdf(file="/scratch/cat7ep/simCline/biosampleresults/pca.pdf")

pc<- pca("/scratch/cat7ep/simCline/biosampleresults/pooledData.geno")

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