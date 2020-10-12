#################LEA###############
library(LEA)

setwd("~/Downloads/GitHub/simCline/biosampleresults/")

vcf.file<- "~/Downloads/GitHub/simCline/biosampleresults/pooledData.vcf"
lfmm.file<- "~/Downloads/GitHub/simCline/biosampleresults/pooledData.lfmm"
vcf2lfmm(vcf.file, lfmm.file, force = TRUE)
pc<- pca("~/Downloads/GitHub/simCline/biosampleresults/pooledData.lfmm")

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
