#!/usr/bin/env Rscript

#Script requires customisation throughout, including population structure tested, ploidy, and filenames.

library(pegas)
library(adegenet)
library(poppr)

## Attempt at a DAPC structure detection and AMOVA analysis
#Read in all the loci in the file
info <- VCFloci("SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_filtered.vcf")
SNP <- is.snp(info)
x <- length(which(SNP))
x <- read.vcf("SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_filtered.vcf", from = 1, to = x)
y <- loci2genind(x)

#Change ploidy to 1 for V. inaequalis data
ploidy(y) <- 1

#Number of individuals for maximum number of considered clusters
a <- nInd(y) - 23

#Find the true number of genetic clusters.
grp <- find.clusters(y, n.pca=200, max.n.clust=a)
dapc1 <- dapc(y, grp$grp, n.pca=200)
scatter(dapc1)


#AMOVA analysis
#Note: classification into pathogens, non-pathogens and INTERMEDIATES not more informative than just a binary one
#Classify according to if pathogenic or not
other(y)$cultivar <- c("Bramley", "Bramley", "Bramley", "Bramley", "Bramley", "Bramley", "Bramley", "Cox", "Cox", "Cox", "Cox", "Cox", "Cox", "Cox", "Cox", "Worcester", "Worcester", "Worcester", "Worcester", "Worcester", "Worcester", "Worcester", "Worcester")


strata_df <- data.frame(other(y))
strata(y) <- strata_df
sink(file = "amova.txt", append = T, type = c("output", "message"), split = F)
amova.pegas <- poppr.amova(y, ~pathogen, method = "pegas")
amova.pegas

