#!/usr/bin/env /afs/ics.muni.cz/software/R/3.4.0/gcc/bin/Rscript
#
# Originally the line above was !/usr/bin/env Rscript; we modified it to fit MetaCentrum; https://stackoverflow.com/questions/7935708/rscript-issue-using-different-version-of-r
#
# Copied from https://github.com/ewels/NGI-RNAseq/blob/master/bin/dupRadar.r
#

# Command line argument processing
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("Usage: dupRadar.r <input.bam> <annotation.gtf> <paired/single> <strandness> <R-package-location (optional)>", call.=FALSE)
}
input_bam <- args[1]
annotation_gtf <- args[2]
paired_end <- if(args[3]=='paired') TRUE else FALSE
stranded <- as.numeric(args[4]) # '0' (unstranded), '1' (stranded) and '2' (reversely stranded)
input_bam_basename <- strsplit(input_bam, "\\.")[[1]][1]

# Load / install packages
if (length(args) > 4) { .libPaths( c( args[5], .libPaths() ) ) }
if (!require("dupRadar")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("dupRadar", suppressUpdates=TRUE)
  library("dupRadar")
}
if (!require("parallel")) {
  install.packages("parallel", dependencies=TRUE, repos='http://cloud.r-project.org/')
  library("parallel")
}

# Duplicate stats
#stranded <- 2
threads <- detectCores()
dm <- analyzeDuprates(input_bam, annotation_gtf, stranded, paired_end, threads)
write.table(dm, file=paste(input_bam_basename, "_dupMatrix.txt", sep=""), quote=F, row.name=F, sep="\t")

# 2D density scatter plot
pdf(paste0(input_bam_basename, "_duprateExpDens.pdf"))
duprateExpDensPlot(DupMat=dm)
title("Density scatter plot")
mtext(input_bam_basename, side=3)
dev.off()
fit <- duprateExpFit(DupMat=dm)
cat(
  paste("- dupRadar Int (duprate at low read counts):", fit$intercept),
  paste("- dupRadar Sl (progression of the duplication rate):", fit$slope),
  fill=TRUE, labels=input_bam_basename,
  file=paste0(input_bam_basename, "_intercept_slope.txt"), append=FALSE
)

# Get numbers from dupRadar GLM
curve_x <- sort(log10(dm$RPK))
curve_y = 100*predict(fit$glm, data.frame(x=curve_x), type="response")
# Remove all of the infinite values
infs = which(curve_x %in% c(-Inf,Inf))
curve_x = curve_x[-infs]
curve_y = curve_y[-infs]
# Reduce number of data points
curve_x <- curve_x[seq(1, length(curve_x), 10)]
curve_y <- curve_y[seq(1, length(curve_y), 10)]
# Convert x values back to real counts
curve_x = 10^curve_x
# Write to file
write.table(
  cbind(curve_x, curve_y),
  file=paste0(input_bam_basename, "_duprateExpDensCurve.txt"),
  quote=FALSE, row.names=FALSE
)

# Get information about genes covered by multimapping reads
dm$mhRate <- (dm$allCountsMulti - dm$allCounts) / dm$allCountsMulti # calculate the fraction of multimappers per gene

sink(paste0(input_bam_basename, "_multimapPerGene.txt"))
print("Genes covered exclusively by multimapped reads")
sum(dm$mhRate == 1, na.rm=TRUE) # how many genes are exclusively covered by multimappers
print("Genes covered exclusively by multimapped reads with RPKM > 5")
sum(dm$mhRate==1 & dm$RPKMMulti > 5, na.rm=TRUE)
dm[dm$mhRate==1 & dm$RPKMMulti > 5, "ID"]
sink()

# Multimapped reads per gene frequency
pdf(paste0(input_bam_basename, "_multimapPerGene.pdf"))
hist(dm$mhRate, 
     breaks=50, 
     col="red", 
     main="Frequency of multimapping rates", 
     xlab="Multimapping rate per gene", 
     ylab="Frequency"
#     ylim=c(0, max(dm$mhRate, na.rm=T)*1.3)
#     ylim=c(0, 100)
)
#title("Frequency of multimapping rates")
mtext(input_bam_basename, side=3)
dev.off()

# Distribution of expression box plot
pdf(paste0(input_bam_basename, "_duprateExpBoxplot.pdf"))
duprateExpBoxplot(DupMat=dm)
title("Percent Duplication by Expression")
mtext(input_bam_basename, side=3)
dev.off()

# Distribution of RPK values per gene
pdf(paste0(input_bam_basename, "_expressionHist.pdf"))
expressionHist(DupMat=dm)
title("Distribution of RPK values per gene")
mtext(input_bam_basename, side=3)
dev.off()

# Distribution of reads per expression bins
pdf(paste0(input_bam_basename, "_readDist.pdf"))
readcountExpBoxplot(DupMat=dm)
#title("Distribution of reads per expression bin")
mtext(input_bam_basename, side=3)
dev.off()

# Print sessioninfo to standard out
print(input_bam_basename)
citation("dupRadar")
sessionInfo()
