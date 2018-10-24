#
# Import counts from RSEM using tximport http://dx.doi.org/10.12688/f1000research.7563.1
# http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#rsem
#

#Note: there are two suggested ways of importing estimates for use with gene-level differential expression 
# methods. The first method, which we show below for edgeR and for DESeq2, is to use the estimated counts 
# from the quantification tools, and additionally to use the transcript-level abundance estimates to 
# calculate an offset that corrects for changes to the average transcript length across samples. The code 
# examples below accomplish these steps for you, keeping track of appropriate matrices and calculating these offsets.

library("tximport")
# library("dplyr")

dir<-"/home/jan/Data/projects/blazek_rnaseq/results/gene_counts/rsem" # Set input dir

removeSamples<-c("UHRR_S9", "UHRR_S10") # Can be empty

setwd(dir)

files <- list.files(path = dir, pattern = "rsem.genes.results$", include.dirs = T) # List input files
names(files) <- unlist(lapply(files, function(x) strsplit(x, "_R1_")[[1]][1])) # Assign sample names
names(files) <- unlist(lapply(files, function(x) strsplit(x, "_1_")[[1]][1])) # Assign sample names

files<-files[!(names(files) %in% c(removeSamples))]

txi.rsem <- tximport(files, type = "rsem") # Read gene counts
head(txi.rsem$counts) # See header

# For edgeR

txi<-txi.rsem
cts <- txi$counts 
txi$length[txi$length == 0] <- 1 # If gene has length 0 replace it with 1 to avoid error later, might be source of bias and/or error; https://support.bioconductor.org/p/84304/
anyNA(cts)
normMat <- txi$length
normMat <- as.data.frame(normMat/exp(rowMeans(log(normMat))))
anyNA(normMat)
#normMat[is.na(normMat)]<-0 # Might be source of error and/or bias!

# For DESeq2
library("DESeq2")
conds<-as.factor(seq(1:ncol(cts)))

coldata<-as.data.frame(t(t(conds)))
colnames(coldata)<-"condition"
rownames(coldata) <- colnames(cts)
coldata<-as.data.frame(coldata)

dds <- DESeqDataSetFromTximport(txi, coldata, ~condition)

library("edgeR")
o <- log(calcNormFactors(cts/normMat, na.rm = T)) + log(colSums(cts/normMat, na.rm = T)) # Must not contain NA values
y <- DGEList(cts)
y$offset <- t(t(log(normMat)) + o)
# y is now ready for estimate dispersion functions see edgeR User's Guide