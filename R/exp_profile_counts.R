#
# Check expression profiles for RNA-Seq experiment
# Needs count table as the input from featureCounts, htseq-counts, STAR counts
# Be careful! This script DOES NOT remove rRNA and other contaminations - this would have to remove in input count tables
#

rm(list=ls(all=TRUE))

INPUT_COUNTS<-"/home/jan/Projects/pernisova_rnaseq/results/counts/featureCounts/rev.featureCounts" # Gene count table
OUTPUT_DIR="/home/jan/Projects/pernisova_rnaseq/results/qc/exp_profile"
dir.create(OUTPUT_DIR, recursive = T)
setwd(OUTPUT_DIR)

removeSamples<-c("")

mrcounts<-read.table(INPUT_COUNTS, header=TRUE, row.names=1)
mrcounts<-mrcounts[!rownames(mrcounts) %in% c("N_ambiguous", "N_multimapping", "N_noFeature", "N_unmapped"), ] # In case we have STAR counts
mrcounts<-mrcounts[,!colnames(mrcounts) %in% c("Chr", "Start", "End", "Strand", "Length")]# In case we have there featureCounts counts

colnames(mrcounts)<-gsub("^X", "", colnames(mrcounts)) # Rename of columns if necessary
colnames(mrcounts)<-gsub("\\..*", "", colnames(mrcounts)) # Rename of columns if necessary
colnames(mrcounts)<-gsub("_R1.*", "", colnames(mrcounts)) # Rename of columns if necessary

# Any samples to remove?
colnames(mrcounts)
SAMPLES_REMOVE<-c("") # can be empty
mrcounts<-mrcounts[,!(colnames(mrcounts) %in% SAMPLES_REMOVE)] # Remove the samples

# Normalize by DESeq2
library("DESeq2")
library("RColorBrewer")
coldata<-as.data.frame(as.factor(seq(1:ncol(mrcounts))))
rownames(coldata)<-colnames(mrcounts)
colnames(coldata)<-"condition"
dds<-DESeqDataSetFromMatrix(countData = mrcounts, colData = coldata, design = ~condition)
dds$condition<-factor(dds$condition, levels=levels(coldata$condition))# Make sure the factors are correct
dds<-estimateSizeFactors(dds)

counts_main_both<-counts(dds, normalized=TRUE) # Save normalized counts

color<-rainbow(n = ncol(counts_main_both))

logcounts <- log(counts_main_both[,1], 10) 
d <- density(logcounts)

# Get min and max for the plot
minForPlotX<-min(d$x)
maxForPlotX<-min(d$x)
maxForPlotY<-min(d$y)

# TODO Replace with recode() from dplyr
for (s in 2:ncol(counts_main_both)){
  logcounts <- log(counts_main_both[,s],10) 
  d <- density(logcounts)
  if(min(d$x)<minForPlotX){
    minForPlotX<-min(d$x)
  }
  if(max(d$x)>maxForPlotX){
    maxForPlotX<-max(d$x)
  }
  if(max(d$y)>maxForPlotY){
    maxForPlotY<-max(d$y)
  }
}

# Plot expression profiles
pdf(paste0("normalized_gene_expression_check.pdf"), width=10)
  plot(d, main="Expression profiles", xlim=c(minForPlotX*1.3, maxForPlotX*1.3), ylim=c(0, maxForPlotY*1.3), xlab=paste0("Log10 of normalized 
         expression per gene (DESeq2)"), ylab="Density", col=color[1])
  
  for (s in 2:ncol(counts_main_both)){
    logcounts <- log(counts_main_both[,s],10) 
    d <- density(logcounts)
    lines(d, col=color[s])
  }
  
  legend("topright", legend=colnames(counts_main_both), col = color, cex = .5, fill=color)
dev.off()

### Plot correlations between the samples

library("corrplot")
Mp<-cor(counts_main_both, method="pearson") # Spearman = rank, pearson = value
Ms<-cor(counts_main_both, method="spearman") # Spearman = rank, pearson = value
pdf("samples_correlation.pdf", width = 12, height = 10)
  corrplot(Ms, method="number", main="Spearman correlation of DESeq2 norm counts")
  corrplot(Mp, method="number", main="Pearson correlation of DESeq2 norm counts")
dev.off()

pdf("samples_correlation_mixedHclust.pdf", width = 12, height = 10)
  corrplot.mixed(Ms, order ="hclust", main="Spearman correlation of DESeq2 norm counts")
  corrplot.mixed(Mp, order ="hclust", main="Pearson correlation of DESeq2 norm counts")
dev.off()

### Heatmaps
rawcounts<-counts(dds, normalized=FALSE) # Save raw counts
normcounts<-counts(dds, normalized=TRUE) # Save normalized counts
log2counts<-log2(normcounts+1) # Save log2 of normalized counts

vsd<-varianceStabilizingTransformation(dds) # Save counts tranformed with variance Stabilizing Transformation
vstcounts<-assay(vsd)

# Heatmaps
cond_colours<-rainbow(n = ncol(rawcounts))[as.factor(coldata$condition)]
names(cond_colours)<-coldata$condition
hmcol<-colorRampPalette(brewer.pal(9, "GnBu"))(100)

library("gplots")

pdf(file="heatmaps_samples.pdf")
  heatmap.2(cor(vstcounts), trace="none", col=hmcol, main="Sample to Sample Correlation (VST)", margins=c(9,7)) # RowSideColors=cond_colours, 
  heatmap.2(cor(log2counts), trace="none", col=hmcol, main="Sample to Sample Correlation (Log2)", margins=c(9,7)) # RowSideColors=cond_colours, 
  heatmap.2(cor(rawcounts), trace="none", col=hmcol, main="Sample to Sample Correlation (Raw Counts)", margins=c(9,7)) # RowSideColors=cond_colours, 
dev.off()

# Get expression boxplots
library("pcaExplorer")
pdf("expression_boxplot.pdf")
  distro_expr(rlogTransformation(dds), plot_type = "boxplot")
dev.off()

# png("scatterplot_norm.png", width = 1000, height = 1000)
#   pairs(counts_main_both, main="Normalized counts")
# dev.off()

png("scatterplot_log10norm.png", width = 1000, height = 1000)
  pairs(log(counts_main_both+1, 10), main="log10 Normalized counts+1")
dev.off()

graphics.off()

