#
# Check expression profiles for RNA-Seq experiment
# Needs count table as the input from RSEM
# Be careful! This script DOES NOT remove rRNA and other contaminations - this would have to remove in input count tables
#

rm(list=ls(all=TRUE))

dir<-"/home/jan/Data/projects/blazek_rnaseq/results/gene_counts/rsem" # Set input dir
OUTPUT_DIR="/home/jan/Data/projects/blazek_rnaseq/results/exp_profile"
dir.create(OUTPUT_DIR, recursive = T)


SAMPLES_REMOVE<-c("")

library("tximport")

setwd(dir)
files <- list.files(path = dir, pattern = "rsem.genes.results$", include.dirs = T) # List input files
names(files) <- unlist(lapply(files, function(x) strsplit(x, "_R1_")[[1]][1])) # Assign sample names
names(files) <- unlist(lapply(files, function(x) strsplit(x, "_1_")[[1]][1])) # Assign sample names

names(files)
files<-files[!(names(files) %in% c(SAMPLES_REMOVE))]

txi.rsem <- tximport(files, type = "rsem") # Read gene counts
head(txi.rsem$counts) # See header

setwd(OUTPUT_DIR)

# Normalize by DESeq2
library("DESeq2")
library("RColorBrewer")

txi<-txi.rsem
cts <- txi$counts 
txi$length[txi$length == 0] <- 1 # If gene has length 0 replace it with 1 to avoid error later, might be source of bias and/or error; https://support.bioconductor.org/p/84304/
anyNA(cts)
normMat <- txi$length
normMat <- as.data.frame(normMat/exp(rowMeans(log(normMat))))
anyNA(normMat)

conds<-as.factor(seq(1:ncol(cts)))

coldata<-as.data.frame(t(t(conds)))
colnames(coldata)<-"condition"
rownames(coldata) <- colnames(cts)
coldata<-as.data.frame(coldata)

dds <- DESeqDataSetFromTximport(txi, coldata, ~condition)

#dds<-DESeqDataSetFromMatrix(countData = mrcounts, colData = coldata, design = ~condition)
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

