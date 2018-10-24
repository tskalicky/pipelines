#
# Create a pieplot of gene expression showing with more than 5% expression of the total
# The rest is summed to "others"
# Designed mainly for miRNA/sRNA analysis - normalizes the counts just by percentage of the total expression in the table
#

INPUT_DIR<-"/home/jan/Desktop/test/2013" # Chimira input gene.count table
OUTPUT_DIR<-INPUT_DIR
INPUT_FILE<-paste(INPUT_DIR, "gene.counts", sep="/")

library("rio")
gene_counts<-rio::import(INPUT_FILE, format = "tsv")

rownames(gene_counts)<-gene_counts[,1]
gene_counts[,1]<-NULL
colnames(gene_counts)<-gsub("_R1_001.gz", "", colnames(gene_counts))

dir.create(OUTPUT_DIR, recursive = T)

pdf(file=paste0(OUTPUT_DIR, "/", "expression_pie.pdf"), width = 10, height = 6)
par(mfrow=c(2,2))

gene_counts[is.na(gene_counts)]<-0

for(featureColumnNum in 1:(ncol(gene_counts))){
  featureColumnMain<-as.matrix(gene_counts[, featureColumnNum]) # Get one column from the input
  pctMain<-round(featureColumnMain/(sum(featureColumnMain)/100), 2) # Make labels with percentages
  
  # Merge features < 5% to others, otherwise to many lables 
  pct<-pctMain[pctMain>=5] # Get abundant features
  pct<-c(pct, sum(pctMain[pctMain<5])) # Other features - added together to one
  lbls<-rownames(gene_counts)[pctMain>=5] # Labels
  lbls<-c(lbls, "other")
  lbls<-paste(lbls, pct) 
  lbls <- paste(lbls, "%", sep="") # ad % to labels     
  
  featureColumn<-featureColumnMain[pctMain>=5] # Get subset of the main table for the abundant
  featureColumn<-c(featureColumn, sum(featureColumnMain[pctMain<5])) # Get the others
  
  pie(featureColumn, labels = lbls, col=rainbow(length(lbls)), main=paste("Pie Chart of Gene Expression", 
    gsub("\\..*", "", colnames(gene_counts)[featureColumnNum]), sum(featureColumnMain), sep="\n"))
}
dev.off()

# Write frequencies table
frequencies<-prop.table(as.matrix(gene_counts), margin=2)*100
write.table(frequencies, file=paste(OUTPUT_DIR, "mirna_freq.tsv", sep="/"), sep="\t", col.names=NA)
