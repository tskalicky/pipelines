#
# Visualize featureCounts gene biotypes output - need parsed featureCounts output where each column represents sample and contains number of reads assigned to each feature
# 

for(inputCounts in list.files(path = ".", pattern = "biotype.biotype_counts.txt$")){
  featureName<-strsplit(inputCounts, ".biotype")[[1]][1] # Get input file name
  featureCounts<-read.table(inputCounts, row.names=1, header=T, stringsAsFactors = F) # Read featureCounts extracted counts
  colnames(featureCounts)<-as.data.frame(strsplit(colnames(featureCounts), "_R1_"), stringsAsFactors=FALSE)[1,]  # Rename columns - split after first "_"
  
  # Plot Pie plot for each sample
  pdf(file=paste0(featureName, ".gene_biotypes", ".pdf"), width = 10, height = 6)
  par(mfrow=c(2,2))
  for(featureColumnNum in 1:(ncol(featureCounts))){
    featureColumnMain<-as.matrix(featureCounts[, featureColumnNum]) # Get one column from the input
    pctMain<-round(featureColumnMain/(sum(featureColumnMain)/100), 2) # Make labels with percentages

    # Merge features < 5% to others, otherwise to many lables 
    pct<-pctMain[pctMain>=5] # Get abundant features
    pct<-c(pct, sum(pctMain[pctMain<5])) # Other features - added together to one
    lbls<-rownames(featureCounts)[pctMain>=5] # Labels
    lbls<-c(lbls, "other")
    lbls<-paste(lbls, pct) 
    lbls <- paste(lbls, "%", sep="") # ad % to labels     

    featureColumn<-featureColumnMain[pctMain>=5] # Get subset of the main table for the abundant
    featureColumn<-c(featureColumn, sum(featureColumnMain[pctMain<5])) # Get the others
    
    pie(featureColumn, labels = lbls, col=rainbow(length(lbls)), main=paste("Pie Chart of Gene Biotypes", 
        colnames(featureCounts)[featureColumnNum], sum(featureColumnMain), sep="\n"))
  }
  dev.off()
}
