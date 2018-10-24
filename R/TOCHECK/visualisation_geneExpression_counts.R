#####################################################################
# R script to visualize simple statistics from STAR gene count output - DESIGNED FOR ILLUMINA STRANDED (REVERSE STRAND) PROTOCOL
# FOR OTHER PROTOCOL CHANGE [, 4] IN GENE COUNTS TO APPROPRIATE COLUMN - [, 2] FOR UNSTRANDED, [,3] FOR STRANDED FORWARD STRAND PROTOCOL
# This is usefull for basic comparison of samples while analyzing RNA-Seq data
#####################################################################

# Define extra functions
# Because we will use ggplot2 library and we will want to output more graphs next to each other we have to create a function can do it. ggplot2 library doesn't work with "classic" par(mfrow=c(number,number2)


# Define multiplot function
# Taken from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, cols) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
}
#####################################################################
#Actual script
#####################################################################

# First specify input directory with all input files
input_dir<-"/home/jan/Desktop/test_STAR"

#loop for loading all input files
counts <- NULL; counts_info <- NULL; names <- NULL # "allocate" empty objects loop for reading all input files and merging them to one matrix
for (i in dir(input_dir)) # Speficify in which dir contains input files and store the name of the input file to variable i
{
  counts_1 <- read.table(paste(input_dir, i, sep="/")) # Load input file to temporary object "counts_1"
  counts_info_1 <- counts_1[(1:4), ] # Store last five rows to new temporary object "counts_info_1"
                                                      #HTSeq-count output files have extra statistics on five last rows, STAR on first 4
  counts_info <- cbind(counts_info, counts_info_1[, 4]) # Merge all statistics from all input files to one object "counts_info"
  counts_1 <- counts_1[5:(nrow(counts_1)), ] # Remove statistic rows from temporary object "counts_1" keeping only read counts for each gene in input file
  counts <- cbind(counts, counts_1[,4]) # Merge temporary "counts_1" object to final object "counts" which will later contain all samples from input directory one per column
  names_1<-sub("^([^.]*).*", "\\1", i)#store name of current input file "i" without all extensions after first "."
  names<-cbind(names,names_1) # Merge all input file names to one final object "names"
}

counts <- data.frame(ensemble=as.character(counts_1[,1]), counts) # Transform "counts" to data.frame "counts" and add names of genes stored in first column of "counts_1" and name it as "ensemble"

str(counts) # Visual check of current data.frame "counts"
colnames(counts)<-c("ensemble",names) # Rename column names in "counts" object

counts_info <- data.frame(t(counts_info)) # Transpose "counts_info" object
colnames(counts_info) <- as.character(counts_info_1[,1]) # Rename column names in "count_info"
counts_info <- cbind(Samples=1:nrow(counts_info), counts_info) # Add column "Samples" to "counts_info" depending on number of input samples 
counts_info[,1]<-names[1,] # Replace numbers in column "Samples" in "counts_info" with actual names of input samples

# Creating the plots using ggplot2 library
library(ggplot2) # Load ggplot2 library
require(reshape) # Make sure library reshape is loaded

melted_counts = melt(counts, id.vars= "ensemble") # Melt "counts" to form suitable for casting
colnames(melted_counts)[2]<-"Samples" # Rename second column in "melted_counts" 

# Plot to compare number of reads for present genes and their density
# Plots all samples on top of each other making it easy to visually compare them
p1<-ggplot(melted_counts, aes(x=value)) + geom_density(aes(color=Samples, fill=Samples),  alpha=0.01) + scale_x_log10() + theme_bw() + 
  ggtitle("Read number distribution per sample") + xlab("Number of reads") + ylab("Density") # Create first plot and store it under "p1"

# Plot to compare one samples versus another one
#p2<-ggplot(counts, aes(x=X1, y=X2)) + geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + theme_bw()

# Plot all statistics (last 5 rows from each input file) in separate columns
p3<-ggplot(melt(counts_info, id.vars=1), aes(x=variable, y=value, group=as.factor(Samples))) + 
  geom_bar(aes(fill=Samples, color=Samples), stat="identity", position="dodge") + 
  theme_bw() + ggtitle("Extra statistics per sample") + xlab("Statistics") + ylab ("Value")

# Plot all graphs in one picture
p1 # Gene count statistics
p3 # Other STAR statistics
multiplot(p1, p3, cols=2) # You can specify more than two plots (here "p1" and "p3") by just adding more stored plots
                          # You can modify "cols=2" to create more/less columns in final graph