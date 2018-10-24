#####################################################################
###R script to visualize simple statistics from HTSeq-count output###
###This is usefull for basic comparison of samples while analyzing###
###RNA-Seq data                                                   ###
#####################################################################
#Define extra functions
#####################################################################
###Because we will use ggplot2 library and we will want to output
###more graphs next to each other we have to create a function
###can do it. ggplot2 library doesn't work with "classic"
###par(mfrow=c(number,number2)
#
###Define multiplot function
###taken from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
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
#
#first specify input directory with all input files
input_dir<-"/home/jan/Desktop/test_STAR"
#
#loop for loading all input files
counts <- NULL; counts_info <- NULL; names <- NULL #"allocate" empty objects loop for reading all input files and merging them to one matrix
for (i in dir(input_dir))#speficify in which dir contains input files and store the name of the input file to variable i
  {
  counts_1 <- read.table(paste(input_dir, i, sep=""))#load input file to temporary object "counts_1"
  counts_info_1 <- counts_1[-(1:(nrow(counts_1)-5)), ]#store last five rows to new temporary object "counts_info_1"
                                                      #HTSeq-count output files have extra statistics on five last rows
  counts_info <- cbind(counts_info, counts_info_1[, 2])#merge all statistics from all input files to one object "counts_info"
  counts_1 <- counts_1[1:(nrow(counts_1)-5), ]#remove statistic rows from temporary object "counts_1" keeping only read counts for each gene in
                                              #input file
  counts <- cbind(counts, counts_1[,2])#merge temporary "counts_1" object to final object "counts" which will later contain all samples from input
                                       #directory one per column
  names_1<-sub("^([^.]*).*", "\\1", i)#store name of current input file "i" without all extensions after first "."
  names<-cbind(names,names_1)#merge all input file names to one final object "names"
}
#
counts <- data.frame(ensemble=as.character(counts_1[,1]), counts)#transform "counts" to data.frame "counts" and add names of genes stored in first
                                                                 #column of "counts_1" and name it as "ensemble"
str(counts)#visual check of current data.frame "counts"
colnames(counts)<-c("ensemble",names)#rename column names in "counts" object
#
counts_info <- data.frame(t(counts_info))#transpose "counts_info" object
colnames(counts_info) <- as.character(counts_info_1[,1])#rename column names in "count_info"
counts_info <- cbind(Samples=1:nrow(counts_info), counts_info)#add column "Samples" to "counts_info" depending on number of input samples 
counts_info[,1]<-names[1,]#replace numbers in column "Samples" in "counts_info" with actual names of input samples
#
#creating the plots using ggplot2 library
library(ggplot2)#load ggplot2 library
require(reshape)#make sure library reshape is loaded
#
melted_counts = melt(counts, id.vars= "ensemble")#melt "counts" to form suitable for casting
colnames(melted_counts)[2]<-"Samples"#rename second column in "melted_counts" 
#
#plot to compare number of reads for present genes and their density
#plots all samples on top of each other making it easy to visually compare them
p1<-ggplot(melted_counts, aes(x=value)) + geom_density(aes(color=Samples, fill=Samples),  alpha=0.01) + scale_x_log10() + theme_bw() + 
  ggtitle("Read number distribution per sample") + xlab("Number of reads") + ylab("Density")#create first plot and store it under "p1"
#
#plot to compare one samples versus another one
#p2<-ggplot(counts, aes(x=X1, y=X2)) + geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() + theme_bw()
#
#plot all statistics (last 5 rows from each input file) in separate columns
p3<-ggplot(melt(counts_info, id.vars=1), aes(x=variable, y=value, group=as.factor(Samples))) + 
  geom_bar(aes(fill=Samples, color=Samples), stat="identity", position="dodge") + 
  theme_bw() + ggtitle("Extra statistics per sample") + xlab("Statistics") + ylab ("Value")
#
#plot all graphs in one picture
multiplot(p1, p3, cols=2)#you can specify more than two plots (here "p1" and "p3") by just adding more stored plots
                         #you can modify "cols=2" to create more/less columns in final graph