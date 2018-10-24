#
# Visualize preseq results
#
# Needs preseq_lc and preseq_c results
#
# Plots real data in solid color and estimates in dashed line
# https://davetang.org/muse/2013/07/10/how-deep-should-we-sequence/

INPUT_DIR<-"."
setwd(INPUT_DIR)

yield_files<-list.files(pattern = "\\.yield_estimates.txt$") # yield estimates = preseq_lc

merged_preseq<-NULL # Merged dataset

### Short results - up to 100M reads
for(i in yield_files){
  sample<-strsplit(i,"_trim")[[1]][1] # Get sample name; might be better to use "_R1_"
  lc<-read.table(file=i, header=T) # Read preseq_lc
  c<-read.table(file=list.files(pattern=glob2rx(paste0(sample, "*.estimates.txt")) ), header=T) # Read preseq_c by sample name - real measuerd data (estimates.txt)
  colnames(lc)[colnames(lc)=="TOTAL_READS"]<-"total_reads" # Rename column
#  lc<-lc[1:700,] # Get only first 700 rows - overall output from preseq is way too big
  lc<-lc[1:101,] # Get only first 101 rows (100M) reads - overall output from preseq is way too big
##  c<-c[1:700,] # Get only first 700 rows - overall output from preseq is way too big
  merged_lcc<-merge(lc, c, all.x=T, by="total_reads") # Merge preseq_lc and preseq_c
  merged_lcc$sample<-sample # Add sample name
  merged_preseq<-rbind(merged_preseq, merged_lcc) # Merge to the final table
}

# Plot using ggplot
library("ggplot2")
g1<-ggplot(merged_preseq, aes(x=total_reads, y=EXPECTED_DISTINCT, color=sample)) + # total reads vs expected distinct reads
  geom_line(linetype="dotted", size=0.75) + # Add estimates (preseq_lc) (dotted line)
  geom_line(aes(x=total_reads, y=distinct_reads, color=sample), size=1) # Add real measured data (preseq_c)(solid line)

g2<-g1 + geom_ribbon(aes(x=total_reads, ymin=LOWER_0.95CI, ymax=UPPER_0.95CI, color=sample, fill=sample), alpha=0.3) # Add intervals of confidence

# Add titles
g1 <- g1 + ggtitle("Preseq library complexity estimate - close look") +
  labs(x="Total number of reads", y="Expected distinct reads") +
  theme(plot.title = element_text(hjust = 0.5))

g2 <- g2 + ggtitle("Preseq library complexity estimate - close look with CI") +
  labs(x="Total number of reads", y="Expected distinct reads") +
  theme(plot.title = element_text(hjust = 0.5))


pdf("preseq_estimates_short.pdf")
  print(g1)
  print(g2)
dev.off()

### Longe results - up to 700M reads
for(i in yield_files){
  sample<-strsplit(i,"_trim")[[1]][1] # Get sample name; might be better to use "_R1_"
  lc<-read.table(file=i, header=T) # Read preseq_lc
  c<-read.table(file=list.files(pattern=glob2rx(paste0(sample, "*.estimates.txt")) ), header=T) # Read preseq_c by sample name - real measuerd data (estimates.txt)
  colnames(lc)[colnames(lc)=="TOTAL_READS"]<-"total_reads" # Rename column
  lc<-lc[1:701,] # Get only first 501 rows (500M) reads - overall output from preseq is way too big
##  c<-c[1:700,] # Get only first 700 rows - overall output from preseq is way too big
  merged_lcc<-merge(lc, c, all.x=T, by="total_reads") # Merge preseq_lc and preseq_c
  merged_lcc$sample<-sample # Add sample name
  merged_preseq<-rbind(merged_preseq, merged_lcc) # Merge to the final table
}

# Plot using ggplot
library("ggplot2")
g1 <- ggplot(merged_preseq, aes(x=total_reads, y=EXPECTED_DISTINCT, color=sample)) + # total reads vs expected distinct reads
  geom_line(linetype="dotted", size=0.75) + # Add estimates (preseq_lc) (dotted line)
  geom_line(aes(x=total_reads, y=distinct_reads, color=sample), size=1) # Add real measured data (preseq_c)(solid line)

g2 <- g1 + geom_ribbon(aes(x=total_reads, ymin=LOWER_0.95CI, ymax=UPPER_0.95CI, color=sample, fill=sample), alpha=0.3) # Add intervals of confidence

# Add titles
g1 <- g1 + ggtitle("Preseq library complexity estimate - overall") +
  labs(x="Total number of reads", y="Expected distinct reads") +
  theme(plot.title = element_text(hjust = 0.5))

g2 <- g2 + ggtitle("Preseq library complexity estimate - overall with CI") +
  labs(x="Total number of reads", y="Expected distinct reads") +
  theme(plot.title = element_text(hjust = 0.5))

pdf("preseq_estimates_long.pdf")
  print(g1)
  print(g2)
dev.off()

