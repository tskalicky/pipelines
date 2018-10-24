#
# Plot rRNA results from fastq_screen in one plot
# Uses multiQC summary of the results
#

INPUT_DIR<-"."
OUTPUT_DIR<-INPUT_DIR
INPUT_FILE<-paste(INPUT_DIR, "multiqc_data", "multiqc_fastq_screen.txt", sep="/")

#library("rio")
#rRNA<-rio::import(INPUT_FILE, format = "tsv")
rRNA<-read.table(INPUT_FILE, sep="\t", stringsAsFactors = F, header=T)

rRNA$Sample<-gsub("_001_trim_screen", "", rRNA$Sample)
rRNA$Sample<-gsub("_trim_screen", "", rRNA$Sample)

rownames(rRNA)<-rRNA$Sample

pdf(paste(OUTPUT_DIR, "rRNA_estimate_fastqScreen.pdf", sep="/"))
  par(mar = c(10,4,4,2) + 0.1) # default is c(5,4,4,2) + 0.1; (1=bottom, 2=left, 3=top and 4=right)
  rRNA_perc<-rRNA[, grep("rRNA.*percentage", colnames(rRNA))]
  barplot(rRNA_perc, ylim=c(0,max(rRNA_perc)*1.3), ylab="rRNA content (%)", 
        main="rRNA estimate", names.arg=rownames(rRNA), las=2)
dev.off()
