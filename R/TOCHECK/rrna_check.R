# 
# rRNA content check from from featureCounts and normalize
#

rm(list=ls(all=TRUE))

library(DESeq2)
library(ensembldb)

INPUT_DIR="/home/jan/Data/projects/jitka_quantSeq/results/star_selected/gene_counts/featureCounts"
INPUT_FILE="quantSeq.fwdStrand.featureCounts"

OUTPUT_DIR="/home/jan/Projects/jitka_quantSeq/results/qc/rRNA_check/star_selected"

#ensembl_parsed<-read.table("/home/jan/Data/projects/jitka_quantSeq/other/Homo_sapiens.GRCh38.87.parsed.txt", sep="\t", header=T, stringsAsFactors = F, quote="")
gene_counts<-read.table(paste(INPUT_DIR, INPUT_FILE, sep="/"), sep="\t", header=T, stringsAsFactors = F, quote="")

#ensembl_parsed[grep("5S_rRNA", ensembl_parsed$Gene.name),]
# Prepare Ensembl annotation
ref_dir<-"/home/jan/Data/projects/jitka_quantSeq/other"
gtf.file<-file.path(ref_dir, "Homo_sapiens.GRCh38.87.gtf")
sqlite_file<-"Homo_sapiens.GRCh38.87.sqlite"
sqlite_path <- file.path(ref_dir, sqlite_file)

if(!file.exists(sqlite_path)) {
  ## generate the SQLite database file
  ensembldb::ensDbFromGtf(gtf=gtf.file, path=ref_dir, outfile=sqlite_file)
}

EnsDb.Hsapiens.v87 <- ensembldb::EnsDb(sqlite_path)
ens_rRNA <- ensembldb::genes(EnsDb.Hsapiens.v87, filter=list(GenebiotypeFilter(c('rRNA_pseudogene', 'rRNA', 'Mt_rRNA'))))
ens_rRNA <- as.data.frame(ens_rRNA)

rownames(gene_counts)<-gene_counts$Geneid
gene_counts<-gene_counts[, c(grep(".bam$", colnames(gene_counts)))]
colnames(gene_counts)<-gsub("\\_.*","", colnames(gene_counts))

mrcounts<-gene_counts

#mrcounts <- mrcounts[ rowSums(mrcounts) > 1, ]

# Count rRNAs
mrcounts_rRNA<-mrcounts[rownames(mrcounts) %in% rownames(ens_rRNA),]

dir.create(OUTPUT_DIR, recursive = T)

pdf(paste0(OUTPUT_DIR, "/rRNA_estimate.pdf"))
barplot(colSums(mrcounts_rRNA)/(colSums(mrcounts)/100), ylim=c(0, max(colSums(mrcounts_rRNA)/(colSums(mrcounts)/100))*1.3), main="rRNA content estimate (rRNA, rRNA pseudogene, Mt rRNA") # Get rRNA percentage
dev.off()
