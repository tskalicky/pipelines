# variables
ref_dir<-"/home/skalda/ownCloud/CEITEC_lab/genomes/human/ensembl91" # Path to the reference GTF
gtf_file<-"Homo_sapiens.GRCh38.91.sorted.gtf" # "Homo_sapiens.GRCh38.91.gtf" # Reference GTF in ref_dir
biotypes<-read.table("/home/skalda/ownCloud/git/pipelines/CLIP_general/annotation/Honza_Oppelt/biotypes_info_ens.txt", header=T, sep="\t", stringsAsFactors = F)
#
featuresToAnalyze<-c("rRNA", "rRNA_pseudogene",); name<-featuresToAnalyze # could be multiple from - "coding", "lnoncoding", "snoncoding", "mnoncoding", "pseudogene", "LRG", "undefined", "no_group". For more info see biotypes table
featuresToAnalyze2<-c(biotypes[biotypes$biotype_group %in% featuresToAnalyze,1])

library("ensembldb")

gtf.file<-file.path(ref_dir, gtf_file)
sqlite_file<-gsub(".gtf", ".sqlite", gtf_file, fixed=T) # Name of the SQL file
sqlite_path <- file.path(ref_dir, sqlite_file)

# This part sometimes does not put sqlite database to the 'path', needs to be transfered manually
if(!file.exists(sqlite_path)) {
  ## generate the SQLite database file
  ensembldb::ensDbFromGtf(gtf=gtf.file, outfile=sqlite_path)
}
# create ensembl database
EnsDb <- ensembldb::EnsDb(sqlite_path)
#ensembldb::listColumns(EnsDb.Hsapiens.v91) # See all available columns
#parsedEnsembl <- ensembldb::genes(EnsDb, filter=list(GeneBiotypeFilter(featuresToAnalyze2))) # Get only genes with specificed biotype(s)
parsedEnsembl <- ensembldb::genes(EnsDb)
parsedEnsembl <- as.data.frame(parsedEnsembl)
parsedEnsembl <- parsedEnsembl[,c("gene_id", "gene_name", "gene_biotype")]
parsedEnsembl[is.na(parsedEnsembl$gene_name), "gene_name"]<-parsedEnsembl[is.na(parsedEnsembl$gene_name), "gene_id"] # Replace missing gene names by gene ids
rownames(parsedEnsembl)<-parsedEnsembl$gene_id

# Get only selected genes for a biotype
mrcounts<-mrcounts[rownames(mrcounts) %in% parsedEnsembl[parsedEnsembl$gene_biotype %in% featuresToAnalyze2, 
                                                         "gene_id"],]