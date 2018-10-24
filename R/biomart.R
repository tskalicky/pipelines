#
# Add information from Biomatr to the DE results (DESeq2, edgeR)
# This one is based on Ensembl results but works with any Biomart compatible results
#

rm(list=ls(all=TRUE))

INPUT_TABLE<-"/home/jan/Projects/anil_quantSeq/G1S/results/BCLAF1/ENCODE/DE/HepG2_shBCLAF1vsHepG2_control_HepG2/prot_coding/DESeq2.csv"
#GO_INTEREST<-"cell cycle" # In what are we interested?
GO_INTEREST<-"G1/S" #
GO_MANUAL<-c("GO:0070318", "GO:0070317", "GO:0045023") # Do you want to add some GO groups manually? Can be empty

res2<-read.table(INPUT_TABLE, sep="\t", header=T, row.names=1, stringsAsFactors = F, quote = "\"")
res2[,c("goslim_goa_accession", "goslim_goa_description")]<-NULL # Remove columns we don't need

# Add Gene Ontology (GO) information
library("biomaRt") # https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html#using-archived-versions-of-ensembl or http://www.ensembl.org/info/data/biomart/biomart_r_package.html; atnoher option would be ftp://cran.r-project.org/pub/R/web/packages/biomartr/vignettes/Functional_Annotation.html
library("dplyr")
#listMarts() # Find gene annotation you used
listMarts(host='dec2016.archive.ensembl.org')
#ensembl=useMart("ensembl")
ensembl87=useMart(host='dec2016.archive.ensembl.org', 
                  biomart='ENSEMBL_MART_ENSEMBL', 
                  dataset='hsapiens_gene_ensembl') # Load it

#listAttributes(ensembl87)
#listFilters(ensembl87) 

# Get full GO
goids <- getBM(attributes = c('ensembl_gene_id', 'go_id', 'name_1006', 'definition_1006'),
              filters = 'ensembl_gene_id',
              values = rownames(res2),
              mart = ensembl87)
head(goids)
select_goids<-NULL
select_goids<-unique(goids[grep(GO_INTEREST, goids$name_1006, ignore.case = T), ]) # Get by name
select_goids<-unique(rbind(select_goids, unique(goids[grep(GO_INTEREST, goids$definition_1006, ignore.case = T), ]))) # Get by description
select_goids<-rbind(select_goids, dplyr::filter(goids, go_id %in% GO_MANUAL))

# Get only slimmed GO
#goids_slim <- getBM(attributes = c('ensembl_gene_id', 'goslim_goa_accession', 'goslim_goa_description'), 
#                    filters = 'ensembl_gene_id', 
#                    values = rownames(res2), 
#                    mart = ensembl87)
# head(goids_slim)
#goids_slim[grep(GO_INTEREST, goids_slim$goslim_goa_description, ignore.case = T, fixed = T), ])

### Set chosen GO for merging
go_to_merge<-select_goids

# Merge with the main table
#rownames(go_to_merge)<-go_to_merge$ensembl_gene_id
#go_to_merge$ensembl_gene_id<-NULL
res2_bckp<-res2 # Make backup
res2$tmp_id<-1:nrow(res2)# Temporary ID to sort by, original order is destroyed by merging
res2$ensembl_gene_id<-rownames(res2) # Temporary copy rownames to new column
res2<-merge(res2, go_to_merge, by="ensembl_gene_id", all=TRUE)
res2<-arrange(res2, tmp_id) # Sort by temporary ID; discrads row.names!
res2$tmp_id<-NULL
#rownames(res2)<-res2$Row.names
#res2$Row.names<-NULL

# Merge all GO for one gene to one cell and separate them with ";". Do the same for other columns as well
head(res2)

res2_mergedGO<-as.data.frame(matrix(ncol=ncol(res2)))
colnames(res2_mergedGO)<-colnames(res2)

for(uniqGene in unique(res2$ensembl_gene_id)){
  res2_tmp<-filter(res2, ensembl_gene_id==uniqGene)

  if(nrow(res2_tmp)>1){ # If there are more rows than one (multiple GO annotations)
    res2_tmp2<-res2_tmp[1,] # Copy only first row
    res2_tmp2[,c("go_id", "name_1006", "definition_1006")]<-"NA" # Delete values in GO columns (those are possibly there multiple times (go_id, name_1006, definition_1006) but keep the rest
    res2_tmp2$go_id<-paste(res2_tmp$go_id, collapse = "; ") # Merge cells to one
    res2_tmp2$name_1006<-paste(res2_tmp$name_1006, collapse = "; ") # Merge cells to one
    res2_tmp2$definition_1006<-paste(res2_tmp$definition_1006, collapse = "; ") # Merge cells to one
    res2_tmp<-res2_tmp2
  }
  
  res2_mergedGO<-rbind(res2_mergedGO, res2_tmp) # Merge with main table
}

res2_mergedGO<-res2_mergedGO[!(rowSums(is.na(res2_mergedGO))==ncol(res2_mergedGO)),] # Remove initial empty row
res2_mergedGO<-arrange(res2_mergedGO, padj, pvalue, log2FoldChange)

#write.table(x = res2, file = paste0(gsub(".csv", "_GO", INPUT_TABLE), "_", GO_INTEREST, ".csv"), sep = ",", col.names = T, row.names=F)
write.table(x = res2_mergedGO, file = paste0(gsub(".csv", "_GO", INPUT_TABLE), ".csv"), sep = ",", col.names = T, 
            row.names=F)

res2<-res2_bckp # Copy it back

# Get all Offsprings for selected GO https://www.biostars.org/p/13746/
#source("http://bioconductor.org/biocLite.R")
#biocLite("GO.db")
#library("GO.db")
#offspring = GOMFOFFSPRING[["GO:0015036"]]
