#!/bin/bash
#PBS -l nodes=1:ppn=30
#PBS -l mem=200gb
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -N 3v06_prep_STAR_index_human_STAR2.6c
#PBS -M tomas.skalicky@seznam.cz
#PBS -m abe
#
# PARU KRTECEK server is using TORQUE scheduling system !!!
#
## initialize the required application
#
# Prepare STAR alignment index
#
# Requires STAR, RSEM, pigz, bowtie or alernatively bbmap aligner
#
# STAR options help:
#--------------------
# "For small genomes, the parameter --genomeSAindexNbases must to be scaled down, with a typical
#       value of min(14, log2(GenomeLength)/2 - 1). For example, for 1 megaBase genome, this is equal
#       to 9 (8.965784), for 100 kiloBase genome, this is equal to 7."
############################################################################################
# Read length for sjdbOverhang; --sjdbOverhang 100 should work fine for most of the data, but more specific setting based on the real read length should be more sensitive https://groups.google.com/forum/#!msg/rna-star/h9oh10UlvhI/BfSPGivUHmsJ
RD_LENGTH=100 # Read length from the sequencing. Illumina sometimes reports N+1 http://seqanswers.com/forums/archive/index.php/t-31154.html; in case you change this value uncomment next line as well
#RD_LENGTH=$[$RD_LENGTH-$TRIM_LEFT-$TRIM_RIGHT-1] # Modified read length for the index creation and mapping - should be read length -1
#RD_LENGTH changed to 124 because FTO PE libs after trimminf are in range of 17-125 bases
############################################################################################
## initialize the required application
# Adding STAR binaries into the PATH
export PATH="/home/users/tskalicky/anaconda2/bin:$PATH"
#
# Binaries
# RSEM-PREP=$(which rsem-prepare-reference)
STAR=$(which STAR)
# BOWTIE2="$(which bowtie2)"
# Check the tools versions
# echo $RSEM-PREP
echo $STAR
# echo $BOWTIE2
############################################################################################
### Variables
THREADS=$PBS_NUM_PPN # Number of threads to use
#
GENOME="/home/users/tskalicky/CEITEC/genomes/human/ensembl91/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
DATADIR="/home/users/tskalicky/CEITEC/genomes/human/ensembl91"
INDEX_OUTPUT_DIR="/home/users/tskalicky/CEITEC/genomes/human/ensembl91/STAR_index"
GTF="/home/users/tskalicky/CEITEC/genomes/human/ensembl91/GRCh38.91.gtf.gz"
GEN_DIR="GRCh38_RD_LENGTH_100_STAR2.6c"
#
#############################################################################################
# copy input data using local storage
# use cp -avr when copying directories
mkdir $INDEX_OUTPUT_DIR/$GEN_DIR
cd $DATADIR
cp -av $GENOME $INDEX_OUTPUT_DIR
cp -av $GTF $INDEX_OUTPUT_DIR
#
####################################################################################################
### Genome and annotation preparation
# "basename" will print file NAME with any leading directory components removed
# "dirname" will print DIRECTORY PATH without filename
GENOME_NAME=$(basename $GENOME)
GTF_NAME=$(basename $GTF)
# 
gzip -vd $GENOME_NAME
gzip -vd $GTF_NAME
#
# UNPIGZ is not working on Krtecek...
# unpigz -p $THREADS $GENOME_NAME
# unpigz -v -p $THREADS -d $GTF_NAME
# 
GENOME2=${GENOME_NAME%.*} # use when files in *.gz
GTF2=${GTF_NAME%.*} # use when files in *.gz
# GEN_DIR=${GENOME_NAME%.*.*}"_RD_LENGTH_79"
# 
## FOR RSEM ONLY: 
# ENSEMBL changed the gtf format and RSEM has problem with missing "transcript_id" patern. 
# We need to filter the *.gtf file using awk for example

#awk '$3 == "exon"' $GTF_NAME > $GTF2_modif.gtf

#GTF3=$GTF2_modif.gtf

### Creating gtf annotation file including predicted tRNAs:
# # Need to modify the tRNA.gtf file first in order that it has "exon" value in 3rd field and each line has to have "transcript_id" and "gene_id" values
# # will first write all comments (lines beginning with #) to file_sorted.gtf, 
# # then all other lines sorted first by the first column, alphanumerically, 
# # second by the fourth column, numerically. This works as well for files that do not contain comments.
# cat Homo_sapiens.GRCh38.91.sorted.gtf gencode.v27.tRNAs_modif_for_star_index.gtf > GRCh38.91_tRNA_added.gtf
# 
# (grep "^#" GRCh38.91_tRNA_added.gtf; grep -v "^#" GRCh38.91_tRNA_added.gtf | sort -k1,1 -k4,4n) > GRCh38.91_tRNA_sorted.gtf
# 
# # for compressed gtf files use:
# # (zgrep "^#" file.gtf; zgrep -v "^#" file.gtf | sort -k1,1 -k4,4n) | gzip > file_sorted.gtf.gz
#
### COMMANDS
#
############################################################################################
### Indexing - RUN ONLY ONCE FOR ONE GENOME AND STAR VERSION AND READ LENGTH OR USE DEFAULT INDEX
# for very small genomes the --genomeSAindexNbases has to ge set = log2(GenomeLength)/2-1
# when using GFF instead of GTF annotation the "--sjdbGTFtagExonParentTranscript Parent" has to be set
# --sjdbOverhang 100 should work fine for most of the data, but more specific setting 
# based on the real read length should be more sensitive https://groups.google.com/forum/#!msg/rna-star/h9oh10UlvhI/BfSPGivUHmsJ
# For small genomes like Scerevisiae use this command:
# STAR --runMode genomeGenerate --runThreadN $THREADS \
# --genomeDir "$INDEX_OUTPUT_DIR/$GEN_DIR" --genomeFastaFiles "$DATADIR/Scerevisiae_R64_genomic.fna" \
# --sjdbGTFfile "$DATADIR/Scerevisiae_R64_genomic.gff" --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang $RD_LENGTH \
#  --genomeSAindexNbases 11
#
cd $INDEX_OUTPUT_DIR
STAR --runMode genomeGenerate --runThreadN $THREADS --limitGenomeGenerateRAM 150000000000 \
--genomeDir "$INDEX_OUTPUT_DIR/$GEN_DIR" --genomeFastaFiles "$GENOME2" \
--sjdbGTFfile "$DATADIR/GTF2" --sjdbOverhang $RD_LENGTH # --sjdbOverhang 100 should work fine for most of the data, but more specific setting based on the real read length should be more sensitive https://groups.google.com/forum/#!msg/rna-star/h9oh10UlvhI/BfSPGivUHmsJ
#
############################################################################################
### Copy data from scratch back to home dir and clean scratch
rm -v $INDEX_OUTPUT_DIR/*.fa $INDEX_OUTPUT_DIR/*.gtf
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
