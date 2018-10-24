#!/bin/bash
#PBS -l nodes=1:ppn=30
#PBS -l mem=200gb
#PBS -k oe
#PBS -N 3v05_prepare_STAR_index_Scerevisiae_STAR2.6.0c
#PBS -M tomas.skalicky@seznam.cz
#PBS -m abe
#
#PARU KRTECEK server is using TORQUE scheduling system !!!
#
## initialize the required application
#
# Prepare STAR alignment index
#
# Requires STAR, RSEM, pigz, bowtie or alernatively bbmap aligner
#
## initialize the required application
# module add rsem-1.2.8 # Metacentrum has only an old version that doesn't support STAR aligner
# module add star-2.5.2b # METAcentrum has only old v2.5.2b which has problems with 2nd pass mapping
# module add bowtie-1.0.0
# module add bowtie2-2.3.0
############################################################################################
### Variables
GENOME="/home/users/tskalicky/CEITEC/genomes/Saccharomyces_cerevisiae/Scerevisiae_R64_genomic.fna.gz"
DATADIR="/home/users/tskalicky/CEITEC/genomes/Saccharomyces_cerevisiae"
INDEX_OUTPUT_DIR="/home/users/tskalicky/CEITEC/genomes/Saccharomyces_cerevisiae/STAR_index"
GTF_MODIF="/home/users/tskalicky/CEITEC/genomes/Saccharomyces_cerevisiae/Scerevisiae_R64_genomic.gff.gz"
#
THREADS=$PBS_NUM_PPN # Number of threads to use
# Adding RSEM binaries into the PATH
# export PATH="/storage/brno3-cerit/home/tskalicky/tools/RSEM-1.3.0/bin:$PATH"
# export PATH="/storage/brno3-cerit/home/tskalicky/tools/STAR-2.6.0a/bin/Linux_x86_64_static:$PATH" # this version is not working
export PATH="/home/users/tskalicky/anaconda2/bin"
#
# Binaries
STAR_BIN="$(which STAR)"
# RSEM-PREP=$(which rsem-prepare-reference)
STAR=$(which STAR)
# BOWTIE2="$(which bowtie2)"
# Check the tools versions
# echo $RSEM-PREP
echo $STAR
# echo $BOWTIE2
#
# Prepare STAR alignment index
#
# Requires STAR, pigz
#
# STAR options help:
#--------------------
# "For small genomes, the parameter --genomeSAindexNbases must to be scaled down, with a typical
#       value of min(14, log2(GenomeLength)/2 - 1). For example, for 1 megaBase genome, this is equal
#       to 9 (8.965784), for 100 kiloBase genome, this is equal to 7."
############################################################################################
# Read length for sjdbOverhang; --sjdbOverhang 100 should work fine for most of the data, but more specific setting based on the real read length should be more sensitive https://groups.google.com/forum/#!msg/rna-star/h9oh10UlvhI/BfSPGivUHmsJ
RD_LENGTH=79 # Read length from the sequencing. Illumina sometimes reports N+1 http://seqanswers.com/forums/archive/index.php/t-31154.html; in case you change this value uncomment next line as well
#RD_LENGTH=$[$RD_LENGTH-$TRIM_LEFT-$TRIM_RIGHT-1] # Modified read length for the index creation and mapping - should be read length -1
#RD_LENGTH changed to 124 because FTO PE libs after trimminf are in range of 17-125 bases
############################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
# cp -av $GENOME $INDEX_OUTPUT_DIR
# cp -av $GTF_MODIF $INDEX_OUTPUT_DIR
cd $INDEX_OUTPUT_DIR

if [ ! -d "$INDEX_OUTPUT_DIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "INDEX_OUTPUT_DIR path is:" $INDEX_OUTPUT_DIR
echo "Following reads were copied to scratch:"
ls -Rc1
####################################################################################################
### Genome and annotation preparation
# "basename" will print file NAME with any leading directory components removed
# "dirname" will print DIRECTORY PATH without filename
# GENOME_NAME=$(basename $GENOME)
# GTF_NAME=$(basename $GTF_MODIF)
# 
# unpigz -p $THREADS $GENOME_NAME
# unpigz -v -p $THREADS -d $GTF_NAME
# 
# GENOME2=${GENOME_NAME%.*} # use when files in *.gz
# GTF2=${GTF_NAME%.*} # use when files in *.gz
# GEN_DIR=${GENOME_NAME%.*.*}"_RD_LENGTH_79"
GEN_DIR="Scerevisiae_R64_RD_LENGTH_79_STAR2.6c"
# ENSEMBL changed the gtf format and RSEM has problem with missing "transcript_id" patern. We need to filter the *.gtf file using awk for example

#awk '$3 == "exon"' $GTF_NAME > $GTF2_modif.gtf

#GTF3=$GTF2_modif.gtf

# # Creating gtf annotation file including predicted tRNAs
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

############################################################################################
### Indexing - RUN ONLY ONCE FOR ONE GENOME AND STAR VERSION AND READ LENGTH OR USE DEFAULT INDEX
mkdir "/home/users/tskalicky/CEITEC/genomes/Saccharomyces_cerevisiae/STAR_index/Scerevisiae_R64_RD_LENGTH_79_STAR2.6c"
# for very small genomes the --genomeSAindexNbases has to ge set = log2(GenomeLength)/2-1
# when using GFF instead of GTF annotation the "--sjdbGTFtagExonParentTranscript Parent" has to be set
# --sjdbOverhang 100 should work fine for most of the data, but more specific setting 
# based on the real read length should be more sensitive https://groups.google.com/forum/#!msg/rna-star/h9oh10UlvhI/BfSPGivUHmsJ
STAR --runMode genomeGenerate --runThreadN $THREADS \
--genomeDir "$INDEX_OUTPUT_DIR/$GEN_DIR" --genomeFastaFiles "$DATADIR/Scerevisiae_R64_genomic.fna" \
--sjdbGTFfile "$DATADIR/Scerevisiae_R64_genomic.gff" --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang $RD_LENGTH \
 --genomeSAindexNbases 11

############################################################################################
### Copy data from scratch back to home dir and clean scratch
# rm $INDEX_OUTPUT_DIR/*.fna $INDEX_OUTPUT_DIR/*.gff
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
