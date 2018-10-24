#!/bin/bash
#PBS -l walltime=24:0:0 
#PBS -q uv@wagap-pro.cerit-sc.cz 
#PBS -l select=1:ncpus=20:mem=150gb:scratch_local=80gb
#PBS -j oe
#PBS -N 03v3_STAR_index_RDlength79_RSEM_genome_ref
#
# Prepare STAR alignment index
#
# Requires STAR, RSEM, pigz, bowtie or alernatively STAR aligner
#
## initialize the required application
# module add rsem-1.2.8 # Metacentrum has only an old version that doesn't support STAR aligner
module add star-2.5.2b # METAcentrum has only old v2.5.2b which has problems with 2nd pass mapping
module add bowtie-1.0.0
module add bowtie2-2.3.0
############################################################################################
### Variables
GENOME="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/GRCh38.91.dna.primary_assembly.fa.gz"
DATADIR="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91"
INDEX_OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/STAR_index"
GTF="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/Homo_sapiens.GRCh38.91.sorted.gtf.gz"
#
THREADS=$PBS_NUM_PPN # Number of threads to use
# Adding RSEM binaries into the PATH
# export PATH="/storage/brno3-cerit/home/tskalicky/tools/RSEM-1.3.0/bin:$PATH"
# export PATH="/storage/brno3-cerit/home/tskalicky/tools/STAR-2.6.0a/bin/Linux_x86_64_static:$PATH" # this version is not working
#
# Binaries
STAR_BIN="$(which STAR)"
RSEM-PREP=$(which rsem-prepare-reference)
STAR=$(which STAR)
BOWTIE2="$(which bowtie2)"
# Check the tools versions
echo $RSEM-PREP
echo $STAR
echo $BOWTIE2
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
RD_LENGTH=100 # Read length from the sequencing. Illumina sometimes reports N+1 http://seqanswers.com/forums/archive/index.php/t-31154.html; in case you change this value uncomment next line as well
#RD_LENGTH=$[$RD_LENGTH-$TRIM_LEFT-$TRIM_RIGHT-1] # Modified read length for the index creation and mapping - should be read length -1
#RD_LENGTH changed to 124 because FTO PE libs after trimminf are in range of 17-125 bases
############################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
cd $DATADIR
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
cp -av $GENOME $SCRATCHDIR
cp -av $GTF $SCRATCHDIR
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR
echo "Following reads were copied to scratch:"
ls -Rc1
echo "Following reads were copied to scratch:" >> /storage/brno3-cerit/home/tskalicky/debug_03v3_STAR_index_RDlength79_RSEM_genome_ref.log # will be removed
ls -Rc1 >> /storage/brno3-cerit/home/tskalicky/debug_03v3_STAR_index_RDlength79_RSEM_genome_ref.log # will be removed
####################################################################################################
### Genome and annotation preparation
# "basename" will print file NAME with any leading directory components removed
# "dirname" will print DIRECTORY PATH without filename
GENOME_NAME=$(basename $GENOME)
GTF_NAME=$(basename $GTF)

unpigz -p $THREADS $GENOME_NAME
unpigz -v -p $THREADS -d $GTF_NAME

GENOME2=${GENOME_NAME%.*} # use when files in *.gz
GTF2=${GTF_NAME%.*} # use when files in *.gz
# GEN_DIR=${GENOME_NAME%.*.*}"_RD_LENGTH_79"
GEN_DIR="GRCh38_with_RD_LENGTH_100"
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
mkdir $SCRATCHDIR/$GEN_DIR

STAR --runMode genomeGenerate --runThreadN $THREADS --limitGenomeGenerateRAM 150000000000 \
--genomeDir $SCRATCHDIR/$GEN_DIR --genomeFastaFiles $SCRATCHDIR/$GENOME2 \
--sjdbGTFfile $SCRATCHDIR/$GTF2 --sjdbOverhang $RD_LENGTH # --sjdbOverhang 100 should work fine for most of the data, but more specific setting based on the real read length should be more sensitive https://groups.google.com/forum/#!msg/rna-star/h9oh10UlvhI/BfSPGivUHmsJ

############################################################################################
### Copy data from scratch back to home dir and clean scratch
mkdir -p $INDEX_OUTPUT_DIR
cp -avr $SCRATCHDIR/$GEN_DIR/ $INDEX_OUTPUT_DIR/ || export CLEAN_SCRATCH=false
rm -r $SCRATCHDIR/*
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
