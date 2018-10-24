#!/bin/bash
#PBS -l nodes=1:ppn=20
#PBS -l mem=100gb
#PBS -l walltime=96:00:00
#PBS -k oe
#PBS -N star_prep_STAR_index_OAT_U12
#PBS -M tomas.skalicky@seznam.cz
#PBS -m abe
#
#PARU KRTECEK server is using TORQUE scheduling system !!!
#
# Prepare STAR alignment index
#
# Requires STAR, RSEM, pigz, bowtie or alernatively STAR aligner
#
############################################################################################
### Variables
export PATH="$PATH:/home/users/tskalicky/anaconda2/bin"
echo "PATH is:"
echo $PATH | tr " " "\n" | nl
#
DATADIR="/home/users/tskalicky/CEITEC/genomes/human/ensembl91"
GENOME="/home/users/tskalicky/CEITEC/genomes/human/OAT_U12_only.fa.gz"
GTF="/home/users/tskalicky/CEITEC/genomes/human/OAT_U12_only.gtf.gz"
INDEX_OUTPUT_DIR="/home/users/tskalicky/CEITEC/genomes/human/ensembl91/STAR_index"

THREADS=$PBS_NUM_PPN # Number of threads to use
# Binaries
STAR=$(which STAR)
BOWTIE2="$(which bowtie2)"
# Check the tools versions
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
#RD_LENGTH changed to 124 because FTO PE libs after trimming are in range of 17-125 bases
#RD_LENGTH changed to 79 because FTO PE libs after trimming are in range of 17-80 bases
############################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
cd $DATADIR
cp -av $GENOME $GTF $INDEX_OUTPUT_DIR
cd $INDEX_OUTPUT_DIR

if [ ! -d "$INDEX_OUTPUT_DIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "INDEX_OUTPUT_DIR path is:" $INDEX_OUTPUT_DIR
####################################################################################################
### Genome and annotation preparation
# "basename" will print file NAME with any leading directory components removed
# "dirname" will print DIRECTORY PATH without filename
GENOME_NAME=$(basename $GENOME)
GTF_NAME=$(basename $GTF)

gzip -dv $GENOME_NAME
gzip -dv $GTF_NAME

GENOME2=${GENOME_NAME%.*}
GTF2=${GTF_NAME%.*}
GEN_DIR="OAT_RNU12_RD_LENGTH_100_STAR_2.6.1"
# ENSEMBL changed the gtf format and RSEM has problem with missing "transcript_id" patern. We need to filter the *.gtf file using awk for example

#awk '$3 == "exon"' $GTF_NAME > $GTF2_modif.gtf

#GTF3=$GTF2_modif.gtf

############################################################################################
### Indexing - RUN ONLY ONCE FOR ONE GENOME AND STAR VERSION AND READ LENGTH OR USE DEFAULT INDEX
mkdir $INDEX_OUTPUT_DIR/$GEN_DIR

$STAR --runMode genomeGenerate --runThreadN $THREADS --genomeDir $INDEX_OUTPUT_DIR/$GEN_DIR --genomeFastaFiles $INDEX_OUTPUT_DIR/$GENOME2 \
--sjdbGTFfile $INDEX_OUTPUT_DIR/$GTF2 --sjdbOverhang $RD_LENGTH --genomeSAindexNbases 8 # --sjdbOverhang 100 should work fine for most of the data, but more specific setting based on the real read length should be more sensitive https://groups.google.com/forum/#!msg/rna-star/h9oh10UlvhI/BfSPGivUHmsJ

############################################################################################
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
