#!/bin/bash
#PBS -l walltime=24:0:0 
#PBS -q default@wagap-pro.cerit-sc.cz 
#PBS -l select=1:ncpus=20:mem=50gb:scratch_local=20gb
#PBS -j oe
#PBS -N 03v045_STAR_index_RDlength100
#
# Prepare STAR alignment index
#
# Requires STAR, RSEM, pigz, bowtie or alernatively STAR aligner
#
## initialize the required application
module add rsem-1.2.8
module add star-2.5.2b
module add bowtie-1.0.0
module add bowtie2-2.3.0
############################################################################################
### Variables
DATADIR="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91"
GENOME="/storage/brno3-cerit/home/tskalicky/genomes/human/OAT_U12/OAT_U12_only.fa.gz"
# GTF="/storage/brno3-cerit/home/tskalicky/genomes/human/OAT_U12/OAT_U12_mapping.gtf.gz"
INDEX_OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/STAR_index"

THREADS=$PBS_NUM_PPN # Number of threads to use
# Binaries
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
#RD_LENGTH changed to 124 because FTO PE libs after trimming are in range of 17-125 bases
#RD_LENGTH changed to 79 because FTO PE libs after trimming are in range of 17-80 bases
############################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
cd $DATADIR
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
cp -av $GENOME $SCRATCHDIR
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR
echo "Following reads were copied to scratch:"
ls -Rc1
####################################################################################################
### Genome and annotation preparation
# "basename" will print file NAME with any leading directory components removed
# "dirname" will print DIRECTORY PATH without filename
GENOME_NAME=$(basename $GENOME)
# GTF_NAME=$(basename $GTF)

unpigz -p $THREADS $GENOME_NAME
# unpigz -p $THREADS $GTF_NAME

GENOME2=${GENOME_NAME%.*}
# GTF2=${GTF_NAME%.*}
GEN_DIR="OAT_RNU12_RD_LENGTH_100_NO_ANNOT_STAR2.5.2"
# ENSEMBL changed the gtf format and RSEM has problem with missing "transcript_id" patern. We need to filter the *.gtf file using awk for example

#awk '$3 == "exon"' $GTF_NAME > $GTF2_modif.gtf

#GTF3=$GTF2_modif.gtf

############################################################################################
### Indexing - RUN ONLY ONCE FOR ONE GENOME AND STAR VERSION AND READ LENGTH OR USE DEFAULT INDEX
mkdir $SCRATCHDIR/$GEN_DIR

STAR --runMode genomeGenerate --runThreadN $THREADS --genomeDir $SCRATCHDIR/$GEN_DIR --genomeFastaFiles $SCRATCHDIR/$GENOME2 \
--genomeSAindexNbases 8 # --sjdbOverhang 100 should work fine for most of the data, but more specific setting based on the real read length should be more sensitive https://groups.google.com/forum/#!msg/rna-star/h9oh10UlvhI/BfSPGivUHmsJ

############################################################################################
### Copy data from scratch back to home dir and clean scratch
mkdir -p $INDEX_OUTPUT_DIR
cp -avr $SCRATCHDIR/$GEN_DIR/ $INDEX_OUTPUT_DIR/ || export CLEAN_SCRATCH=false
rm -r $SCRATCHDIR/*
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
