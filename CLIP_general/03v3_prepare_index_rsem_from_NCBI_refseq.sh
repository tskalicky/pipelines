#!/bin/bash
#PBS -l walltime=24:0:0 
#PBS -q uv@wagap-pro.cerit-sc.cz 
#PBS -l select=1:ncpus=1:mem=40gb:scratch_local=80gb
#PBS -j oe
#PBS -N 03v3_RSEM_index_from_NCBI_refseq
#
# Prepare RSEM alignment index
#
# Requires RSEM, pigz, bowtie or alernatively STAR aligner
#
## initialize the required application
module add rsem-1.2.8
module add star-2.5.2b
module add bowtie-1.0.0
module add bowtie2-2.3.0
############################################################################################
### Variables
DATADIR="/storage/brno3-cerit/home/tskalicky/genomes/human"
GENOME="/storage/brno3-cerit/home/tskalicky/genomes/human/Homo_sapiens_GCF_000001405.37_GRCh38.p11_genomic.fna.gz"
GFF3="/storage/brno3-cerit/home/tskalicky/genomes/human/Homo_sapiens_GCF_000001405.37_GRCh38.p11_genomic.gff.gz"
INDEX_OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/genomes/human/RSEM_index"

THREADS=$PBS_NUM_PPN

# Binaries
RSEM_BIN="/storage/brno3-cerit/home/tskalicky/tools/RSEM-1.3.0/bin"

############################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
cd $DATADIR
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
cp -av $GENOME $GTF Homo_sapiens.GRCh38.91_modif.gtf $SCRATCHDIR
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR

####################################################################################################
### Genome and annotation preparation
# "basename" will print file NAME with any leading directory components removed
# "dirname" will print DIRECTORY PATH without filename
GENOME_NAME=$(basename $GENOME)
GTF_NAME=$(basename $GTF)

unpigz -p $THREADS $GENOME_NAME
unpigz -p $THREADS $GTF_NAME

GENOME2=${GENOME_NAME%.*}
GFF3=${GTF_NAME%.*}
GEN_DIR=${GENOME_NAME%.fa*}
# ENSEMBL changed the gtf format and RSEM has problem with missing "transcript_id" patern. We need to filter the *.gtf file using awk for example

#awk '$3 == "exon"' $GTF_NAME > $GFF3_modif.gtf

#GTF3=$GFF3_modif.gtf
## Extract primary assembly for RSEM indexing
rsem-refseq-extract-primary-assembly Homo_sapiens_GCF_000001405.37_GRCh38.p11_genomic.fna Homo_sapiens_GCF_000001405.37_GRCh38.p11_genomic.primary_assembly.fna
############################################################################################
### Indexing - RUN ONLY ONCE FOR ONE GENOME AND RSEM VERSION 
mkdir $SCRATCHDIR/rsem_index
#rsem-prepare-reference --gtf $SCRATCHDIR/$GTF3 --bowtie --bowtie2 $GENOME2 $SCRATCHDIR/rsem_index/$GEN_DIR
rsem-prepare-reference --gff3 $SCRATCHDIR/Homo_sapiens_GCF_000001405.37_GRCh38.p11_genomic.gff \
               --trusted-sources BestRefSeq,Curated\ Genomic \
               --bowtie2 $GENOME2 $SCRATCHDIR/rsem_index/$GEN_DIR
############################################################################################
### Copy data from scratch back to home dir and clean scratch
mkdir -p $INDEX_OUTPUT_DIR
cp -avr $SCRATCHDIR/rsem_index/ $INDEX_OUTPUT_DIR/ || export CLEAN_SCRATCH=false
rm -r $SCRATCHDIR/*
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"


