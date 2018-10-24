#!/bin/bash
#
# Prepare RSEM alignment index
#
# Requires RSEM, pigz
#

############################################################################################
### Variables
GENOME=/storage/brno2/home/opplatek/genomes/human/ensembl87/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
GTF=/storage/brno2/home/opplatek/genomes/human/ensembl87/Homo_sapiens.GRCh38.87.gtf.gz
INDEX_OUTPUT_DIR=/storage/brno2/home/opplatek/genomes/human/ensembl87/RSEM_index

THREADS=$PBS_NUM_PPN

# Binaries
RSEM_BIN=/storage/brno2/home/opplatek/tools/RSEM-1.3.0

############################################################################################
### Copy inputs
cp $GTF $SCRATCH/ &
cp $GENOME $SCRATCH/ &

wait

####################################################################################################
### Genome and annotation preparation

GENOME=$(basename $GENOME)
GTF=$(basename $GTF)

unpigz -p $THREADS $GENOME
unpigz -p $THREADS $GTF

GENOME=${GENOME%.gz*}
GTF=${GTF%.gz*}
GEN_DIR=${GENOME%.fa*}

############################################################################################
### Indexing - RUN ONLY ONCE FOR ONE GENOME AND RSEM VERSION 
mkdir $SCRATCH/rsem_index
$RSEM_BIN/rsem-prepare-reference --gtf $SCRATCH/$GTF $GENOME $SCRATCH/rsem_index/$GEN_DIR

############################################################################################
### Copying
mkdir -p $INDEX_OUTPUT_DIR
cp $SCRATCH/rsem_index/$GEN_DIR* $INDEX_OUTPUT_DIR/

rm -r $SCRATCH/*
