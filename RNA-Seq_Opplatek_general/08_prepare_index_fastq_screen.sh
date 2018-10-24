#!/bin/bash
#
# Prepare rRNA and tRNA sequences for fastq_screen 
# Downloads annotation from NCBI (RefSeq) .gff, extracts selected gene biotypes and extracts .fasta sequences for fastq_screen
#
# Requires: cufflings (gffread)

NCBI_GENOME="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz"
NCBI_GFF="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff.gz"

OUTPUT_DIR=/storage/brno2/home/opplatek/genomes/drosophila/fastq_screen/

THREADS=$PBS_NUM_PPN

module add cufflinks-2.2.1
GFFREAD=$(which gffread)
module add bowtie2-2.3.0
BOWTIE2_BUILD=$(which bowtie2-build)
module add bowtie-1.0.0
BOWTIE_BUILD=$(which bowtie-build)

##########################################################################################################################################
cd $SCRATCH

wget $NCBI_GENOME &
wget $NCBI_GFF

unpigz -c -p $THREADS $(basename $NCBI_GFF) | grep "gbkey=rRNA" | grep -v "ribosomal RNA protein" > $(basename $NCBI_GFF .gz).rRNA.gff
unpigz -c -p $THREADS $(basename $NCBI_GFF) | grep "gbkey=tRNA" > $(basename $NCBI_GFF .gz).tRNA.gff

wait

unpigz -p $THREADS $(basename $NCBI_GENOME)

$GFFREAD $(basename $NCBI_GFF .gz).rRNA.gff -g $(basename $NCBI_GENOME .gz) -w $(basename $NCBI_GFF .gz).rRNA.fasta &
$GFFREAD $(basename $NCBI_GFF .gz).tRNA.gff -g $(basename $NCBI_GENOME .gz) -w $(basename $NCBI_GFF .gz).tRNA.fasta

wait

##########################################################################################################################################
mkdir -p $SCRATCH/index

# Make bowtie2 index
$BOWTIE2_BUILD $(basename $NCBI_GENOME .gz) $SCRATCH/index/$(basename $NCBI_GENOME .gz)
$BOWTIE2_BUILD $(basename $NCBI_GFF .gz).rRNA.fasta $SCRATCH/index/$(basename $NCBI_GFF .gz).rRNA.fasta
$BOWTIE2_BUILD $(basename $NCBI_GFF .gz).tRNA.fasta $SCRATCH/index/$(basename $NCBI_GFF .gz).tRNA.fasta

# Make bowtie1 index
$BOWTIE_BUILD $(basename $NCBI_GENOME .gz) $SCRATCH/index/$(basename $NCBI_GENOME .gz)
$BOWTIE_BUILD $(basename $NCBI_GFF .gz).rRNA.fasta $SCRATCH/index/$(basename $NCBI_GFF .gz).rRNA.fasta
$BOWTIE_BUILD $(basename $NCBI_GFF .gz).tRNA.fasta $SCRATCH/index/$(basename $NCBI_GFF .gz).tRNA.fasta

##########################################################################################################################################
# Copy outputs
mkdir $OUTPUT_DIR

rm $(basename $NCBI_GFF)
rm $(basename $NCBI_GENOME .gz)

cp -r $SCRATCH/* $OUTPUT_DIR

rm -r $SCRATCH/*