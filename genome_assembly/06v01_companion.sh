#!/bin/bash
###########################################################################################
# This script will be running on Trypka server 
###########################################################################################
#
## initialize the required application
#
# Requires docker, companion
#
############################################################################################
### Variables
INPUT="/home/tomas/Hinxton/CANU_assembly/wallacemonas-25m-erate-0.075/Wallacemonas.contigs.fasta"
OUTPUT_DIR="/home/tomas/Hinxton/annotation_comppanion/NEW/Mbr04_Wallacemonas"
# CONTIGS="/home/users/tskalicky/Hinxton/SPADES_assemblies/illumina_only/contigs.fa"
# PacBio_assembly="/home/users/tskalicky/Hinxton/PacBio_assembly/HGAP_25_polished_assembly.fa"
# PacBio_assembly="/home/users/tskalicky/Hinxton/PacBio_assembly/HGAP_25_filtered_longreads.fasta"
#
# Adding SPAdes binaries into the PATH
# Will work ONLY on NFS4 connected servers
export PATH="/home/tomas/miniconda2/bin:$PATH"
export PATH="/home/tomas/bioinfo_software:$PATH"
export PATH="/home/tomas/bioinfo_software/companion:$PATH"
export PATH="/home/tomas/bioinfo_software/companion/bin:$PATH"
# export PATH="/storage/brno3-cerit/home/tskalicky/tools/STAR-2.6.0a/bin/Linux_x86_64_static:$PATH" # this version is broken
# Binaries
nextflow=$(which nextflow)
# Check the tools versions
which $nextflow
#
MY_RAM=300 # Max RAM memory for Spades and Samtools sort
#THREADS=$PBS_NUM_PPN
THREADS=30
#
################################################################################
## Commands ##
# mkdir $OUTPUT_DIR
cd $OUTPUT_DIR
# date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# echo "Copy short reads"
# cp $INPUT $OUTPUT_DIR
# date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# echo "Starting canu assembly of Wallacemonas corrected reads"
# #
nextflow run sanger-pathogens/companion -c 06v01_companion_params_Wallacemonas.config -profile docker
