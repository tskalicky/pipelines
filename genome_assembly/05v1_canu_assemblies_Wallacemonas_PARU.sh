#!/bin/bash
#PBS -l nodes=1:ppn=30
#PBS -l mem=300gb
#PBS -l walltime=168:00:00
#PBS -k oe
#PBS -N 05v01_canu_assembly_wallacemonas_krtecek
#PBS -M tomas.skalicky@seznam.cz
#PBS -m abe
### MUST READ ###
# CANU is detecting and creating jobs for PBS schedulding environment by itself!!! 
# So just run canu with required options and it will compose and submit the scripts for you.
# If you wanna see which settings will canu use add "useGrid=remote" option
#PARU KRTECEK server is using TORQUE scheduling system !!!
###########################################################################################
# This script will be running on Krtecek server 
###########################################################################################
#
## initialize the required application
#
# Requires HALC, blasr, lordec, canu
#
# NOTE for HALC:
# Mandatory Inputs are:
# 
# 1) Long reads in FASTA format.
# 2) mContigs assembled from the corresponding short reads in FASTA format.
# 3) The initial short reads in FASTA format (only for -ordinary mode; 
# obtained with cat left_reads.fa >short_reads.fa and then cat right_reads.fa >>short_reads.fa).
# Output:
# - Error corrected full long reads.
# - Error corrected trimmed long reads.
# - Error corrected split long reads.
#
############################################################################################
### Variables
INPUT="/home/users/tskalicky/Hinxton/error_corrections/output/Wallacemonas_HALC.trim.fa"
OUTPUT_DIR="/home/users/tskalicky/Hinxton/CANU_assembly"
# CONTIGS="/home/users/tskalicky/Hinxton/SPADES_assemblies/illumina_only/contigs.fa"
# PacBio_assembly="/home/users/tskalicky/Hinxton/PacBio_assembly/HGAP_25_polished_assembly.fa"
# PacBio_assembly="/home/users/tskalicky/Hinxton/PacBio_assembly/HGAP_25_filtered_longreads.fasta"
#
# Adding SPAdes binaries into the PATH
# Will work ONLY on NFS4 connected servers
export PATH="/usr/bin:$PATH"
export PATH="/home/users/tskalicky/anaconda2/bin:$PATH"
export PATH="/home/users/tskalicky/software/halc/bin:$PATH"
export PATH="/home/users/tskalicky/software/halc:$PATH"
# export PATH="/storage/brno3-cerit/home/tskalicky/tools/STAR-2.6.0a/bin/Linux_x86_64_static:$PATH" # this version is broken
# Binaries
canu=$(which canu)
# Check the tools versions
which $canu
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
# canu -trim-assemble \
# useGrid=false -maxMemory=$MY_RAM -maxThreads=$THREADS \
# -p Wallacemonas -d wallacemonas-25m-erate-0.045 \
# genomeSize=25m \
# correctedErrorRate=0.045 \
# -pacbio-corrected Wallacemonas_HALC.trim.fa
# date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# echo "Finnished canu assembly of Wallacemonas size 25M erate 0.045"
#
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Starting canu assembly of Wallacemonas size 25M erate 0.075"
canu -trim-assemble \
useGrid=false -maxMemory=$MY_RAM -maxThreads=$THREADS \
-p Wallacemonas -d wallacemonas-25m-erate-0.075 \
genomeSize=25m \
correctedErrorRate=0.075 \
-pacbio-corrected Wallacemonas_HALC.trim.fa
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Finnished canu assembly of Wallacemonas size 25M erate 0.075"
#
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Starting canu assembly of Wallacemonas size 50M erate 0.045"
canu -trim-assemble \
useGrid=false -maxMemory=$MY_RAM -maxThreads=$THREADS \
-p Wallacemonas -d wallacemonas-50m-erate-0.045 \
genomeSize=50m \
correctedErrorRate=0.045 \
-pacbio-corrected Wallacemonas_HALC.trim.fa
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Finnished canu assembly of Wallacemonas size 50M erate 0.045"
#
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Finnished canu assembly of Wallacemonas size 50M erate 0.075"
canu -trim-assemble \
useGrid=false -maxMemory=$MY_RAM -maxThreads=$THREADS \
-p Wallacemonas -d wallacemonas-50m-erate-0.075 \
genomeSize=50m \
correctedErrorRate=0.075 \
-pacbio-corrected Wallacemonas_HALC.trim.fa
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Finnished canu assembly of Wallacemonas size 50M erate 0.075"
#
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Finnished all canu assemblies of Wallacemonas corrected reads"
# #
# #
# canu -trim-assemble useGrid=false -p Wallacemonas -d wallacemonas-25m-erate-0.045 genomeSize=25m correctedErrorRate=0.045 -d /home/users/tskalicky/Hinxton/CANU_assembly -pacbio-corrected Wallacemonas_HALC.trim.fa
# canu -trim-assemble useGrid=false -maxMemory=$MY_RAM -maxThreads=$THREADS -p Wallacemonas -d wallacemonas-25m-erate-0.045 genomeSize=25m correctedErrorRate=0.045 -pacbio-corrected Wallacemonas_HALC.trim.fa
