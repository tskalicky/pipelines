#!/bin/bash
#PBS -l nodes=1:ppn=30
#PBS -l mem=300gb
#PBS -l walltime=168:00:00
#PBS -k oe
#PBS -N 03v04_assembly_wallacemonas_spades_carefull_krtecek
#PBS -M tomas.skalicky@seznam.cz
#PBS -m abe
#
#PARU KRTECEK server is using TORQUE scheduling system !!!
###########################################################################################
# This script will be running on Krtecek server 
###########################################################################################
#
## initialize the required application
#
# Requires Spades
# Assembly on pre-corrected reads from HALC software
#
############################################################################################
### Variables
INPUT="/home/users/tskalicky/Hinxton/error_corrections/output/Wallacemonas_HALC.trim.fa"
# WALL_MISEQ1_PE1="/home/users/tskalicky/Hinxton/trimming/18098_1_6_1_trimmed_Paired.fq.gz"
# WALL_MISEQ1_PE2="/home/users/tskalicky/Hinxton/trimming/18098_1_6_2_Paired.fq.gz"
# WALL_MISEQ1_UP1="/home/users/tskalicky/Hinxton/trimming/18098_1_6_1_trimmed_Unpaired.fq.gz"
# WALL_MISEQ1_UP2="/home/users/tskalicky/Hinxton/trimming/18098_1_6_2_Unpaired.fq.gz"
# WALL_MISEQ1_UP3="/home/users/tskalicky/Hinxton/trimming/18098_1_6_both_Unpaired.fq.gz"
# #
# WALL_MISEQ2_PE1="/home/users/tskalicky/Hinxton/trimming/18021_1_6_1_trimmed_Paired.fq.gz"
# WALL_MISEQ2_PE2="/home/users/tskalicky/Hinxton/trimming/18021_1_6_2_Paired.fq.gz"
# WALL_MISEQ2_UP1="/home/users/tskalicky/Hinxton/trimming/18021_1_6_1_trimmed_Unpaired.fq.gz"
# WALL_MISEQ2_UP2="/home/users/tskalicky/Hinxton/trimming/18021_1_6_2_Unpaired.fq.gz"
# WALL_MISEQ2_UP3="/home/users/tskalicky/Hinxton/trimming/18021_1_6_both_Unpaired.fq"
# #
# WALL_HISEQ1_1="/home/users/tskalicky/Hinxton/trimming/Mbr04_Wallacemonas_both_6_1_trimmed_Paired.fq.gz"
# WALL_HISEQ1_2="/home/users/tskalicky/Hinxton/trimming/Mbr04_Wallacemonas_both_6_2_trimmed_Paired.fq.gz"
# WALL_HISEQ1_UP1="/home/users/tskalicky/Hinxton/trimming/Mbr04_Wallacemonas_both_6_1_trimmed_Unpaired.fq.gz"
# WALL_HISEQ1_UP2="/home/users/tskalicky/Hinxton/trimming/Mbr04_Wallacemonas_both_6_2_trimmed_Unpaired.fq.gz"
# WALL_HISEQ1_UP3="/home/users/tskalicky/Hinxton/trimming/Mbr04_Wallacemonas_both_6_both_trimmed_Unpaired.fq.gz"
# #
WALL_PACBIO="/home/users/tskalicky/Hinxton/PacBio_assembly/HGAP_25_filtered_longreads.fasta"
# TRUSTED_CONTIGS="/home/users/tskalicky/Hinxton/PacBio_assembly/HGAP_25/polished_assembly.fasta"
#
OUTPUT_DIR="/home/users/tskalicky/Hinxton/SPADES_assemblies/HALC_precorrected/trusted"
#
MY_RAM=300 # Max RAM memory for Spades and Samtools sort
#THREADS=$PBS_NUM_PPN
THREADS=30
# Adding SPAdes binaries into the PATH
# Will work ONLY on NFS4 connected servers
export PATH="/home/users/tskalicky/anaconda2/bin:$PATH"
# export PATH="/storage/brno3-cerit/home/tskalicky/tools/STAR-2.6.0a/bin/Linux_x86_64_static:$PATH" # this version is broken
# Binaries
SPADES=$(which spades.py)
# Check the tools versions
which $SPADES
#
################################################################################
## Commands ##
mkdir $OUTPUT_DIR
cd $OUTPUT_DIR
echo "Starting SPAdes assembly of Wallacemonas using pre-corrected reads from HALC:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
#
$SPADES -t $THREADS -m $MY_RAM --only-assembler --careful \
--s1 $INPUT \
--trusted-contigs $WALL_PACBIO \
-o $OUTPUT_DIR
# $SPADES -t $THREADS -m $MY_RAM  \
# --pe1-1 $WALL_HISEQ1_1 --pe1-2 $WALL_HISEQ1_2 --pe1-s $WALL_HISEQ1_UP3 \
# --pe2-1 $WALL_MISEQ1_PE1 --pe2-2 $WALL_MISEQ1_PE2 --pe2-s $WALL_MISEQ1_UP3 \
# --pe3-1 $WALL_MISEQ2_PE1 --pe3-2 $WALL_MISEQ2_PE2 --pe3-s $WALL_MISEQ2_UP3 \
# --pacbio $WALL_PACBIO --trusted-contigs $TRUSTED_CONTIGS -o $OUTPUT_DIR
#
echo "Assembly finnished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"