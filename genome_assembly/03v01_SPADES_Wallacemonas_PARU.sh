#!/bin/bash
#PBS -l nodes=1:ppn=40
#PBS -l mem=300gb
#PBS -k oe
#PBS -N 03v01_assembly_wallacemonas_spades_paru
#PBS -M tomas.skalicky@seznam.cz
#PBS -m abe
#
#PARU KRTECEK server is using TORQUE scheduling system !!!
###########################################################################################
# This script will be running on Trypka server 
###########################################################################################
#
## initialize the required application
#
# Requires Spades, samtools, multiQC
#
############################################################################################
### Variables
INPUT_DIR="/media/4TB3/Tomas_data/Hinxton/Trimming/"
WALL_MISEQ1_PE1="/media/4TB3/Tomas_data/Hinxton/Trimming/18098_1_6_1_trimmed_Paired.fq.gz"
WALL_MISEQ1_PE2="/media/4TB3/Tomas_data/Hinxton/Trimming/18098_1_6_2_Paired.fq.gz"
WALL_MISEQ1_UP1="/media/4TB3/Tomas_data/Hinxton/Trimming/18098_1_6_1_trimmed_Unpaired.fq.gz"
WALL_MISEQ1_UP2="/media/4TB3/Tomas_data/Hinxton/Trimming/18098_1_6_2_Unpaired.fq.gz"
WALL_MISEQ1_UP3="/media/4TB3/Tomas_data/Hinxton/Trimming/18098_1_6_both_Unpaired.fq.gz"
#
WALL_MISEQ2_PE1="/media/4TB3/Tomas_data/Hinxton/Trimming/18021_1_6_1_trimmed_Paired.fq.gz"
WALL_MISEQ2_PE2="/media/4TB3/Tomas_data/Hinxton/Trimming/18021_1_6_2_Paired.fq.gz"
WALL_MISEQ2_UP1="/media/4TB3/Tomas_data/Hinxton/Trimming/18021_1_6_1_trimmed_Unpaired.fq.gz"
WALL_MISEQ2_UP2="/media/4TB3/Tomas_data/Hinxton/Trimming/18021_1_6_2_Unpaired.fq.gz"
WALL_MISEQ2_UP3="/media/4TB3/Tomas_data/Hinxton/Trimming/18021_1_6_both_Unpaired.fq"
#
WALL_HISEQ1_1="/media/4TB3/Tomas_data/Hinxton/Trimming/Mbr04_Wallacemonas_both_6_1_trimmed_Paired.fq.gz"
WALL_HISEQ1_2="/media/4TB3/Tomas_data/Hinxton/Trimming/Mbr04_Wallacemonas_both_6_2_trimmed_Paired.fq.gz"
WALL_HISEQ1_UP1="/media/4TB3/Tomas_data/Hinxton/Trimming/Mbr04_Wallacemonas_both_6_1_trimmed_Unpaired.fq.gz"
WALL_HISEQ1_UP2="/media/4TB3/Tomas_data/Hinxton/Trimming/Mbr04_Wallacemonas_both_6_2_trimmed_Unpaired.fq.gz"
WALL_HISEQ1_UP3="/media/4TB3/Tomas_data/Hinxton/Trimming/Mbr04_Wallacemonas_both_6_both_trimmed_Unpaired.fq.gz"
#
WALL_PACBIO="/media/4TB3/Tomas_data/Hinxton/Assembly/PacBio/HGAP/Mbr04_Wallacemonas/Mbr04_filtered_longreads.fasta"
TRUSTED_CONTIGS="/media/4TB3/Tomas_data/Hinxton/Assembly/PacBio/HGAP/Mbr04_Wallacemonas/HGAP_25/polished_assembly.fasta"
#
OUTPUT_DIR="/home/tomas/Hinxton_SPADES_assemblies/no_trusted_contigs"
#
MY_RAM=120 # Max RAM memory for Samtools sort
#THREADS=$PBS_NUM_PPN
THREADS=29
# Adding SPAdes binaries into the PATH
# Will work ONLY on NFS4 connected servers
export PATH="/home/tomas/bioinfo_software/SPAdes-3.12.0-Linux/bin:$PATH"
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
echo "Starting SPAdes assembly of Wallacemonas:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
#
$SPADES -t $THREADS -m $MY_RAM  \
--pe1-1 $WALL_HISEQ1_1 --pe1-2 $WALL_HISEQ1_2 --pe1-s $WALL_HISEQ1_UP3 \
--pe2-1 $WALL_MISEQ1_PE1 --pe2-2 $WALL_MISEQ1_PE2 --pe2-s $WALL_MISEQ1_UP3 \
--pe3-1 $WALL_MISEQ2_PE1 --pe3-2 $WALL_MISEQ2_PE2 --pe3-s $WALL_MISEQ2_UP3 \
--pacbio $WALL_PACBIO -o $OUTPUT_DIR
# $SPADES -t $THREADS -m $MY_RAM  \
# --pe1-1 $WALL_HISEQ1_1 --pe1-2 $WALL_HISEQ1_2 --pe1-s $WALL_HISEQ1_UP3 \
# --pe2-1 $WALL_MISEQ1_PE1 --pe2-2 $WALL_MISEQ1_PE2 --pe2-s $WALL_MISEQ1_UP3 \
# --pe3-1 $WALL_MISEQ2_PE1 --pe3-2 $WALL_MISEQ2_PE2 --pe3-s $WALL_MISEQ2_UP3 \
# --pacbio $WALL_PACBIO --trusted-contigs $TRUSTED_CONTIGS -o $OUTPUT_DIR
#
echo "Assembly finnished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"