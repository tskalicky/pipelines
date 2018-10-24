#!/bin/bash
#PBS -l nodes=1:ppn=30
#PBS -l mem=300gb
#PBS -l walltime=168:00:00
#PBS -k oe
#PBS -N 04v01_error_correction_wallacemonas_HALC_paru
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
# Requires HALC, blasr, lordec
#
# NOTE for HALC:
# Mandatory Inputs are:
# 
# 1) Long reads in FASTA format.
# 2) mContigs assembled from the corresponding short reads in FASTA format.
# 3) The initial short reads in FASTA format (only for -ordinary mode; 
# obtained with cat left_reads.fa >short_reads.fa and then cat right_reads.fa >>short_reads.fa).
#
############################################################################################
### Variables
INPUT_DIR="/home/users/tskalicky/Hinxton/trimming/"
OUTPUT_DIR="/home/users/tskalicky/Hinxton/error_corrections"
CONTIGS="/home/users/tskalicky/Hinxton/SPADES_assemblies/illumina_only/contigs.fa"
# PacBio_assembly="/home/users/tskalicky/Hinxton/PacBio_assembly/HGAP_25_polished_assembly.fa"
PacBio_assembly="/home/users/tskalicky/Hinxton/PacBio_assembly/HGAP_25_filtered_longreads.fasta"
#
WALL_MISEQ1_PE1="/home/users/tskalicky/Hinxton/trimming/18098_1_6_1_trimmed_Paired.fq.gz"
WALL_MISEQ1_PE2="/home/users/tskalicky/Hinxton/trimming/18098_1_6_2_Paired.fq.gz"
WALL_MISEQ1_UP1="/home/users/tskalicky/Hinxton/trimming/18098_1_6_1_trimmed_Unpaired.fq.gz"
WALL_MISEQ1_UP2="/home/users/tskalicky/Hinxton/trimming/18098_1_6_2_Unpaired.fq.gz"
WALL_MISEQ1_UP3="/home/users/tskalicky/Hinxton/trimming/18098_1_6_both_Unpaired.fq.gz"
#
WALL_MISEQ2_PE1="/home/users/tskalicky/Hinxton/trimming/18021_1_6_1_trimmed_Paired.fq.gz"
WALL_MISEQ2_PE2="/home/users/tskalicky/Hinxton/trimming/18021_1_6_2_Paired.fq.gz"
WALL_MISEQ2_UP1="/home/users/tskalicky/Hinxton/trimming/18021_1_6_1_trimmed_Unpaired.fq.gz"
WALL_MISEQ2_UP2="/home/users/tskalicky/Hinxton/trimming/18021_1_6_2_Unpaired.fq.gz"
WALL_MISEQ2_UP3="/home/users/tskalicky/Hinxton/trimming/18021_1_6_both_Unpaired.fq"
#
WALL_HISEQ1_1="/home/users/tskalicky/Hinxton/trimming/Mbr04_Wallacemonas_both_6_1_trimmed_Paired.fq.gz"
WALL_HISEQ1_2="/home/users/tskalicky/Hinxton/trimming/Mbr04_Wallacemonas_both_6_2_trimmed_Paired.fq.gz"
WALL_HISEQ1_UP1="/home/users/tskalicky/Hinxton/trimming/Mbr04_Wallacemonas_both_6_1_trimmed_Unpaired.fq.gz"
WALL_HISEQ1_UP2="/home/users/tskalicky/Hinxton/trimming/Mbr04_Wallacemonas_both_6_2_trimmed_Unpaired.fq.gz"
WALL_HISEQ1_UP3="/home/users/tskalicky/Hinxton/trimming/Mbr04_Wallacemonas_both_6_both_trimmed_Unpaired.fq.gz"
#
WALL_PACBIO="/home/users/tskalicky/Hinxton/PacBio_assembly/Mbr04_filtered_longreads.fasta"
#
MY_RAM=300 # Max RAM memory for Spades and Samtools sort
#THREADS=$PBS_NUM_PPN
THREADS=30
# Adding SPAdes binaries into the PATH
# Will work ONLY on NFS4 connected servers
export PATH="/usr/bin:$PATH"
export PATH="/home/users/tskalicky/anaconda2/bin:$PATH"
export PATH="/home/users/tskalicky/software/halc/bin:$PATH"
export PATH="/home/users/tskalicky/software/halc:$PATH"
# export PATH="/storage/brno3-cerit/home/tskalicky/tools/STAR-2.6.0a/bin/Linux_x86_64_static:$PATH" # this version is broken
# Binaries
runHALC=$(which runHALC.py)
UNPIGZ="/usr/bin/unpigz"
# Check the tools versions
which $runHALC
#
################################################################################
## Commands ##
mkdir $OUTPUT_DIR
cd $OUTPUT_DIR
# date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# echo "Copy short reads"
# cp -av $WALL_MISEQ1_PE1 $WALL_MISEQ1_PE2 $WALL_MISEQ2_PE1 $WALL_MISEQ2_PE2 $WALL_HISEQ1_1 $WALL_HISEQ1_2 $OUTPUT_DIR
# echo "Following files were copied to scratch:"
# ls -Rc1
#
# date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# echo "Decompressing files:"
# for b in *.fq.gz
# do
# 	gzip -dv $b &
# done
# wait
# date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# echo "Finnished decompressing files"
# #
# date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# echo "Converting files from fastq to fasta format"
# #
# for a in *.fq
# do
# 	FILE=$a
# 	FILENAME=${a%.*}
# 	sed -n '1~4s/^@/>/p;2~4p' $FILE > $FILENAME".fasta"
# done
# wait
# #
# date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# echo "Done converting files from fastq to fasta format"
#
# WALL_MISEQ1_PE1=$(basename $WALL_MISEQ1_PE1)
# WALL_MISEQ1_PE2=$(basename $WALL_MISEQ1_PE2)
# #
# WALL_MISEQ2_PE1=$(basename $WALL_MISEQ2_PE1)
# WALL_MISEQ2_PE2=$(basename $WALL_MISEQ2_PE2)
# #
# WALL_HISEQ1_1=$(basename $WALL_HISEQ1_1)
# WALL_HISEQ1_2=$(basename $WALL_HISEQ1_2)
# #
# #
# WALL_MISEQ1_PE1=${WALL_MISEQ1_PE1%.*.*}
# WALL_MISEQ1_PE2=${WALL_MISEQ1_PE2%.*.*}
# #
# WALL_MISEQ2_PE1=${WALL_MISEQ2_PE1%.*.*}
# WALL_MISEQ2_PE2=${WALL_MISEQ2_PE2%.*.*}
# #
# WALL_HISEQ1_1=${WALL_HISEQ1_1%.*.*}
# WALL_HISEQ1_2=${WALL_HISEQ1_2%.*.*}
# #
# date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# echo "Preparing short reads for HALC error corrections of Wallacemonas PacBio assembly"
# #
# cat $WALL_MISEQ1_PE1".fasta" $WALL_MISEQ2_PE1".fasta" $WALL_HISEQ1_1".fasta" >>short_reads.fa 
# cat $WALL_MISEQ1_PE2".fasta" $WALL_MISEQ2_PE2".fasta" $WALL_HISEQ1_2".fasta" >>short_reads.fa
# #
# date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# echo "Finnished preparing short reads for HALC error corrections of Wallacemonas PacBio assembly"
#
echo "Starting HALC error corrections of Wallacemonas PacBio assembly:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
#
# python runHALC.py long_reads.fa contigs.fa [-options|-options]
python $runHALC $PacBio_assembly $CONTIGS --threads $THREADS -log -o short_reads.fa
#
echo "Finnished HALC error corrections of Wallacemonas PacBio assembly:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"