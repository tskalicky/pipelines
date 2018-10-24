#!/bin/bash
#PBS -l nodes=krtecek2:ppn=5
#PBS -l mem=300gb
#PBS -j oe
#PBS -N 13v01_awk_renaming_AHB8_annotated
#
#PARU KRTECEK server is using TORQUE scheduling system !!!
#
## initialize the required application
# module add gcc-4.9.2 # required libraries for bedtools on Zewura
# module add bedtools-2.25.0 # version 2.26.0 is BROKEN!
# module add bamtools
############################################################################################
# CLIP Seq annotation pipeline
# 1) Changing names of bed intervals to match identified gene names using awk
# 2) Sorting and merging renamed intervals using Bedtools
#
############################################################################################
# Check if tools are installed
BEDTOOLS=$(which bedtools)
SORTBED=$(which sortBed)
MERGEBED=$(which mergeBed)
COMPLEMENTBED=$(which complementBed)
INTERSECTBED=$(which intersectBed)
BAMTOBED=$(which bamToBed)
#
which $BEDTOOLS
which $SORTBED
which $MERGEBED
which $COMPLEMENTBED
which $INTERSECTBED
#
# variables
CLIP_SES_log2="temp_2_tRNA.ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_log2_RPKM_smooth.ren.bed"
CLIP_RPKM_log2="temp_2_tRNA.ABH8_CLIP_vs_RNAseq_bamCompare_bin10_log2_RPKM_smooth.ren.bed"
OUTPUT_DIR="/home/users/tskalicky/CEITEC/ABH8/anotace/final/INTERSECT/renamed"
### Commands
####################################################################################################
cd $OUTPUT_DIR
# for a in *.bed
# do
#     FILE=$a
#     FILENAME=${a%.*}
#     # date +"%d/%m/%Y %H:%M:%S"
#     # echo "Renaming reads in file $FILE"
#     # awk -F'[\t]' -v OFS='\t' '{print $1, $2, $3, $8, $10, $4}' $FILE > "$SCRATCHDIR/results"/"$FILENAME".ren.bed
#     # date +"%d/%m/%Y %H:%M:%S"
#     # echo "Finnished renaming reads in file $FILE"
#     #
#     date +"%d/%m/%Y %H:%M:%S"
#     echo "Sorting and merging reads in file $FILE"
#     sortBed -i $FILE > "$FILENAME".sorted.bed
#     mergeBed -s -c 4,5,6 -o distinct,distinct,distinct -i "$FILENAME".sorted.bed > "$FILENAME".sorted.merged.bed
#     # sortBed -i "$SCRATCHDIR/results"/"$FILENAME".ren.bed | \
#     # mergeBed -s -c 4,5,6 -o distinct,distinct,distinct -i stdin > "$SCRATCHDIR/results"/"$FILENAME".ren.merged.sorted.bed
#     date +"%d/%m/%Y %H:%M:%S"
#     echo "Finnished sorting and merging reads in file $FILE"
# done
mergeBed -s -c 4,5,6 -o distinct,distinct,distinct -i temp_2_tRNA.ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_log2_RPKM_smooth.ren.sorted.bed > tRNA.ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_log2_RPKM_smooth.ren.sorted.merged.bed
mergeBed -s -c 4,5,6 -o distinct,distinct,distinct -i temp_2_tRNA.ABH8_CLIP_vs_RNAseq_bamCompare_bin10_log2_RPKM_smooth.ren.bed > tRNA.ABH8_CLIP_vs_RNAseq_bamCompare_bin10_log2_RPKM_smooth.ren.sorted.merged.bed
wait
#
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
