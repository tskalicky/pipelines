#!/bin/bash
#PBS -l walltime=24:0:0 
#PBS -q uv@wagap-pro.cerit-sc.cz
#PBS -l select=1:ncpus=1:mem=100gb:scratch_local=100gb
#PBS -j oe
#PBS -N 13v01_awk_renaming_AHB8_annotated
#
## initialize the required application
# module add gcc-4.9.2 # required libraries for bedtools on Zewura
module add bedtools-2.25.0 # version 2.26.0 is BROKEN!
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
CLIP_SES_log2="/storage/brno3-cerit/home/tskalicky/ABH8/anotace/bamcompare_annot/ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_log2_RPKM_smooth/temp_2_tRNA.ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_log2_RPKM_smooth.ren.bed"
CLIP_RPKM_log2="/storage/brno3-cerit/home/tskalicky/ABH8/anotace/bamcompare_annot/ABH8_CLIP_vs_RNAseq_bamCompare_bin10_log2_RPKM_smooth/temp_2_tRNA.ABH8_CLIP_vs_RNAseq_bamCompare_bin10_log2_RPKM_smooth.ren.bed"
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/ABH8/anotace/bamcompare_annot"
####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
# cp -av $ABH8/*.bam.gz $ABH8_POOLED/*.bam.gz $ABH8_PEAKS/*.bed $SCRATCHDIR
cp -av $CLIP_SES_log2 $CLIP_RPKM_log2 $SCRATCHDIR
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR
echo "Following files/folders were copied to scratch:"
ls -c1

### Commands
####################################################################################################
for a in *.bed
do
    FILE=$a
    FILENAME=${a%.*}
    # date +"%d/%m/%Y %H:%M:%S"
    # echo "Renaming reads in file $FILE"
    # awk -F'[\t]' -v OFS='\t' '{print $1, $2, $3, $8, $10, $4}' $FILE > "$SCRATCHDIR/results"/"$FILENAME".ren.bed
    # date +"%d/%m/%Y %H:%M:%S"
    # echo "Finnished renaming reads in file $FILE"
    #
    date +"%d/%m/%Y %H:%M:%S"
    echo "Sorting and merging reads in file $FILE"
    sortBed -i $FILE > "$FILENAME".sorted.bed
    mergeBed -s -c 4,5,6 -o distinct,distinct,distinct -i "$FILENAME".sorted.bed > "$FILENAME".sorted.merged.bed
    # sortBed -i "$SCRATCHDIR/results"/"$FILENAME".ren.bed | \
    # mergeBed -s -c 4,5,6 -o distinct,distinct,distinct -i stdin > "$SCRATCHDIR/results"/"$FILENAME".ren.merged.sorted.bed
    date +"%d/%m/%Y %H:%M:%S"
    echo "Finnished sorting and merging reads in file $FILE"
done &
wait
#
############################################################################################
### Copy data from scratch back to home dir and clean scratch
mkdir -p $OUTPUT_DIR
cp -avr $SCRATCHDIR $OUTPUT_DIR || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
