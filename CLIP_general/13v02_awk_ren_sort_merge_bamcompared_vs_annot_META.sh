#!/bin/bash
#PBS -l walltime=4:0:0 
#PBS -q default@wagap-pro.cerit-sc.cz
#PBS -l select=1:ncpus=2:mem=60gb:scratch_local=100gb:os=debian9
#PBS -j oe
#PBS -N 13v02_AHB8_awk_ren_sort_merge_annotated
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
INPUT_DIR="/storage/brno3-cerit/home/tskalicky/ABH8/anotace/bamcompare_annot/renamed/re-intersect"
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/ABH8/anotace/bamcompare_annot/renamed/re-intersect"
ANNOTATION_DB="/mnt/storage-brno3-cerit/nfs4/home/tskalicky/genomes/human/annotation/my_annot_DB/ensembl/bed/final"
####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
# mkdir $SCRATCHDIR/annotation_db
# cp -av $ANNOTATION_DB/*.bed $SCRATCHDIR/annotation_db
cp -av $INPUT_DIR/*.bed $SCRATCHDIR
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR
echo "Following files/folders were copied to scratch:"
ls -c1

### Commands
####################################################################################################
shopt nullglob # If set, Bash allows filename patterns which match no files to expand to a null string, rather than themselves.
# declare an array variables
declare -a SAMPLES=(*.bed)
declare -a ANNOT_DB=("GRCh38.91.rRNA.bed" "GRCh38.91-tRNAs.ren.bed" "GRCh38.91.snoRNA.bed" \
    "GRCh38.91.snRNA.bed" "GRCh38.91.miRNA.bed" "GRCh38.91.mRNA_exons.merged.bed" \
    "GRCh38.91.mRNA_introns.bed" "GRCh38.91.lincRNA.bed" "ucsc_hg38_repetitive_elements.ren.bed" \
    "GRCh38.91.antisense_mRNA.bed" "GRCh38.91.intergenic.bed")

# get length of an array
smp_num="${#SAMPLES[@]}"
annot_num="${#ANNOT_DB[@]}"
#
echo "There are $smp_num samples that will be processed."
echo "Sample names are: ${SAMPLES[@]}"
#
# mkdir $SCRATCHDIR/results
# mkdir $SCRATCHDIR/results/NOT
# # for a in *.bed
# # do
# #     FILE=$a
# #     FILENAME=${a%.*}
# #     date +"%d/%m/%Y %H:%M:%S"
# #     echo "Renaming reads in file $FILE"
# #     awk -F'[\t]' -v OFS='\t' '{print $1, $2, $3, $8, $4, $10}' $FILE > "$SCRATCHDIR/results"/"$FILENAME".ren.bed
# #     date +"%d/%m/%Y %H:%M:%S"
# #     echo "Finnished renaming reads in file $FILE"
# # done
# ##
# # use for loop to read all values and indexes
# ## rRNA ##
# for a in temp_1_rRNA.*.ren.bed
# do
#     FILE=$a
#     FILENAME=${a%.*.*}
#     date +"%d/%m/%Y %H:%M:%S" 
#     echo "Start intersecting of sample $FILE with database ${ANNOT_DB[0]}"
#     intersectBed -wo -s -f 0.50 -split \
#     -a "$SCRATCHDIR/$FILE" \
#     -b "$SCRATCHDIR/annotation_db/${ANNOT_DB[0]}" \
#     > "$SCRATCHDIR/results/re_$FILENAME.bed"
#     #
#     intersectBed -v -s -split \
#     -a "$SCRATCHDIR/$FILE" \
#     -b "$SCRATCHDIR/annotation_db/${ANNOT_DB[0]}" \
#     > "$SCRATCHDIR/results/NOT/re_NOT_$FILENAME.bed"
#     date +"%d/%m/%Y %H:%M:%S" 
#     echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[0]}"
# done &
## snoRNA ##
# for a in temp_3_snoRNA.*.ren.bed
# do
#     FILE=$a
#     FILENAME=${a%.*.*}
#     date +"%d/%m/%Y %H:%M:%S" 
#     echo "Start intersecting of sample $FILE with database ${ANNOT_DB[2]}"
#     intersectBed -wo -s -f 0.50 -split \
#     -a "$SCRATCHDIR/$FILE" \
#     -b "$SCRATCHDIR/annotation_db/${ANNOT_DB[2]}" \
#     > "$SCRATCHDIR/results/re_$FILENAME.bed"
#     #
#     intersectBed -v -s -split \
#     -a "$SCRATCHDIR/$FILE" \
#     -b "$SCRATCHDIR/annotation_db/${ANNOT_DB[2]}" \
#     > "$SCRATCHDIR/results/NOT/re_NOT_$FILENAME.bed"
#     date +"%d/%m/%Y %H:%M:%S" 
#     echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[2]}"
# done
# ## snRNA ##
# for a in temp_4_snRNA.*.ren.bed
# do
#     FILE=$a
#     FILENAME=${a%.*.*}
#     date +"%d/%m/%Y %H:%M:%S" 
#     echo "Start intersecting of sample $FILE with database ${ANNOT_DB[3]}"
#     intersectBed -wo -s -f 0.50 -split \
#     -a "$SCRATCHDIR/$FILE" \
#     -b "$SCRATCHDIR/annotation_db/${ANNOT_DB[3]}" \
#     > "$SCRATCHDIR/results/re_$FILENAME.bed"
#     #
#     intersectBed -v -s -split \
#     -a "$SCRATCHDIR/$FILE" \
#     -b "$SCRATCHDIR/annotation_db/${ANNOT_DB[3]}" \
#     > "$SCRATCHDIR/results/NOT/re_NOT_$FILENAME.bed"
#     date +"%d/%m/%Y %H:%M:%S" 
#     echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[3]}"
# done
# ## miRNA ##
# for a in temp_5_miRNA.*.ren.bed
# do
#     FILE=$a
#     FILENAME=${a%.*.*}
#     date +"%d/%m/%Y %H:%M:%S" 
#     echo "Start intersecting of sample $FILE with database ${ANNOT_DB[4]}"
#     intersectBed -wo -s -f 0.50 -split \
#     -a "$SCRATCHDIR/$FILE" \
#     -b "$SCRATCHDIR/annotation_db/${ANNOT_DB[4]}" \
#     > "$SCRATCHDIR/results/re_$FILENAME.bed"
#     #
#     intersectBed -v -s -split \
#     -a "$SCRATCHDIR/$FILE" \
#     -b "$SCRATCHDIR/annotation_db/${ANNOT_DB[4]}" \
#     > "$SCRATCHDIR/results/NOT/re_NOT_$FILENAME.bed"
#     date +"%d/%m/%Y %H:%M:%S" 
#     echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[4]}"
# done &
# ## mRNA exons##
# for a in temp_6_mRNA_exons.*.ren.bed
# do
#     FILE=$a
#     FILENAME=${a%.*.*}
#     date +"%d/%m/%Y %H:%M:%S" 
#     echo "Start intersecting of sample $FILE with database ${ANNOT_DB[5]}"
#     intersectBed -wo -s -f 0.50 -split \
#     -a "$SCRATCHDIR/$FILE" \
#     -b "$SCRATCHDIR/annotation_db/${ANNOT_DB[5]}" \
#     > "$SCRATCHDIR/results/re_$FILENAME.bed"
#     #
#     intersectBed -v -s -split \
#     -a "$SCRATCHDIR/$FILE" \
#     -b "$SCRATCHDIR/annotation_db/${ANNOT_DB[5]}" \
#     > "$SCRATCHDIR/results/NOT/re_NOT_$FILENAME.bed"
#     date +"%d/%m/%Y %H:%M:%S" 
#     echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[5]}"
# done
# ## mRNA introns ##
# for a in temp_7_mRNA_introns.*.ren.bed
# do
#     FILE=$a
#     FILENAME=${a%.*.*}
#     date +"%d/%m/%Y %H:%M:%S" 
#     echo "Start intersecting of sample $FILE with database ${ANNOT_DB[6]}"
#     intersectBed -wo -s -f 0.50 -split \
#     -a "$SCRATCHDIR/$FILE" \
#     -b "$SCRATCHDIR/annotation_db/${ANNOT_DB[6]}" \
#     > "$SCRATCHDIR/results/re_$FILENAME.bed"
#     #
#     intersectBed -v -s -split \
#     -a "$SCRATCHDIR/$FILE" \
#     -b "$SCRATCHDIR/annotation_db/${ANNOT_DB[6]}" \
#     > "$SCRATCHDIR/results/NOT/re_NOT_$FILENAME.bed"
#     date +"%d/%m/%Y %H:%M:%S" 
#     echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[6]}"
# done &
# ## linc RNA introns ##
# for a in temp_8_lincRNA.*.ren.bed
# do
#     FILE=$a
#     FILENAME=${a%.*.*}
#     date +"%d/%m/%Y %H:%M:%S" 
#     echo "Start intersecting of sample $FILE with database ${ANNOT_DB[7]}"
#     intersectBed -wo -s -f 0.50 -split \
#     -a "$SCRATCHDIR/$FILE" \
#     -b "$SCRATCHDIR/annotation_db/${ANNOT_DB[7]}" \
#     > "$SCRATCHDIR/results/re_$FILENAME.bed"
#     #
#     intersectBed -v -s -split \
#     -a "$SCRATCHDIR/$FILE" \
#     -b "$SCRATCHDIR/annotation_db/${ANNOT_DB[7]}" \
#     > "$SCRATCHDIR/results/NOT/re_NOT_$FILENAME.bed"
#     date +"%d/%m/%Y %H:%M:%S" 
#     echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[7]}"
# done
# ## repetitive elements ##
# for a in temp_9_repetitive_elements.*.ren.bed
# do
#     FILE=$a
#     FILENAME=${a%.*.*}
#     date +"%d/%m/%Y %H:%M:%S" 
#     echo "Start intersecting of sample $FILE with database ${ANNOT_DB[8]}"
#     intersectBed -wo -s -f 0.50 -split \
#     -a "$SCRATCHDIR/$FILE" \
#     -b "$SCRATCHDIR/annotation_db/${ANNOT_DB[8]}" \
#     > "$SCRATCHDIR/results/re_$FILENAME.bed"
#     #
#     intersectBed -v -s -split \
#     -a "$SCRATCHDIR/$FILE" \
#     -b "$SCRATCHDIR/annotation_db/${ANNOT_DB[8]}" \
#     > "$SCRATCHDIR/results/NOT/re_NOT_$FILENAME.bed"
#     date +"%d/%m/%Y %H:%M:%S" 
#     echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[8]}"
# done &
# ## antisense RNA ##
# for a in temp_10_antisense_RNA.*.ren.bed
# do
#     FILE=$a
#     FILENAME=${a%.*.*}
#     date +"%d/%m/%Y %H:%M:%S" 
#     echo "Start intersecting of sample $FILE with database ${ANNOT_DB[9]}"
#     intersectBed -wo -s -f 0.50 -split \
#     -a "$SCRATCHDIR/$FILE" \
#     -b "$SCRATCHDIR/annotation_db/${ANNOT_DB[9]}" \
#     > "$SCRATCHDIR/results/re_$FILENAME.bed"
#     #
#     intersectBed -v -s -split \
#     -a "$SCRATCHDIR/$FILE" \
#     -b "$SCRATCHDIR/annotation_db/${ANNOT_DB[9]}" \
#     > "$SCRATCHDIR/results/NOT/re_NOT_$FILENAME.bed"
#     date +"%d/%m/%Y %H:%M:%S" 
#     echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[9]}"
# done
# ## intergenic RNA ##
# for a in temp_11_intergenic.*.ren.bed
# do
#     FILE=$a
#     FILENAME=${a%.*.*}
#     date +"%d/%m/%Y %H:%M:%S" 
#     echo "Start intersecting of sample $FILE with database ${ANNOT_DB[10]}"
#     intersectBed -wo -s -f 0.50 -split \
#     -a "$SCRATCHDIR/$FILE" \
#     -b "$SCRATCHDIR/annotation_db/${ANNOT_DB[10]}" \
#     > "$SCRATCHDIR/results/re_$FILENAME.bed"
#     #
#     intersectBed -v -s -split \
#     -a "$SCRATCHDIR/$FILE" \
#     -b "$SCRATCHDIR/annotation_db/${ANNOT_DB[10]}" \
#     > "$SCRATCHDIR/results/re_NOT_$FILENAME.bed"
#     date +"%d/%m/%Y %H:%M:%S" 
#     echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[10]}"
# done
# wait
#
# cd $SCRATCHDIR/results
for b in re_*.bed
do
    FILE2=$b
    FILENAME2=${b%.*}
    date +"%d/%m/%Y %H:%M:%S"
    echo "Sorting and merging reads in file $FILE2"
    awk -F'[\t]' -v OFS='\t' '{print $1, $2, $3, $10, $5, $12}' "$FILE2" > "$FILENAME2".ren2.bed
    sortBed -i "$FILENAME2".ren2.bed > "$FILENAME2".ren2.sorted.bed
    mergeBed -s -c 4,5,6 -o distinct,median,distinct -i "$FILENAME2".ren2.sorted.bed > "$FILENAME2".ren2.sorted.merged.bed
    # Need to delete extra field containing strain +,- info on bad position
    awk -F'[\t]' -v OFS='\t' '{print $1, $2, $3, $5, $6, $7}' "$FILENAME2".ren2.sorted.merged.bed > "$FILENAME2".ren3.sorted.merged.bed
    date +"%d/%m/%Y %H:%M:%S"
    echo "Finnished sorting and merging reads in file $FILE2"
done
wait
#
############################################################################################
### Copy data from scratch back to home dir and clean scratch
# mkdir -p $OUTPUT_DIR
cp -avr $SCRATCHDIR $OUTPUT_DIR || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
