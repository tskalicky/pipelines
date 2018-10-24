#!/bin/bash
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
ENSEMBL_DB="$HOME/ownCloud/CEITEC_lab/genomes/human/ensembl91"
INPUT="/run/media/skalda/My_Passport_800GB/CEITEC_lab/ABH8/mapping"
INPUT2="$HOME/ownCloud/CEITEC_lab/ABH8/mapping"
OUTPUT="$HOME/ownCloud/CEITEC_lab/ABH8/anotace"
ANNOTATION_DB="$HOME/ownCloud/CEITEC_lab/genomes/annotation/my_annot_DB/bed"
# commands
# if [ $# -ne 2 ]; then
#   echo "Please set two arguments: sample file"
#   exit 0
# fi
#

filename="${file%.*}"
filename="${filename%.*}"
cd $INPUT
# BAM to BED conversion
# for a in *.bam.gz
# do
#   FILE=$a
#   FILENAME=${a%.*.*}
#   gzip -dc $FILE | $BAMTOBED -i stdin > $FILENAME.bed
# done
# wait
#
shopt -s nullglob # If set, Bash allows filename patterns which match no files to expand to a null string, rather than themselves.
# declare an array variables
declare -a SAMPLES=(*.bed)
declare -a ANNOT_DB=("GRCh38.91.rRNA.sorted.bed" "GRCh38.91-tRNAs.ren.sorted.bed" "GRCh38.91.snoRNA.sorted.bed" \
  "GRCh38.91.snRNA.sorted.bed" "GRCh38.91.miRNA.sorted.bed" "GRCh38.91.mRNA_exons.merged.sorted.bed" 
  "GRCh38.91.mRNA_introns.sorted.bed" "GRCh38.91.lincRNA.sorted.bed" "ucsc_hg38_repetitive_elements.ren.sorted.bed" \
  "GRCh38.91.antisense_mRNA.sorted.bed" "GRCh38.91.intergenic.sorted.bed")
# get length of an array
smp_num="${#SAMPLES[@]}"
annot_num="${#ANNOT_DB[@]}"
#
# use for loop to read all values and indexes
for (( a=0;a<=${smp_num};a++ )); # I know, I know programmers COUNT from zero :-D
do
  SAMPLE="${SAMPLES[$a]}"
  SAMPLENAME="${SAMPLE%.*.*}"
  mkdir $OUTPUT/$SAMPLENAME
  #
  date +"%d/%m/%Y %H:%M:%S"
  echo "Start intersecting of sample $SAMPLE with databases."
  ## rRNA ##
  $INTERSECTBED -wo -s \
  -a "$INPUT/$SAMPLE" \
  -b "$ANNOTATION_DB/GRCh38.91.rRNA.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_1_rRNA.$SAMPLENAME.bed"
  #
  $INTERSECTBED -v -s \
  -a "$INPUT/$SAMPLE" \
  -b "$ANNOTATION_DB/GRCh38.91.rRNA.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_1_NOT.$SAMPLENAME.bed"
  date +"%d/%m/%Y %H:%M:%S" 
  echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[0]}"
  ## tRNA ##
  intersectBed -wo -s \
  -a "$OUTPUT/$SAMPLENAME/temp_1_NOT.$SAMPLENAME.bed" \
  -b "$ANNOTATION_DB/GRCh38.91-tRNAs.ren.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_2_tRNA.$SAMPLENAME.bed"
  #
  intersectBed -v -s \
  -a "$OUTPUT/$SAMPLENAME/temp_1_NOT.$SAMPLENAME.bed" \
  -b "$ANNOTATION_DB/GRCh38.91-tRNAs.ren.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_2_NOT.$SAMPLENAME.bed"
  echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[1]}"
  ## snoRNA
  intersectBed -wo -s \
  -a "$OUTPUT/$SAMPLENAME/temp_2_NOT.$SAMPLENAME.bed" \
  -b "$ANNOTATION_DB/GRCh38.91.snoRNA.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_3_snoRNA.$SAMPLENAME.bed"
  #
  intersectBed -v -s \
  -a "$OUTPUT/$SAMPLENAME/temp_2_NOT.$SAMPLENAME.bed" \
  -b "$ANNOTATION_DB/GRCh38.91.snoRNA.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_3_NOT.$SAMPLENAME.bed"
  echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[2]}"
  ## snRNA
  intersectBed -wo -s \
  -a "$OUTPUT/$SAMPLENAME/temp_3_NOT.$SAMPLENAME.bed" \
  -b "$ANNOTATION_DB/GRCh38.91.snRNA.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_4_snRNA.$SAMPLENAME.bed"
  #
  intersectBed -v -s \
  -a "$OUTPUT/$SAMPLENAME/temp_3_NOT.$SAMPLENAME.bed" \
  -b "$ANNOTATION_DB/GRCh38.91.snRNA.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_4_NOT.$SAMPLENAME.bed"
  echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[3]}"
  ## miRNA
  intersectBed -wo -s \
  -a "$OUTPUT/$SAMPLENAME/temp_4_NOT.$SAMPLENAME.bed" \
  -b "$ANNOTATION_DB/GRCh38.91.miRNA.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_5_miRNA.$SAMPLENAME.bed"
  #
  intersectBed -v -s \
  -a "$OUTPUT/$SAMPLENAME/temp_4_NOT.$SAMPLENAME.bed" \
  -b "$ANNOTATION_DB/GRCh38.91.miRNA.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_5_NOT.$SAMPLENAME.bed"
  echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[4]}"
  ## exons mRNA
  intersectBed -wo -s \
  -a "$OUTPUT/$SAMPLENAME/temp_5_NOT.$SAMPLENAME.bed" \
  -b "$ANNOTATION_DB/GRCh38.91.mRNA_exons.merged.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_6_mRNA_exons.$SAMPLENAME.bed"
  #
  intersectBed -v -s \
  -a "$OUTPUT/$SAMPLENAME/temp_5_NOT.$SAMPLENAME.bed" \
  -b "$ANNOTATION_DB/GRCh38.91.mRNA_exons.merged.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_6_NOT.$SAMPLENAME.bed"
  echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[5]}"
  ## introns mRNA
  intersectBed -wo -s \
  -a "$OUTPUT/$SAMPLENAME/temp_6_NOT.$SAMPLENAME.bed" \
  -b "$ANNOTATION_DB/GRCh38.91.mRNA_introns.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_7_mRNA_introns.$SAMPLENAME.bed"
  #
  intersectBed -v -s \
  -a "$OUTPUT/$SAMPLENAME/temp_6_NOT.$SAMPLENAME.bed" \
  -b "$ANNOTATION_DB/GRCh38.91.mRNA_introns.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_7_NOT.$SAMPLENAME.bed"
  echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[6]}"
  ## lincRNA
  intersectBed -wo -s \
  -a "$OUTPUT/$SAMPLENAME/temp_7_NOT.$SAMPLENAME.bed" \
  -b "$ANNOTATION_DB/GRCh38.91.lincRNA.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_8_lincRNA.$SAMPLENAME.bed"
  #
  intersectBed -v -s \
  -a "$OUTPUT/$SAMPLENAME/temp_7_NOT.$SAMPLENAME.bed" \
  -b "$ANNOTATION_DB/GRCh38.91.lincRNA.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_8_NOT.$SAMPLENAME.bed"
  echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[7]}"
  ## repetitive elements
  intersectBed -wo -s \
  -a "$OUTPUT/$SAMPLENAME/temp_8_NOT.$SAMPLENAME.bed" \
  -b "$ANNOTATION_DB/ucsc_hg38_repetitive_elements.ren.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_9_repetitive_elements.$SAMPLENAME.bed"
  #
  intersectBed -v -s \
  -a "$OUTPUT/$SAMPLENAME/temp_8_NOT.$SAMPLENAME.bed" \
  -b "$ANNOTATION_DB/ucsc_hg38_repetitive_elements.ren.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_9_NOT.$SAMPLENAME.bed"
  echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[8]}"
  ## antisense RNA
  intersectBed -wo -s \
  -a "$OUTPUT/$SAMPLENAME/temp_9_NOT.$SAMPLENAME.bed" \
  -b "$ANNOTATION_DB/GRCh38.91.antisense_mRNA.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_10_antisense_RNA.$SAMPLENAME.bed"
  #
  intersectBed -v -s \
  -a "$OUTPUT/$SAMPLENAME/temp_9_NOT.$SAMPLENAME.bed" \
  -b "$ANNOTATION_DB/GRCh38.91.antisense_mRNA.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_10_NOT.$SAMPLENAME.bed"
  echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[9]}"
  ## intergenic region
  intersectBed -wo -s \
  -a "$OUTPUT/$SAMPLENAME/temp_10_NOT.$SAMPLENAME.bed" \
  -b "$ANNOTATION_DB/GRCh38.91.intergenic.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_11_intergenic.$SAMPLENAME.bed"
  #
  intersectBed -v -s \
  -a "$OUTPUT/$SAMPLENAME/temp_10_NOT.$SAMPLENAME.bed" \
  -b "$ANNOTATION_DB/GRCh38.91.intergenic.sorted.bed" \
  > "$OUTPUT/$SAMPLENAME/temp_11_NOT.$SAMPLENAME.bed"
  echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[10]}"
  date +"%d/%m/%Y %H:%M:%S" 
  echo "Finnished intersecting all samples."
 #
  ## Merge files
  date +"%d/%m/%Y %H:%M:%S" 
  echo "Merging all annotated files."
  cat "$OUTPUT/$SAMPLENAME/temp_1_rRNA.$SAMPLENAME.bed" \
  "$OUTPUT/$SAMPLENAME/temp_2_tRNA.$SAMPLENAME.bed" \
  "$OUTPUT/$SAMPLENAME/temp_3_snoRNA.$SAMPLENAME.bed" \
  "$OUTPUT/$SAMPLENAME/temp_4_snRNA.$SAMPLENAME.bed" \
  "$OUTPUT/$SAMPLENAME/temp_5_miRNA.$SAMPLENAME.bed" \
  "$OUTPUT/$SAMPLENAME/temp_6_mRNA_exons.$SAMPLENAME.bed" \
  "$OUTPUT/$SAMPLENAME/temp_7_mRNA_introns.$SAMPLENAME.bed" \
  "$OUTPUT/$SAMPLENAME/temp_8_lincRNA.$SAMPLENAME.bed" \
  "$OUTPUT/$SAMPLENAME/temp_9_repetitive_elements.$SAMPLENAME.bed" \
  "$OUTPUT/$SAMPLENAME/temp_10_antisense_RNA.$SAMPLENAME.bed" \
  "$OUTPUT/$SAMPLENAME/temp_11_intergenic.$SAMPLENAME.bed" \
  > "$OUTPUT/$SAMPLENAME/$SAMPLEname.all_annot.bed"
done
# Wait for all background jobs to finish before the script continues
wait

# gzip -f "$INPUT/$SAMPLE/$SAMPLEname.ann.bed"

## Remove temp
# rm "$INPUT/$smp"/temp_?_NOT.$SAMPLENAME.bed
date +"%d/%m/%Y %H:%M:%S" 
echo "Annotation script finished."

#####################################################################s
  
