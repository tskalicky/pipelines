#!/bin/bash

homeSmp="$HOME/Projects/dis3l2/samples"
homeDb="$HOME/Projects/dis3l2/annotation/databases/output"

if [ $# -ne 2 ]; then
  echo "Please set two arguments: sample file"
  exit 0
fi

smp=$1   # Sample to analyze
file=$2  # File to analyze
filename="${file%.*}"
filename="${filename%.*}"

echo "Annotating $smp/$file"

## tRNA
intersectBed -wo -s \
  -a "$homeSmp/$smp/$file" \
  -b "$homeDb/db_tRNA.bed.gz" \
   > "$homeSmp/$smp/temp_1_tRNA.$file.bed"

intersectBed -v -s \
  -a "$homeSmp/$smp/$file" \
  -b "$homeDb/db_tRNA.bed.gz" \
   > "$homeSmp/$smp/temp_1_NOT.$file.bed"

## miRNA
intersectBed -wo -s \
  -a "$homeSmp/$smp/temp_1_NOT.$file.bed" \
  -b "$homeDb/db_miRNA.bed.gz" \
   > "$homeSmp/$smp/temp_2_miRNA.$file.bed"

intersectBed -v -s \
  -a "$homeSmp/$smp/temp_1_NOT.$file.bed" \
  -b "$homeDb/db_miRNA.bed.gz" \
   > "$homeSmp/$smp/temp_2_NOT.$file.bed"

## snRNA
intersectBed -wo -s \
  -a "$homeSmp/$smp/temp_2_NOT.$file.bed" \
  -b "$homeDb/db_snRNA.bed.gz" \
   > "$homeSmp/$smp/temp_3_snRNA.$file.bed"

intersectBed -v -s \
  -a "$homeSmp/$smp/temp_2_NOT.$file.bed" \
  -b "$homeDb/db_snRNA.bed.gz" \
   > "$homeSmp/$smp/temp_3_NOT.$file.bed"

## rRNA
intersectBed -wo -s \
  -a "$homeSmp/$smp/temp_3_NOT.$file.bed" \
  -b "$homeDb/db_rRNA.bed.gz" \
   > "$homeSmp/$smp/temp_4_rRNA.$file.bed"

intersectBed -v -s \
  -a "$homeSmp/$smp/temp_3_NOT.$file.bed" \
  -b "$homeDb/db_rRNA.bed.gz" \
   > "$homeSmp/$smp/temp_4_NOT.$file.bed"

## miscRNA
intersectBed -wo -s \
  -a "$homeSmp/$smp/temp_4_NOT.$file.bed" \
  -b "$homeDb/db_miscRNA.bed.gz" \
   > "$homeSmp/$smp/temp_5_miscRNA.$file.bed"

intersectBed -v -s \
  -a "$homeSmp/$smp/temp_4_NOT.$file.bed" \
  -b "$homeDb/db_miscRNA.bed.gz" \
   > "$homeSmp/$smp/temp_5_NOT.$file.bed"

## snoRNA
intersectBed -wo -s \
  -a "$homeSmp/$smp/temp_5_NOT.$file.bed" \
  -b "$homeDb/db_snoRNA.bed.gz" \
   > "$homeSmp/$smp/temp_6_snoRNA.$file.bed"

intersectBed -v -s \
  -a "$homeSmp/$smp/temp_5_NOT.$file.bed" \
  -b "$homeDb/db_snoRNA.bed.gz" \
   > "$homeSmp/$smp/temp_6_NOT.$file.bed"

## mRNA = CDS + UTR
## UTR
intersectBed -wo -s \
  -a "$homeSmp/$smp/temp_6_NOT.$file.bed" \
  -b "$homeDb/db_mRNA_UTR.bed.gz" \
   > "$homeSmp/$smp/temp_7_mRNA_UTR.$file.bed"

intersectBed -v -s \
  -a "$homeSmp/$smp/temp_6_NOT.$file.bed" \
  -b "$homeDb/db_mRNA_UTR.bed.gz" \
   > "$homeSmp/$smp/temp_7_NOT.$file.bed"

## CDS
intersectBed -wo -s \
  -a "$homeSmp/$smp/temp_7_NOT.$file.bed" \
  -b "$homeDb/db_mRNA_CDS.bed.gz" \
   > "$homeSmp/$smp/temp_8_mRNA_CDS.$file.bed"

intersectBed -v -s \
  -a "$homeSmp/$smp/temp_7_NOT.$file.bed" \
  -b "$homeDb/db_mRNA_CDS.bed.gz" \
   > "$homeSmp/$smp/temp_8_NOT.$file.bed"

## mRNA transcripts (only for histone stem-loop heatmap)
intersectBed -wo -s \
  -a "$homeSmp/$smp/temp_6_NOT.$file.bed" \
  -b "$homeDb/db_mRNA_transcript.bed.gz" \
   | gzip \
   > "$homeSmp/$smp/$filename.ann_transcripts.bed.gz"

## Annotate lincRNA
intersectBed -wo -s \
  -a "$homeSmp/$smp/temp_8_NOT.$file.bed" \
  -b "$homeDb/db_lincRNA.bed.gz" \
   > "$homeSmp/$smp/temp_9_lincRNA.$file.bed"

## Merge files
cat "$homeSmp/$smp/temp_1_tRNA.$file.bed"     \
    "$homeSmp/$smp/temp_2_miRNA.$file.bed"    \
    "$homeSmp/$smp/temp_3_snRNA.$file.bed"    \
    "$homeSmp/$smp/temp_4_rRNA.$file.bed"     \
    "$homeSmp/$smp/temp_5_miscRNA.$file.bed"  \
    "$homeSmp/$smp/temp_6_snoRNA.$file.bed"   \
    "$homeSmp/$smp/temp_7_mRNA_UTR.$file.bed" \
    "$homeSmp/$smp/temp_8_mRNA_CDS.$file.bed" \
    "$homeSmp/$smp/temp_9_lincRNA.$file.bed"  \
  > "$homeSmp/$smp/$filename.ann.bed"

gzip -f "$homeSmp/$smp/$filename.ann.bed"

## Remove temp
rm "$homeSmp/$smp"/temp_?_*.$file.bed




