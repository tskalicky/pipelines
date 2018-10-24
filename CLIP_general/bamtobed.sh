#!/bin/bash
for a in *.bam.gz
do
  FILE=$a
  FILENAME=${a%.*.*}
  # For spliced alignments (like STAR RNA) need to create BED12 format
  # otherwise not only spliced reads on exons but the whole range will be 
  # reported as interval !!!
  gzip -dc $FILE | bamToBed -bed12 -i stdin > $FILENAME.bed
  sortBed -i $FILENAME.bed > $FILENAME.sorted.bed
done
wait
#
#bamToBed -bed12 -i ABH8_RNAseq1-3_GRCh38.pcr_dedupl.rg.merged.bam > ABH8_RNAseq1-3_GRCh38.pcr_dedupl.rg.merged.bed
#sortBed -i ABH8_RNAseq1-3_GRCh38.pcr_dedupl.rg.merged.bed > ABH8_RNAseq1-3_GRCh38.pcr_dedupl.rg.merged.sorted.bed
