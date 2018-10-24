#!/bin/bash
#Create MERGED Databses
# -s = force merging only on the same strand
sortBed -i GRCh38.91.mRNA_exons.bed | \
mergeBed -s -c 4,5,6 -o distinct,distinct,distinct -i stdin > GRCh38.91.mRNA_exons.merged.sorted.bed
# Manual Intersecting of not assigned leftover temp_11_NOT.$SAMPLENAME.bed again against exons
intersectBed -wo \
-a "/home/tomas/CEITEC_lab/ABH8/anotace/ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged/temp_11_NOT.ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.bed" \
-b "/home/tomas/ownCloud/CEITEC_lab/genomes/annotation/my_annot_DB/ensembl/bed/GRCh38.91.mRNA_exons.merged.bed" \
> "/home/tomas/CEITEC_lab/ABH8/anotace/ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged/test/temp_12_exons_ABH8_CLIP1-3_GRCh38.bed"
#
intersectBed -v \
-a "/home/tomas/CEITEC_lab/ABH8/anotace/ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged/test/temp_12_exons_ABH8_CLIP1-3_GRCh38.bed" \
-b "/home/tomas/ownCloud/CEITEC_lab/genomes/annotation/my_annot_DB/ensembl/bed/GRCh38.91.mRNA_exons.merged.bed" \
> "/home/tomas/CEITEC_lab/ABH8/anotace/ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged/test/temp_12_NOT_exons.merged_ABH8_CLIP1-3_GRCh38.bed"