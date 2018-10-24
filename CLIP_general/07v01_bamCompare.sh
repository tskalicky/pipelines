#!/bin/bash
# CLIP Seq pipeline
# Sorting and index creation for BED files after PCR colapsing and bedgraph creation
#
# 
# Requires IGVtools
############################################################################################
### Variables
APPENDIX=".rgb.bed"
APPENDIX2=".rgb.sorted.bed"
THREADS=$(nproc)
# Need to export Anaconda installation path because it can NOT be in login shell profiles 
# like .bashrc or .profile. If present in $PATH it will crash KDE !!!
export PATH="$PATH:/home/tomas/anaconda2/bin"
INPUT="/home/tomas/CEITEC_lab/ABH8/mapping"
CLIP="/home/tomas/CEITEC_lab/ABH8/mapping/ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.bam"
RNA-SEQ="/home/tomas/CEITEC_lab/ABH8/mapping/ABH8_RNAseq1-3_GRCh38.pcr_dedupl.rg.merged.bam"
OUTPUTDIR="/home/tomas/CEITEC_lab/ABH8/bamCompare"
# Binaries
IGVTOOLS=$(which igvtools)
BAMCOMPARE=$(which bamCompare)

# Check if tools are installed
which $IGVTOOLS
which $BAMCOMPARE
# Commands
cd $INPUT
# for a in *.bam.gz
# do
# 	FILE=$a
# 	FILENAME=${a%.*.*}
# 	date +"%d/%m/%Y %H:%M:%S" 
# 	echo "Now I am decompressing library $FILE"
# 	unpigz -v -p $THREADS $FILE
# 	echo "Done decompressing PE reads $FILE"
# done
# # wait
# wait
# for b in *.bam
# do
# 	FILE=$b
# 	FILENAME=${b%.*}
# 	date +"%d/%m/%Y %H:%M:%S" 
# 	echo "Now I am sorting file $FILE"
# 	igvtools "sort" $FILE $FILENAME".sorted.bam"
# 	echo "Finnished sorting file $FILE"
# 	date +"%d/%m/%Y %H:%M:%S" 
# 	echo "Started indexing file $FILE"
# 	igvtools "index" $FILENAME".sorted.bam"
# done
# wait
#
## Compare 2 datasets
# mkdir /home/tomas/CEITEC_lab/ABH8/bamCompare
date +"%d/%m/%Y %H:%M:%S" 
echo "Now I am comparing files using bamCompare"
bamCompare \
-b1 ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.sorted.bam \
-b2 ABH8_RNAseq1-3_GRCh38.pcr_dedupl.rg.merged.sorted.bam \
--numberOfProcessors 3 \
--binSize 20 \
--scaleFactorsMethod 'None' \
--operation 'ratio' \
--pseudocount 1.0 \
--normalizeUsing 'RPKM' \
--skipNonCoveredRegions \
--verbose \
-of bedgraph \
-o $OUTPUTDIR/ABH8_CLIP_vs_RNAseq_enrichment_compare_RPKM_nosmooth.bed
date +"%d/%m/%Y %H:%M:%S" 
echo "Now I am comparing files using bamCompare"

# bamCompare \
# -b1 ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.sorted.bam \
# -b2 ABH8_RNAseq1-3_GRCh38.pcr_dedupl.rg.merged.sorted.bam \
# --numberOfProcessors 3 \
# --binSize 20 \
# --scaleFactorsMethod 'None' \
# --operation 'ratio' \
# --pseudocount 1.0 \
# --normalizeUsing 'RPKM' \
# --smoothLength 60 \
# --skipNonCoveredRegions \
# --verbose \
# -of bedgraph \
# -o $OUTPUTDIR/ABH8_CLIP_vs_RNAseq_enrichment_compare_RPKM_normal.bed

# bamCompare \
# -b1 ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.sorted.bam \
# -b2 ABH8_RNAseq1-3_GRCh38.pcr_dedupl.rg.merged.sorted.bam \
# --numberOfProcessors 3 \
# --binSize 20 \
# --scaleFactorsMethod 'SES' \
# --sampleLength 10000 \
# --numberOfSamples 100000 \
# --operation 'log2' \
# --pseudocount 1.0 \
# --smoothLength 60 \
# --skipNonCoveredRegions \
# --verbose \
# -of bedgraph \
# -o $OUTPUTDIR/ABH8_CLIP_vs_RNAseq_log2_enrichment_compare.bed

# bamCompare \
# -b1 ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.sorted.bam \
# -b2 ABH8_RNAseq1-3_GRCh38.pcr_dedupl.rg.merged.sorted.bam \
# --numberOfProcessors 3 \
# --binSize 20 \
# --scaleFactorsMethod 'SES' \
# --sampleLength 10000 \
# --numberOfSamples 100000 \
# --operation 'ratio' \
# --pseudocount 1.0 \
# --smoothLength 60 \
# --skipNonCoveredRegions \
# --verbose \
# -of bedgraph \
# -o $OUTPUTDIR/ABH8_CLIP_vs_RNAseq_enrichment_compare.bed