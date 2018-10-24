#!/bin/bash
# CLIP Seq pipeline
# Sorting and index creation for BED files after PCR colapsing and bedgraph creation
# 
# Requires IGVtools
############################################################################################
### Variables
APPENDIX=".rgb.bed"
APPENDIX2=".rgb.sorted.bed"
INPUT="/home/skalda/ownCloud/CEITEC_lab/METTL16/mapped/collapsed_PCR_duplicates/METTL16_UV1-2/"
INPUT2="/home/skalda/ownCloud/CEITEC_lab/METTL16/mapped/collapsed_PCR_duplicates/METTL16PAR1-2"
# Binaries
IGVTOOLS=$(which igvtools)

# Check if tools are installed
which $IGVTOOLS
# Commands
cd $INPUT
for a in *$APPENDIX
do
	FILE=$a
	FILENAME=${a%.*}
	echo "Now I am sorting file $FILE"
	$IGVTOOLS "sort" $FILE $FILENAME".sorted.bed"
	echo "Finnished sorting file $FILE"
	echo "Started indexing file $FILE"
	$IGVTOOLS "index" $FILENAME".sorted.bed"
done
