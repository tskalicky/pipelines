#!/bin/bash
INPUTDIR="/home/tomas/CEITEC_lab/ABH8/05_anotace/final/ABH8_CLIP1-3"
OUTPUTDIR="/home/tomas/CEITEC_lab/ABH8/05_anotace/final/ABH8_CLIP1-3/renamed_merged"
cd $INPUTDIR
mkdir $OUTPUTDIR
for a in *.bed
do
	FILE=$a
  	FILENAME=${a%.*}
  	# For spliced alignments (like STAR RNA) need to create BED12 format
  	# otherwise not only spliced reads on exons but the whole range will be 
  	# reported as interval !!!
  	date +"%d/%m/%Y %H:%M:%S" 
  	echo "Reshuffling fields in bed file $FILE"
	awk -F'[\t]' -v OFS='\t' '{print $1, $2, $3, $16, $25, $6}' $FILE > $OUTPUTDIR/"$FILENAME".ren.bed
	sortBed -i $OUTPUTDIR/"$FILENAME".ren.bed | \
	mergeBed -s -c 4,5,6 -o distinct,median,distinct -i stdin | \
	awk -F'[\t]' -v OFS='\t' '{print $1, $2, $3, $5, $6, $4}' > $OUTPUTDIR/"$FILENAME".rsm.bed
	date +"%d/%m/%Y %H:%M:%S" 
  	echo "Finnished reshuffling fields in bed file $FILE"
done

