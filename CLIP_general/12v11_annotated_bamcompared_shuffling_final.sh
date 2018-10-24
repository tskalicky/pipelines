#!/bin/bash
## Final sorting, merging and shuffling on local data for statistics
cd /home/tomas/CEITEC_lab/ABH8/05_anotace/annot_on_bamCompared/Intersect_ABH8_bamcompared_and_PIPE-CLIP_enrichedClusters_n50/renamed_merged
for a in *.rsm.bed; # I know, I know programmers COUNT from zero :-D
do
	SAMPLE=$a
	SAMPLENAME=${a%.*}
	date +"%d/%m/%Y %H:%M:%S"
	echo "Sorting and merging renamed file $SAMPLE."
	# sortBed -i $SAMPLE | \
	# mergeBed -s -c 4,5,6,7,8 -o distinct,median,distinct,median,distinct -i stdin > temp_2_tRNA.rsm.bed
	# awk -F'[\t]' -v OFS='\t' '{print $1, $2, $3, $5, $6, $7, $8, $9 }' temp_2_tRNA.rsm.bed > temp_2_tRNA.rsm.shuf3.bed
	awk -F'[\t]' '{print $1 "\t" $2 "\t" $3 "\t" $5 "\t" $6 "," $8 "," $9 "\t" $7 }' $SAMPLE > $SAMPLENAME.shuf3.bed
	date +"%d/%m/%Y %H:%M:%S"
	echo "Finnished Sorting and merging renamed file $SAMPLE."
done
#
wait
