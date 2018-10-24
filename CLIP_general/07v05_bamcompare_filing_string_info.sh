#!/bin/bash
############################################################################################
# CLIPseq bamCompare analysis of ABH8 CLIPseq and RNAseq
# Filing strand info, merging and sorting resulting files for further annotation
### Commands
####################################################################################################
cd /home/tomas/CEITEC_lab/ABH8/04_bamCompare
#
FW_BED=(*.FW.bed)
RV_BED=(*.RV.bed)
if [[ -f ${FW_BED[0]} ]]; then
	for a in *.FW.bed; do
		FILE=$a
		FILENAME=${a%.*}
    	date +"%d/%m/%Y %H:%M:%S"
		echo "Start filling strand info to file $FILE"
		awk -F'[\t]' -v OFS='\t' '{print $1, $2, $3, $4, $4, "+"}' $FILE > "$FILENAME".fill.bed
		date +"%d/%m/%Y %H:%M:%S"
		echo "Finnished filling strand info to file $FILE"
	done
else
	date +"%d/%m/%Y %H:%M:%S"
	echo "There is no *.FW.bed file for strand info filling!" && exit 1
fi
if [[ -f ${RV_BED[0]} ]]; then
	for b in *.RV.bed; do
		FILE2=$b
		FILENAME2=${b%.*}
    	date +"%d/%m/%Y %H:%M:%S"
		echo "Start filling strand info to file $FILE2"
		awk -F'[\t]' -v OFS='\t' '{print $1, $2, $3, $4, $4, "-"}' $FILE2 > "$FILENAME2".fill.bed
		date +"%d/%m/%Y %H:%M:%S"
		echo "Finnished filling strand info to file $FILE2"
	done
else
	date +"%d/%m/%Y %H:%M:%S"
	echo "There is no *.RV.bed file for strand info filling!" && exit 1
fi
#
wait
#
for c in *.FW.fill.bed; do
		FILE3=$c
		FILENAME3=${c%.*.*.*}

		if [[ -a "$FILENAME3".RV.fill.bed ]]; then
			date +"%d/%m/%Y %H:%M:%S"
			echo "Start combining $FILE3 and "$FILENAME3".RV.fill.bed files with filled strand info"
			cat "$FILE3" "$FILENAME3".RV.fill.bed |
			sortBed -i stdin > "$FILENAME3".filled.bed
			cho "Finnished combining $FILE3 and "$FILENAME3".RV.fill.bed files with filled strand info"
		else
			echo "There are no matching FW and RV files for merging!" && exit 1
		fi
done