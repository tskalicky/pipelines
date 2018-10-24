#!/bin/bash
# awk to print only fields from file containing mutations in unique tags .tag.uniq.mutation.txt
INPUTDIR="/home/skalda/Dropbox/CEITEC_lab/FTO/mapped_PCR_collapsed/"

cd $INPUTDIR
for z in *.tag.uniq.mutation.txt
do
	FILE=$z
	FILENAME=${z%.*}
	printf "Counting number of mutations, deletions and insertions for file $FILE \n"
	awk '{print $8,$9,$10,$11}' $FILE >> $FILENAME".summary.temp.txt"
	MUT=$(egrep -c '>' $FILENAME".summary.temp.txt")
	DEL=$(egrep -c '-' $FILENAME".summary.temp.txt")
	INS=$(egrep -c '+' $FILENAME".summary.temp.txt")
	printf "$MUT is summary number of mutattions in file $FILE \n" >> $FILENAME".summary.final.txt"
	printf "$DEL is summary number of deletions in file $FILE \n" >> $FILENAME".summary.final.txt"
	printf "$INS is summary number of insertions in file $FILE \n" >> $FILENAME".summary.final.txt"
done
