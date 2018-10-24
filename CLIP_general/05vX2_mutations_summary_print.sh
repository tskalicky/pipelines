#!/bin/bash
# awk to print only fields from file containing mutations in unique tags .tag.uniq.mutation.txt
INPUTDIR="/home/skalda/Dropbox/CEITEC_lab/FTO/mapped_PCR_collapsed/"

cd $INPUTDIR
for z in *.tag.uniq.mutation.txt
do
	FILE=$z
	FILENAME=${z%.*.*}
	printf "Counting number of mutations, deletions and insertions for file $FILE \n"
	awk '{if($9=="-") {print $0}}' $FILE | cut 1-6 >> $FILENAME".del.mutations.txt"
	awk '{if($9==">") {print $0}}' $FILE | cut 1-6 >> $FILENAME".sub.mutations.txt"
	awk '{if($9=="+") {print $0}}' $FILE | cut 1-6 >> $FILENAME".ins.mutations.txt"
	SUB=$(egrep -c '>' $FILENAME".sub.mutations.txt")
	DEL=$(egrep -c '-' $FILENAME".del.mutations.txt")
	INS=$(egrep -c '+' $FILENAME".ins.mutations.txt")
	#awk '{print $8,$9,$10,$11}' $FILE >> $FILENAME".summary.temp.txt"
	#SUB=$(egrep -c '>' $FILENAME".summary.temp.txt")
	#DEL=$(egrep -c '-' $FILENAME".summary.temp.txt")
	#INS=$(egrep -c '+' $FILENAME".summary.temp.txt")
	printf "$SUB is summary number of mutattions in file $FILE \n" >> $FILENAME".mutation.summary.txt"
	printf "$DEL is summary number of deletions in file $FILE \n" >> $FILENAME".mutation.summary.txt"
	printf "$INS is summary number of insertions in file $FILE \n" >> $FILENAME".mutation.summary.txt"
done
