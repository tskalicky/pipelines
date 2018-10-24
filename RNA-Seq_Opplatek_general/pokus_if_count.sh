#!/bin/bash
cd /home/skalda/ownCloud/CEITEC_lab/FTO_projekt/scripts
echo "Checking if we have biological replicates"
shopt -s nullglob
# declare an array variables
declare -a BEDNAMES=( *.tag.uniq.bed )
# get length of an array
quantity=${#BEDNAMES[@]}
#
if [[ $quantity -gt 1 ]]; then
	echo "Number of biological replicates is:" $quantity
	echo "Start concatenation of biological replicates"
	# declare an array variables
	declare -a COLORarray=("204,0,0" "0,204,0" "0,0,204" "204,204,0" "0,204,204" "204,0,204")
	# get length of an array
	length2=${#COLORarray[@]}
	# use for loop to read all values and indexes
	for (( c=0;c<${quantity}+1;c++ ));
	do
		SAMPLE=${BEDNAMES[$c]}
		echo $SAMPLE
		COLOR=${COLORarray[$c]}
		SAMPLENAME4=${SAMPLE%.tag*}
		#
		echo bed2rgb.pl -v -col "$COLOR" ${BEDNAMES[$c]} $SAMPLENAME4".tag.uniq.rgb.bed" #as one example
		#
		## repeat the above step for all other uniq.bed files to generate rgb.bed files, 
		## but use different rgb colors "x,x,x" (see http://www.rapidtables.com/web/color/RGB_Color.htm or other charts)
		#cat *tag.uniq.rgb.bed > $SAMPLENAME4".pool.tag.uniq.rgb.bed"
		#cat *tag.uniq.mutation.txt > $SAMPLENAME4".pool.tag.uniq.mutation.txt" 
	done
else
    echo "There is $quantity biological replicate! Skipping concatenation and replicate color taging"
fi