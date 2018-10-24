#!/bin/bash
cd /storage/brno3-cerit/home/tskalicky/Dis3L2/Dasa_spikein/data/trimmed
# cd /home/tomas/CEITEC_lab/Dis3L2/Dasa_spikein/data/trimmed
shopt nullglob # If set, Bash allows filename patterns which match no files to expand to a null string, rather than themselves.
# declare an array variables
IFS=$'\n' # split on newline only
set -f    # disable globbing
declare -a OAT_CYTO=($(find . -name '*OAT_cyto*.fastq.gz' -exec basename {} \; | sort -n))
declare -a OAT_NONFRAC=($(find . -name '*OAT_nonfrac*.fastq.gz' -exec basename {} \; | sort -n))
declare -a OAT_NUCL=($(find . -name '*OAT_nucl*.fastq.gz' -exec basename {} \; | sort -n))
declare -a U12_CYTO=($(find . -name '*U12_cyto*.fastq.gz' -exec basename {} \; | sort -n))
declare -a U12_NONFRAC=($(find . -name '*U12_nonfrac*.fastq.gz' -exec basename {} \; | sort -n))
declare -a U12_NUCL=($(find . -name '*U12_nucl*.fastq.gz' -exec basename {} \; | sort -n))
declare -a SAMPLEarray=("OAT_CYTO" "OAT_NONFRAC" "OAT_NUCL" "U12_CYTO" "U12_NONFRAC" "U12_NUCL")
declare -a ALL_ARRAYS=("${OAT_CYTO[@]}" "${OAT_NONFRAC[@]}" "${OAT_NUCL[@]}" "${U12_CYTO[@]}" "${U12_NONFRAC[@]}" "${U12_NUCL[@]}" )
# get length of an array
OAT_cyto_num="${#OAT_CYTO[@]}"
OAT_nonfrac_num="${#OAT_NONFRAC[@]}"
OAT_nucl_num="${#OAT_NUCL[@]}"
U12_cyto_num="${#U12_CYTO[@]}"
U12_nonfrac_num="${#U12_NONFRAC[@]}"
U12_nucl_num="${#U12_NUCL[@]}"
quantity="${#ALL_ARRAYS[@]}"
#
echo "There are $OAT_cyto_num OAT_CYTO samples that will be mapped."
echo "Sample names are: ${OAT_CYTO[@]}"
echo "There are $OAT_nonfrac_num OAT_nonfrac_num samples that will be mapped."
echo "Sample names are: ${OAT_NONFRAC[@]}"
echo "There are $OAT_nucl_num OAT_nucl_num samples that will be mapped."
echo "Sample names are: ${OAT_NUCL[@]}"
echo "There are $U12_cyto_num U12_cyto_num samples that will be mapped."
echo "Sample names are: ${U12_CYTO[@]}"
echo "There are $U12_nonfrac_num U12_nonfrac_num samples that will be mapped."
echo "Sample names are: ${U12_NONFRAC[@]}"
echo "There are $U12_nucl_num U12_nucl_num samples that will be mapped."
echo "Sample names are: ${U12_NUCL[@]}"
### Alignment
for (( w=0;w<${quantity};w +=2 ));
do
	SAMPLENAME="${ALL_ARRAYS[$w]}"
	SAMPLENAME2="${ALL_ARRAYS[$w+1]}"
	echo "$SAMPLENAME"
	echo "$SAMPLENAME2"
	if [[ -f $SAMPLENAME ]]; then
		date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
		echo "Now I am processing PE reads $SAMPLENAME and $SAMPLENAME2 - alignment"
		echo ${SAMPLENAME%"_flexbar_1.fastq.gz"*}.${GTF%.*}
		date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
		echo "Done processing PE reads $SAMPLENAME and $SAMPLENAME2 - alignment"
	else
		date +"%d/%m/%Y %H:%M:%S"
		echo "There is no $SAMPLENAME file for mapping!" && exit 1
	fi
done
wait