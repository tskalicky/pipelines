#!/bin/bash
############################################################################################
# Uridylated RNAseq reads processing.
# This script will take sequence names from reads that were originally uridylated, trimmed and mapped back to human genome,
# and use them as a query to search in the original unmapped and not Ttrimmed data.
# Output will be the original uridylated sequences in fasta format.
# Requires samtools and bbmap
############################################################################################
### Variables
THREADS="6"
UPSCALE_TRIMMED_fasta="/home/tomas/ownCloud/CEITEC_lab/Dis3L2/Dasa_spikein/uridylation/UPSCALE/mapping_trimmedT/extract_mapped/fasta"
NEW_TRIMMED_fasta="/home/tomas/ownCloud/CEITEC_lab/Dis3L2/Dasa_spikein/uridylation/NEW/mapping_trimmedT/extract_mapped/fasta"
URIDYL_NEW_fasta="/home/tomas/ownCloud/CEITEC_lab/Dis3L2/Dasa_spikein/uridylation/NEW/sorted_original"
URIDYL_UPSCALE_fasta="/home/tomas/ownCloud/CEITEC_lab/Dis3L2/Dasa_spikein/uridylation/UPSCALE/sorted_original/fasta"
OUTPUT_DIR="/home/tomas/ownCloud/CEITEC_lab/Dis3L2/Dasa_spikein/uridylation/final_uridylated_untrimmed"
#
export PATH="$PATH:/home/tomas/anaconda2/bin/"
echo "PATH after modification is:"
echo $PATH | tr ":" "\n" | nl
# Binaries
SAMTOOLS=$(which samtools)
# Check the tools versions
which $SAMTOOLS
which python
####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
# mkdir -p $OUTPUT_DIR
# trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
# find "$NEW_TRIMMED_fasta" -maxdepth 1 -name "trimT_DIS3l2*uniq_mapped.fa" -exec cp -vt "$OUTPUT_DIR" {} +
# find "$UPSCALE_TRIMMED_fasta" -maxdepth 1 -name "trimT_upscale_DIS3l2*uniq_mapped.fa" -exec cp -vt "$OUTPUT_DIR" {} +
# find "$URIDYL_NEW_fasta" -maxdepth 1 -name "Uridyl_DIS3l2*unmapped.fa.gz" -exec cp -vt "$OUTPUT_DIR" {} +
# find "$URIDYL_UPSCALE_fasta" -maxdepth 1 -name "Uridyl_upscale_DIS3l2*unmapped.fa.gz" -exec cp -vt "$OUTPUT_DIR" {} +
# cp -av "/home/tomas/ownCloud/git/pipelines/RNA-seq_general/Dis3L2_uridylation_project/13_extract_untrimmed_seq_fasta_v2.py" $OUTPUT_DIR
cd $OUTPUT_DIR
#
if [ ! -d "$OUTPUT_DIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "OUTPUT_DIR path is:" $OUTPUT_DIR
echo "Following files were copied to scratch:"
ls -Rc1
# Compress output files
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Compressing results"
for a in *.fa.gz
do
	SAMPLE=$a
	SAMPLENAMEa=${a%.*}
	echo "Now I am decompressing reads $SAMPLE"
	unpigz -v -p $THREADS $SAMPLE
	echo "Done decompressing PE reads $SAMPLE"
done
#wait
wait
#
shopt nullglob # If set, Bash allows filename patterns which match no files to expand to a null string, rather than themselves.
IFS=$'\n' # split on newline only, needed for filling arrays with output from find
set -f    # disable globbing, needed for filling arrays with output from find
# declare an array variables
declare -a trimT_DIS3l2=($(find . -maxdepth 1 -name 'trimT_DIS3l2*uniq_mapped.fa' -exec basename {} \; | sort -n))
declare -a trimT_UPSCALE=($(find . -maxdepth 1 -name 'trimT_upscale_DIS3l2*uniq_mapped.fa' -exec basename {} \; | sort -n))
declare -a ORIG_DIS3L2=($(find . -maxdepth 1 -name 'Uridyl_DIS3l2*unmapped.fa' -exec basename {} \; | sort -n))
declare -a ORIG_UPSCALE=($(find . -maxdepth 1 -name 'Uridyl_upscale_DIS3l2*unmapped.fa' -exec basename {} \; | sort -n))
declare -a ARRAY_NAMES=("trimT_DIS3l2" "ORIG_DIS3L2" "trimT_UPSCALE" "ORIG_UPSCALE")
declare -a ALL_ARRAYS=("${trimT_DIS3l2[@]}" "${ORIG_DIS3L2[@]}" "${trimT_UPSCALE[@]}" "${ORIG_UPSCALE[@]}")
set +f # enable globbing, because is needed for other parts of the script, like variable expansions !!!
# get length of an array
trimT_DIS3l2_num="${#trimT_DIS3l2[@]}"
ORIG_DIS3L2_num="${#ORIG_DIS3L2[@]}"
trimT_UPSCALE_num="${#trimT_UPSCALE[@]}"
ORIG_UPSCALE_num="${#ORIG_UPSCALE[@]}"
quantity="${#ALL_ARRAYS[@]}"
lib_count="${#ARRAY_NAMES[@]}"
#
echo "There are $trimT_DIS3l2_num trimT_DIS3L2 samples that will be mapped."
echo Sample_names_are: ${trimT_DIS3l2[@]} | tr " " "\n" | nl
echo "There are $ORIG_DIS3L2_num ORIG_DIS3L2 samples that will be mapped."
echo Sample_names_are: ${ORIG_DIS3L2[@]} | tr " " "\n" | nl
echo "There are $trimT_UPSCALE_num trimT_UPSCALE samples that will be mapped."
echo Sample_names_are: ${trimT_UPSCALE[@]} | tr " " "\n" | nl
echo "There are $ORIG_UPSCALE_num ORIG_UPSCALE samples that will be mapped."
echo Sample_names_are: ${ORIG_UPSCALE[@]} | tr " " "\n" | nl
echo "There are $quantity SAMPLES that will be processed."
echo Sample_names_are: ${ALL_ARRAYS[@]} | tr " " "\n" | nl

for (( w=0;w<=${trimT_DIS3l2_num};w +=1 ));
do
	SAMPLE="${trimT_DIS3l2[$w]}"
	SAMPLE2="${ORIG_DIS3L2[$w]}"
	if [[ -f $SAMPLE ]]; then
		date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
		echo "Now I am processing uridylated samples $SAMPLE vs $SAMPLE2"
		python 13_extract_untrimmed_seq_fasta_v2.py $SAMPLE $SAMPLE2 
		date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
		echo "Done processing uridylated samples $SAMPLE vs $SAMPLE2"
	else
		date +"%d/%m/%Y %H:%M:%S"
		echo "There is no $SAMPLE file for processing!" && exit 1
	fi
done &
#
for (( w=0;w<=${trimT_UPSCALE_num};w +=1 ));
do
	SAMPL3="${trimT_UPSCALE[$w]}"
	SAMPLE4="${ORIG_UPSCALE[$w]}"
	if [[ -f $SAMPLE3 ]]; then
		date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
		echo "Now I am processing uridylated samples $SAMPLE3 vs $SAMPLE4"
		python 13_extract_untrimmed_seq_fasta_v2.py $SAMPLE3 $SAMPLE4 
		date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
		echo "Done processing uridylated samples $SAMPLE3 vs $SAMPLE4"
	else
		date +"%d/%m/%Y %H:%M:%S"
		echo "There is no $SAMPLE3 file for processing!" && exit 1
	fi
done
wait
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"