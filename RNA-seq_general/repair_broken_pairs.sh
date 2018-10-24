#!/bin/bash
#PBS -l walltime=24:0:0 
#PBS -q default@wagap-pro.cerit-sc.cz
#PBS -l select=1:ncpus=12:mem=200gb:scratch_local=200gb
#PBS -j oe
#PBS -N repair_broken_pairs_01
#
## initialize the required application
module add samtools-1.4 # Samtools v1.6 are broken! 
module add python27-modules-gcc #required by multiQC
#
# Sorting and index creation for BAM files
#
# Requires samtools, multiQC
#
############################################################################################
### Variables
INPUT="/storage/brno3-cerit/home/tskalicky/Dis3L2/data/trimmed"
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/Dis3L2/data/trimmed/repaired_pairs"
MY_RAM=200 # Max RAM memory for Samtools sort
THREADS=$PBS_NUM_PPN
# Adding bbmap binaries into the PATH
# Will work ONLY on NFS4 connected servers
export PATH="/storage/brno3-cerit/home/tskalicky/tools/bbmap:$PATH"
#
# Binaries
SAMTOOLS=$(which samtools)
REPAIR=$(which repair.sh)
MULTIQC=$(which multiqc)
# Check the tools versions
which $SAMTOOLS
which $MULTIQC
which $REPAIR
#
####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
cd $INPUT
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
find "$INPUT" -maxdepth 2 -name "*.fastq.gz" -exec cp -vt "$SCRATCHDIR" {} +
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR
echo "Following files were copied to scratch:"
ls -Rc1
echo "Following files were copied to scratch:" >> /storage/brno3-cerit/home/tskalicky/Dis3L2/data/trimmed/log.txt #Debugging
ls -Rc1 >> /storage/brno3-cerit/home/tskalicky/Dis3L2/data/trimmed/log.txt #Debugging
####################################################################################################
## Commands

# for a in *.fastq.gz
# do
# 	SAMPLE=$a
# 	SAMPLENAME=${a%.*}
# 	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# 	echo "Now I am decompressing PE reads $SAMPLE"
# 	unpigz -dv -p $THREADS $SAMPLE &
# done
# # wait
# wait
# date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# echo "Done decompressing all PE reads"
# for a in *.fastq.gz; do unpigz -dv -p 4 $a; done

#
shopt nullglob # If set, Bash allows filename patterns which match no files to expand to a null string, rather than themselves.
IFS=$'\n' # split on newline only, needed for filling arrays with output from find
set -f    # disable globbing, needed for filling arrays with output from find
declare an array variables
declare -a KO1=($(find . -maxdepth 1 -name 'Dis3L2_KO1*.fastq.gz' -exec basename {} \; | sort -n))
declare -a KO2=($(find . -maxdepth 1 -name 'Dis3L2_KO2*.fastq.gz' -exec basename {} \; | sort -n))
declare -a KO3=($(find . -maxdepth 1 -name 'Dis3L2_KO3*.fastq.gz' -exec basename {} \; | sort -n))
declare -a WT1=($(find . -maxdepth 1 -name 'Dis3L2_WT1*.fastq.gz' -exec basename {} \; | sort -n))
declare -a WT2=($(find . -maxdepth 1 -name 'Dis3L2_WT2*.fastq.gz' -exec basename {} \; | sort -n))
declare -a WT3=($(find . -maxdepth 1 -name 'Dis3L2_WT3*.fastq.gz' -exec basename {} \; | sort -n))
declare -a ARRAY_NAMES=("KO1" "KO2" "KO3" "WT1" "WT2" "WT3")
declare -a ALL_ARRAYS=("${KO1[@]}" "${KO2[@]}" "${KO3[@]}" "${WT1[@]}" "${WT2[@]}" "${WT3[@]}")
set +f # enable globbing, because is needed for other parts of the script, like variable expansions !!!
# get length of an array
KO1_num="${#KO1[@]}"
KO2_num="${#KO2[@]}"
KO3_num="${#KO3[@]}"
WT1_num="${#WT1[@]}"
WT2_num="${#WT2[@]}"
WT3_num="${#WT3[@]}"
quantity="${#ALL_ARRAYS[@]}"
#
echo "There are $KO1_num KO1 samples that will be mapped."
echo "Sample names are: ${KO1[@]}"
echo "There are $KO2_num KO2_num samples that will be mapped."
echo "Sample names are: ${KO2[@]}"
echo "There are $KO3_num KO3_num samples that will be mapped."
echo "Sample names are: ${KO3[@]}"
echo "There are $WT1_num WT1_num samples that will be mapped."
echo "Sample names are: ${WT1[@]}"
echo "There are $WT2_num WT2_num samples that will be mapped."
echo "Sample names are: ${WT2[@]}"
echo "There are $WT3_num WT3_num samples that will be mapped."
echo "Sample names are: ${WT3[@]}"
### bbmap repair broken pairs
# bbmap repair.sh can handle gzip compressed files
for (( w=0;w<${quantity};w +=2 ));
do
	SAMPLE1="${ALL_ARRAYS[$w]}"
	SAMPLE2="${ALL_ARRAYS[$w+1]}"
	SAMPLENAME=${SAMPLE1%"_flexbar_"*}
	SAMPLENAME2=${SAMPLE2%"_flexbar_"*}
	if [[ -f $SAMPLE1 ]]; then
		date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
		echo "Now I am repairing pairs in PE reads $SAMPLE1 and $SAMPLE2 - alignment"
		$REPAIR in1=$SAMPLE1 in2=$SAMPLE2 out1="$SAMPLENAME"_fixed_1.fastq out2="$SAMPLENAME2"_fixed_2.fastq outs="$SAMPLENAME"_all_singletons.fastq repair
		echo "Done repairing pairs in PE reads $SAMPLE"
	else
		date +"%d/%m/%Y %H:%M:%S"
		echo "There is no $SAMPLE1 file for repairing pairs!" && exit 1
	fi
done
wait
# deleting original files
rm -v ${ALL_ARRAYS[@]}
#
# Packing results
for c in *.fastq
do
	SAMPLE3=$c
	SAMPLENAME3=${c%.*}
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Now I am compressing PE reads $SAMPLE3"
	pigz -v -p $THREADS $SAMPLE3
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Done compressing PE reads $SAMPLE3"
done
# wait
wait
####################################################################################################
### Finalize and copy results
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Finalize and copy results"
mkdir -p $OUTPUT_DIR
mv -v *.fastq.gz $OUTPUT_DIR
rm -rv $SCRATCHDIR/*
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"