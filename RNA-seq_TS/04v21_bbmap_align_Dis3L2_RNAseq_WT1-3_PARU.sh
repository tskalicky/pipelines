#!/bin/bash
#PBS -l nodes=1:ppn=20
#PBS -l mem=200gb
#PBS -l walltime=96:00:00
#PBS -k oe
#PBS -N star_Dis3L2_RNAseq_WT1-3
#PBS -M tomas.skalicky@seznam.cz
#PBS -m abe
#
#PARU KRTECEK server is using TORQUE scheduling system !!!
#
## initialize the required application
# For 1 library you need cca 35GB RAM
#
# Alignment - PE RNA-Seq
# Requires STAR, pigz, samtools, multiQC, bedGraphToBigWig
#
# Note: do not use modified GTF (added features) to alignment as it can cause issues later on
############################################################################################
### Variables
INPUT_DIR="/home/users/tskalicky/CEITEC/Dis3L2/data/trimmed/WT"
MY_GENOME="/home/users/tskalicky/CEITEC/genomes/human/ensembl91/GRCh38.dna.primary_assembly.fa.gz"
MY_GTF="/home/users/tskalicky/CEITEC/genomes/human/ensembl91/GRCh38.91.gtf.gz"
OUTPUT_DIR="/home/users/tskalicky/CEITEC/Dis3L2/mapping/WT"
# Adding RSEM and bbmap binaries into the PATH
# Need to export system PATH on Krtecek server or bash commands wil NOT work111 :-/
# Will work ONLY on NFS4 connected servers
# export PATH="/storage/brno3-cerit/home/tskalicky/tools/RSEM-1.3.0/bin:$PATH"
echo "Default PATH is:"
echo "$PATH"
export PATH="/bin/:$PATH"
export PATH="/usr/bin/:$PATH"
export PATH="$PATH:/home/users/tskalicky/anaconda2/bin"
echo "PATH after modification is:"
echo "$PATH"
#
# Set number of CPU
THREADS=$PBS_NUM_PPN
# Set MAX ram
MY_RAM="200" # SET THE AMOUNT OF RAM HERE !!!
RAM=$[$MY_RAM-1] # lower it by one to leave some space
MEM_LIMIT=$RAM # $[$RAM/$THREADS] # In case we need GB per thread
#
# Binaries
UNPIGZ=$(which unpigz)
STAR=$(which STAR)
BBMAP=$(which bbmap.sh)
# SAMTOOLS=$(which samtools)
SAMTOOLS="/home/users/tskalicky/anaconda2/bin/samtools"
# BEDGRAPH2BIGWIG=$(which bedGraphToBigWig)
MULTIQC=$(which multiqc)
# FASTQ2COLLAPSE=$(which fastq2collapse.pl)
#
# Check the tools versions
echo "Checking required tools:"
#which $STAR
which $BBMAP
which $SAMTOOLS
#which $BEDGRAPH2BIGWIG
which $MULTIQC
which $UNPIGZ
#
####################################################################################################
# copy input data
# use cp -avr when copying directories
mkdir -p $OUTPUT_DIR
#
cp -av $MY_GENOME $OUTPUT_DIR
find "$INPUT_DIR" -maxdepth 1 -name *.fastq.gz -exec cp -avt "$OUTPUT_DIR" {} +
# cp -avr $GENOME_INDEX $OUTPUT_DIR/ 
cd $OUTPUT_DIR
#
if [ ! -d "$OUTPUT_DIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if output directory is created
echo "OUTPUT_DIR path is:" $OUTPUT_DIR
echo "Following files were copied to scratch:"
ls -Rc1
#
####################################################################################################
###
# Setting Variables
shopt nullglob # If set, Bash allows filename patterns which match no files to expand to a null string, rather than themselves.
IFS=$'\n' # split on newline only, needed for filling arrays with output from find
set -f    # disable globbing, ONLY needed for filling arrays with output from find
# declare an array variables
# declare -a KO1=($(find . -maxdepth 1 -name  'Dis3L2_KO1*.fastq.gz' -exec basename {} \; | sort -n))
# declare -a KO2=($(find . -maxdepth 1 -name  'Dis3L2_KO2*.fq.gz' -exec basename {} \; | sort -n))
# declare -a KO3=($(find . -maxdepth 1 -name  'Dis3L2_KO3*.fq.gz' -exec basename {} \; | sort -n))
declare -a WT1=($(find . -maxdepth 1 -name  'Dis3L2_WT1*.fq.gz' -exec basename {} \; | sort -n))
declare -a WT2=($(find . -maxdepth 1 -name  'Dis3L2_WT2*.fq.gz' -exec basename {} \; | sort -n))
declare -a WT3=($(find . -maxdepth 1 -name  'Dis3L2_WT3*.fq.gz' -exec basename {} \; | sort -n))
# declare -a ARRAY_NAMES=("KO1" "KO2" "KO3" "WT1" "WT2" "WT3")
# declare -a ALL_ARRAYS=("${KO1[@]}" "${KO2[@]}" "${KO3[@]}" "${WT1[@]}" "${WT2[@]}" "${WT3[@]}" )
declare -a ARRAY_NAMES=("WT1" "WT2" "WT3")
declare -a ALL_ARRAYS=("${WT1[@]}" "${WT2[@]}" "${WT3[@]}")
set +f # enable globbing, because is needed for other parts of the script, like variable expansions !!!
# get length of an array
# KO1_num="${#KO1[@]}"
# KO2_num="${#KO2[@]}"
# KO3_num="${#KO3[@]}"
WT1_num="${#WT1[@]}"
WT2_num="${#WT2[@]}"
WT3_num="${#WT3[@]}"
quantity="${#ALL_ARRAYS[@]}"
#
####################################################################################################
### Genome and annotation preparation
GENOME_NAME=$(basename $MY_GENOME)
GTF_NAME=$(basename $MY_GTF)

GENOME=${GENOME_NAME%.*}
GTF=${GTF_NAME%.*}
#
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Decompressing files:"
# $UNPIGZ -v -p $THREADS -d $GENOME_NAME
# unpigz -v -p $THREADS -d $GTF_NAME
# gzip -vd $GENOME_NAME
#
### Commands 
## BBmap - Alignment
# echo "There are $KO1_num KO1 samples that will be mapped."
# echo "Sample names are: ${KO1[@]}"
# echo "There are $KO2_num KO2_num samples that will be mapped."
# echo "Sample names are: ${KO2[@]}"
# echo "There are $KO3_num KO3_num samples that will be mapped."
# echo "Sample names are: ${KO3[@]}"
echo "There are $WT1_num WT1_num samples that will be mapped."
echo "Sample names are: ${WT1[@]}"
echo "There are $WT2_num WT2_num samples that will be mapped."
echo "Sample names are: ${WT2[@]}"
echo "There are $WT3_num WT3_num samples that will be mapped."
echo "Sample names are: ${WT3[@]}"
#
for (( w=0;w<${quantity};w +=2 )); do
	SAMPLE="${ALL_ARRAYS[$w]}"
	SAMPLE2="${ALL_ARRAYS[$w+1]}"
	SAMPLENAME=${SAMPLE%"_flexbar_"*}
	if [[ -f $SAMPLE ]]; then
		date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
		echo "Now I am BBmap mapping PE reads $SAMPLE and $SAMPLE2 - alignment"
		# To map quickly with very high precision and lower sensitivity, as when removing contaminant reads 
		# specific to a genome without risking false-positives with unmapped reads in fastaq and statistics:
		# To map vertebrate RNA-seq reads to a genome:
		# bbmap.sh in=reads.fq out=mapped.sam maxindel=200k ambig=random intronlen=20 xstag=us
		$BBMAP -Xmx"$RAM"g -threads="$THREADS" maxindel=200k ambig=random intronlen=20 xstag=us \
		in1="$SAMPLE" in2="$SAMPLE2" ref="$GENOME_NAME" out=$SAMPLENAME"_bbmap.bam" \
		outm1=$SAMPLENAME"_mapped_only_1.fq" outm2=$SAMPLENAME"_mapped_only_2.fq" \
		outu1=$SAMPLENAME"_unmapped_only_1.fq" outu2=$SAMPLENAME"_unmapped_only_2.fq" \
		covstats=$SAMPLENAME"_covstats.txt" covhist=$SAMPLENAME"_covhist.txt" \
		basecov=$SAMPLENAME"_basecov.txt" bincov=$SAMPLENAME"_bincov.txt"
		date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
		echo "Done BBmap aligning PE reads for sample $SAMPLENAME"
	else
		date +"%d/%m/%Y %H:%M:%S"
		echo "There is no $SAMPLE file for mapping!" && exit 1
	fi
done
wait
for a in *.bam
do
	FILE=$a
	FILENAME=${a%.*}
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Sorting and indexing BAM file $FILE"
	$SAMTOOLS sort -@ $THREADS --output-fmt BAM -o "$FILENAME.sorted.bam" "$FILE"
    $SAMTOOLS index -@ $THREADS "$FILENAME.sorted.bam" "$FILENAME.sorted.bai"
done
#wait
wait
##
#
mkdir -p $OUTPUT_DIR/{multiqc_log,bbstats,alignment,mapped_reads,unmapped_reads}
#
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Start creating MultiQC report"
$MULTIQC -o $OUTPUT_DIR/multiqc_log $OUTPUT_DIR/
#
echo "Removing genome and original data"
rm -v *.fa.gz
rm -v *.fastq.gz
echo "Finnished removing genome and original data"
#
for a in *.fq
do
	FILE=$a
	FILENAME=${a%.*}
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Start packing read library $FILE"
	gzip -v $FILE & # unpigz is not working on Krtecek server :-/
done
wait
#
echo "Moving mapping statistics to bbstats folder"
mv -v $OUTPUT_DIR/*.txt $OUTPUT_DIR/bbstats/
echo "Moving aligned reads to alignment folder"
mv -v $OUTPUT_DIR/*.bam $OUTPUT_DIR/alignment/
mv -v $OUTPUT_DIR/*.bai $OUTPUT_DIR/alignment/
echo "Moving mapped reads to mapped_reads folder"
mv -v $OUTPUT_DIR/*_mapped_only*.fq.gz $OUTPUT_DIR/mapped_reads/
echo "Moving unmapped reads to unmapped_reads folder"
mv -v $OUTPUT_DIR/*_unmapped_only*.fq.gz $OUTPUT_DIR/unmapped_reads/
##
############################################################################################
### Copy data from scratch back to home dir and clean scratch
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
