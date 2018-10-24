#!/bin/bash
#PBS -l nodes=1:ppn=20
#PBS -l walltime=96:00:00
#PBS -l mem=80gb
#PBS -k oe
#PBS -N bbmap_Dis3L2_spikein_OAT_U12_only
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
INPUT_DIR="/home/users/tskalicky/CEITEC/Dis3L2/Dasa_SpikeIn/mapping/Scerevisiae_genome/unmapped_reads"
MY_GENOME="/home/users/tskalicky/CEITEC/genomes/human/OAT_U12.fa.gz"
MY_GTF="/home/users/tskalicky/CEITEC/genomes/human/OAT_U12_mapping.gtf.gz"
OUTPUT_DIR="/home/users/tskalicky/CEITEC/Dis3L2/Dasa_SpikeIn/mapping/OAT_U12"
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
echo $PATH | tr ":" "\n" | nl
#
# Set number of CPU
THREADS=$PBS_NUM_PPN
# Set MAX ram
MY_RAM="80" # SET THE AMOUNT OF RAM HERE !!!
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
find "$INPUT_DIR" -maxdepth 1 -name *.fq.gz -exec cp -avt "$OUTPUT_DIR" {} +
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
declare -a OAT_CYTO=($(find . -maxdepth 1 -name  '*OAT_cyto*.fq.gz' -exec basename {} \; | sort -n))
declare -a OAT_NONFRAC=($(find . -maxdepth 1 -name  '*OAT_nonfrac*.fq.gz' -exec basename {} \; | sort -n))
declare -a OAT_NUCL=($(find . -maxdepth 1 -name  '*OAT_nucl*.fq.gz' -exec basename {} \; | sort -n))
declare -a U12_CYTO=($(find . -maxdepth 1 -name  '*U12_cyto*.fq.gz' -exec basename {} \; | sort -n))
declare -a U12_NONFRAC=($(find . -maxdepth 1 -name  '*U12_nonfrac*.fq.gz' -exec basename {} \; | sort -n))
declare -a U12_NUCL=($(find . -maxdepth 1 -name  '*U12_nucl*.fq.gz' -exec basename {} \; | sort -n))
declare -a ARRAY_NAMES=("OAT_CYTO" "OAT_NONFRAC" "OAT_NUCL" "U12_CYTO" "U12_NONFRAC" "U12_NUCL")
declare -a ALL_ARRAYS=("${OAT_CYTO[@]}" "${OAT_NONFRAC[@]}" "${OAT_NUCL[@]}" "${U12_CYTO[@]}" "${U12_NONFRAC[@]}" "${U12_NUCL[@]}" )
# declare -a ARRAY_NAMES=("OAT_NONFRAC" "OAT_NUCL" "U12_CYTO" "U12_NONFRAC" "U12_NUCL")
# declare -a ALL_ARRAYS=("${OAT_NONFRAC[@]}" "${OAT_NUCL[@]}" "${U12_CYTO[@]}" "${U12_NONFRAC[@]}" "${U12_NUCL[@]}" )
set +f # enable globbing, because is needed for other parts of the script, like variable expansions !!!
# get length of an array
OAT_cyto_num="${#OAT_CYTO[@]}"
OAT_nonfrac_num="${#OAT_NONFRAC[@]}"
OAT_nucl_num="${#OAT_NUCL[@]}"
U12_cyto_num="${#U12_CYTO[@]}"
U12_nonfrac_num="${#U12_NONFRAC[@]}"
U12_nucl_num="${#U12_NUCL[@]}"
quantity="${#ALL_ARRAYS[@]}"
#
####################################################################################################
### Genome and annotation preparation
GENOME_NAME=$(basename $MY_GENOME)
GTF_NAME=$(basename $MY_GTF)

GENOME=${GENOME_NAME%.*}
GTF=${GTF_NAME%.*}
#
# date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# echo "Decompressing files:"
# $UNPIGZ -v -p $THREADS -d $GENOME_NAME
# unpigz -v -p $THREADS -d $GTF_NAME
# gzip -vd $GENOME_NAME
#
### Commands 
## BBmap - Alignment
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
#
for (( w=0;w<${quantity};w +=2 )); do
	SAMPLE="${ALL_ARRAYS[$w]}"
	SAMPLE2="${ALL_ARRAYS[$w+1]}"
	SAMPLENAME=${SAMPLE%"_unmapped_only_"*}
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
rm -v *.fa
rm -v *.fq.gz
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
