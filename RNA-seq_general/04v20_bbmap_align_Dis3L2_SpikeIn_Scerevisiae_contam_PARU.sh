#!/bin/bash
#PBS -l nodes=1:ppn=30
#PBS -l mem=80gb
#PBS -j oe
#PBS -N bbmap_Dis3L2_yeast_contam
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
INPUT_DIR="/home/users/tskalicky/CEITEC/Dis3L2/Dasa_SpikeIn/data/trimmed"
MY_GENOME="/home/users/tskalicky/CEITEC/genomes/Saccharomyces_cerevisiae/Scerevisiae_R64_genomic.fa"
MY_GTF="/home/users/tskalicky/CEITEC/genomes/Saccharomyces_cerevisiae/Scerevisiae_R64_genomic.gff.gz"
OUTPUT_DIR="/home/users/tskalicky/CEITEC/Dis3L2/Dasa_SpikeIn/mapping/Scerevisiae_genome"
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
MY_RAM="80"
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
find "$INPUT_DIR" -name "*.fastq.gz" -exec cp -vt "$OUTPUT_DIR" {} +
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
set -f    # disable globbing, needed for filling arrays with output from find
# declare an array variables
declare -a OAT_CYTO=($(find . -name '*OAT_cyto*.fastq.gz' -exec basename {} \; | sort -n))
declare -a OAT_NONFRAC=($(find . -name '*OAT_nonfrac*.fastq.gz' -exec basename {} \; | sort -n))
declare -a OAT_NUCL=($(find . -name '*OAT_nucl*.fastq.gz' -exec basename {} \; | sort -n))
declare -a U12_CYTO=($(find . -name '*U12_cyto*.fastq.gz' -exec basename {} \; | sort -n))
declare -a U12_NONFRAC=($(find . -name '*U12_nonfrac*.fastq.gz' -exec basename {} \; | sort -n))
declare -a U12_NUCL=($(find . -name '*U12_nucl*.fastq.gz' -exec basename {} \; | sort -n))
declare -a ARRAY_NAMES=("OAT_CYTO" "OAT_NONFRAC" "OAT_NUCL" "U12_CYTO" "U12_NONFRAC" "U12_NUCL")
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
####################################################################################################
### Genome and annotation preparation
# GENOME_NAME=$(basename $MY_GENOME) # for normal alignment, not for RSEM analysis prep
# GTF_NAME=$(basename $RSEM_GTF) # for normal alignment, not for RSEM analysis prep
GENOME_NAME=$(basename $MY_GENOME)
GTF_NAME=$(basename $MY_GTF)

# unpigz -p $THREADS $GENOME_NAME # for normal alignment, not for RSEM analysis prep
# unpigz -p $THREADS $GTF_NAME # for normal alignment, not for RSEM analysis prep

# GENOME=${GENOME_NAME%.*} # for normal alignment, not for RSEM analysis prep
# GTF=${GTF_NAME%.*} # for normal alignment, not for RSEM analysis prep

GENOME=${GENOME_NAME%.*}
GTF=${GTF_NAME%.*}
#
# date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# echo "Decompressing files:"
# $UNPIGZ -v -p $THREADS -d $GENOME_NAME
# unpigz -v -p $THREADS -d $GTF_NAME
#
### Commands 
#
# Genome index preparation
# $BBMAP -threads="$THREADS" ref="$GENOME_NAME"
#
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
for (( w=0;w<${quantity};w +=2 ));
do
	SAMPLE="${ALL_ARRAYS[$w]}"
	SAMPLE2="${ALL_ARRAYS[$w+1]}"
	SAMPLENAME=${SAMPLE%"_flexbar_1.fastq.gz"*}
	if [[ -f $SAMPLE ]]; then
		date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
		echo "Now I am BBmap mapping PE reads $SAMPLE and $SAMPLE2 - alignment"
		# To map quickly with very high precision and lower sensitivity, as when removing contaminant reads 
		# specific to a genome without risking false-positives with unmapped reads in fastaq and statistics:
		# To map vertebrate RNA-seq reads to a genome:
		# bbmap.sh in=reads.fq out=mapped.sam maxindel=200k ambig=random intronlen=20 xstag=us
		$BBMAP -Xmx"$RAM"g -threads="$THREADS" minratio=0.9 maxindel=3 bwr=0.16 bw=12 fast minhits=2 qtrim=r trimq=10 \
		untrim idtag printunmappedcount kfilter=25 maxsites=1 k=14 \
		in1="$SAMPLE" in2="$SAMPLE2" ref="$GENOME_NAME" nodisk out=$SAMPLENAME"_bbmap.bam" \
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
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Sorting and indexing BAM files"
for a in *.bam;
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
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Start creating MultiQC report"
mkdir -p $OUTPUT_DIR/{multiqc_log,bbstats,alignment,mapped_reads,unmapped_reads}
#
$MULTIQC -o $OUTPUT_DIR/multiqc_log $OUTPUT_DIR/
#
echo "Moving mapping statistics to bbstats folder"
mv -v $OUTPUT_DIR/*.txt $OUTPUT_DIR/bbstats/
echo "Moving aligned reads to alignment folder"
mv -v $OUTPUT_DIR/*.bam $OUTPUT_DIR/alignment/
mv -v $OUTPUT_DIR/*.bai $OUTPUT_DIR/alignment/
echo "Moving mapped reads to mapped_reads folder"
mv -v $OUTPUT_DIR/*_mapped_only*.fq $OUTPUT_DIR/mapped_reads/
echo "Moving unmapped reads to unmapped_reads folder"
mv -v $OUTPUT_DIR/*_unmapped_only*.fq $OUTPUT_DIR/unmapped_reads/
##
############################################################################################
### Copy data from scratch back to home dir and clean scratch
rm -v $OUTPUT_DIR/*.fastq.gz
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
