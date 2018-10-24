#!/bin/bash
#
# Preprocessing - Lexogen QuanSeq SE (RNA-Seq)
#
# Should handle both Lexogen QuantSeq FWD (TRIM_LE="true") and REV (TRIM_LE="false")
#
# Requires FastQC, BBMap (java), multiQC
#
# There is a difference between preprocessing of Lexogen QuantSeq FWD and REV! Briefly, FWD does need trimming of first bases and REV does not.
# For more information see https://www.lexogen.com/quantseq-data-analysis/ for FWD and https://www.lexogen.com/quantseq-data-analysis-rev/ for REV
#
# IN PROGRESS
#
############################################################################################
### Variables
INPUT_DIR=/storage/brno3-cerit/home/opplatek/biocore/karin_quantSeq/data/raw
OUTPUT_DIR=/storage/brno3-cerit/home/opplatek/biocore/karin_quantSeq/data/preprocessed
OUTPUT_DIR_QC=/storage/brno3-cerit/home/opplatek/biocore/karin_quantSeq/data/raw/results/qc

TRIM_AD="true" # Trim also adapters? true or false
TRIM_LE="true" # Force trim 5' end? true (QuantSeq FWD) or false (QuantSeq REV)
TRIM_RI="false" # Force trim also 3' end?
ADAPTERS=/storage/brno2/home/opplatek/tools/bbmap-37.25/resources/truseq.fa.gz # Adapters to trim, required if $TRIM_AD is true

RD_LENGTH=65 # Read length from the sequencing OR trimming down to this length. Illumina sometimes reports N+1 http://seqanswers.com/forums/archive/index.php/t-31154.html

TRIM_LEFT=10 # Applied only if trim lefth is true; adapter should be 12 bp long but from the nuc.distribution it seems trimming 11 bp might be enough and at web they recommend to trim 10 from R1 and 7 from R2...
TRIM_RIGHT=3 # Applied only if trim right is true

APPENDIX=".fastq.gz"

PHRED_TRIM=5 # Trim the 3' end of read if four consequent bases have average base quality smaller than this value
LEN_FILTER=35 # Filter sequences shorter than this value after quality trimming

MIN_OVERLAP=3 # Minimal overlap for matching adapter

THREADS=$PBS_NUM_PPN # Number of threads to use

module add fastQC-0.11.5
FASTQC=fastqc
module add jdk-8
BBMAP=/storage/brno2/home/opplatek/tools/bbmap-37.25
module add python27-modules-gcc
PYTHONPATH=$PYTHONPATH:/storage/brno2/home/opplatek/tools/MultiQC-1.0/lib/python2.7/site-packages/
MULTIQC=/storage/brno2/home/opplatek/tools/MultiQC-1.0/bin/multiqc

############################################################################################
### Copy inputs
cp $INPUT_DIR/*$APPENDIX $SCRATCH/

cd $SCRATCH/

############################################################################################
### Preparation and extra flags for BBDuk
msg="Trimming quality" # Default
extra_flags=""
if [ "$TRIM_AD" == "true" ]; then
	cp $ADAPTERS $SCRATCH/
	ADAPTERS=$(basename $ADAPTERS)

	msg="$msg, adapter"
	extra_flags="ktrim=r mink=$MIN_OVERLAP ref=$SCRATCH/$ADAPTERS"
else
	msg="$msg, no adapter"
fi

if [ "$TRIM_LE" == "true" ]; then
	msg="$msg, $TRIM_LEFT left"
	extra_flags="$extra_flags forcetrimleft=$TRIM_LEFT"
else
	msg="$msg, nothing on left"
fi

if [ "$TRIM_RI" == "true" ]; then
	msg="$msg, $TRIM_RIGHT right"
	extra_flags="$extra_flags forcetrimright2=$TRIM_RIGHT"
else
	msg="$msg, nothing on right"
fi

echo "$msg"

RD_LENGTH=$[$RD_LENGTH -1] # BBDuk is 0-based 
############################################################################################
### Preprocess and qc
mkdir $SCRATCH/bbduk

for i in *$APPENDIX
do
	READ_FOR=$i

	echo "Now I am processing SE reads $READ_FOR - Trimming"

	unpigz -c -p $PBS_NUM_PPN ${READ_FOR} | $BBMAP/bbduk.sh in=stdin.fq out=${READ_FOR%.fastq*}_trim.fastq.gz threads=$PBS_NUM_PPN forcetrimright=$RD_LENGTH qtrim=rl trimq=$QT minlength=$LEN_FILTER useshortkmers=t literal=GGGGGGGGGG,AAAAAAAAAA k=10 $extra_flags &>$SCRATCH/bbduk/${READ_FOR%.fastq*}.bbduk.out # useshortkmers=t - probably old version of bbduk but doesn't give error and is in the code; sliding window can be enabled by "qtrim=window,4 trimq=$QT" but it is not recommended by Brian Bushnell http://seqanswers.com/forums/showthread.php?t=42776&page=7

	echo "Done processing SE reads $READ_FOR - Trimming"

done

mkdir -p $SCRATCH/fastqc/post

$FASTQC --outdir $SCRATCH/fastqc/post --format fastq --threads $THREADS $SCRATCH/*_trim.fastq.gz

cd $SCRATCH/fastqc/post

echo "Number of reads after preprocessing is in fastqc/post/total_sequences.txt"
for a in *zip
do
	unzip -q $a
	grep "Total Sequences" ${a%.zip*}/*.txt
done > total_sequences.txt
rm -r *fastqc

$MULTIQC -o $SCRATCH/fastqc/post $SCRATCH/fastqc/post/

############################################################################################
### Copy results
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR_QC/fastqc

cp -r $SCRATCH/fastqc/post $OUTPUT_DIR_QC/fastqc/ &
cp -r $SCRATCH/bbduk $OUTPUT_DIR_QC/ &
cp $SCRATCH/*_trim.fastq.gz $OUTPUT_DIR/

rm -r $SCRATCH/*
