#!/bin/bash
### Variables
INPUT_DIR=/Users/tomas/Data/METTL16/raw
OUTPUT_DIR=/Users/tomas/Data/METTL16/data/trimmed
OUTPUT_DIR_QC=//Users/tomas/Data/METTL16/fastqc/trimmed

ADAPTERS=/Users/tomas/Data/METTL16/adapters.fasta
APPENDIX=".fastq.gz"
THREADS=3 # Number of threads to use

FASTQC="/Users/tomas/bin/FastQC.app/Contents/MacOS/fastqc"
FLEXBAR=$(which flexbar)
#TRIMMOMATIC="java -jar /storage/brno9-ceitec/home/tskalicky/tools/Trimmomatic-0.36"
# Check the tools versions
echo $FASTQC
echo $FLEXBAR
#echo $TRIMMOMATIC


############################################################################################
cd $INPUT_DIR
# commands
date +"%m/%d/%Y %H:%M:%S $HOSTNAME"
mkdir -p /Users/tomas/Data/METTL16/fastqc/trimmed
mkdir -p /Users/tomas/Data/METTL16/data/trimmed
#trim reads witch flexbar
for i in *$APPENDIX
do
	READ_FOR=$i
	echo "Now I am processing SE reads $READ_FOR - Flexbar Trimming"
	
	$FLEXBAR -t $READ_FOR"_flexbar" -r $READ_FOR -a $ADAPTERS -ac -ao 3 -u 2 -m 17 -n $THREADS -z GZ 

	echo "Done processing SE reads $READ_FOR - Flexbar Trimming"

done

#$FASTQC --outdir $OUTPUT_DIR_QC --format fastq --threads $THREADS $INPUT_DIR/*_trim.fastq.gz
#
#cd $OUTPUT_DIR_QC
#
#echo "Number of reads after preprocessing is in fastqc/trimmed/total_sequences.txt"
#for a in *zip
#do
#	unzip -q $a
#	grep "Total Sequences" ${a%.zip*}/*.txt
#done > total_sequences.txt
##rm -r *fastqc
############################################################################################
#echo "Script finished on:"
#date +"%m/%d/%Y %H:%M:%S $HOSTNAME"






