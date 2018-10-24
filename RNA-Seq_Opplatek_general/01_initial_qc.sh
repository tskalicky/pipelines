#!/bin/bash
#PBS -l select=1:ncpus=4:mem=8gb:scratch_local=200gb
#PBS -l walltime=04:00:00
#PBS -N 01_qc_initial_run
#
# Primary QC
#
# Requires FastQC, Kraken (minion, swan), multiQC
#

############################################################################################
### Variables
INPUT_DIR=/storage/brno3-cerit/home/opplatek/biocore/tichy_rna-seq_workshop/data
OUTPUT_DIR=/storage/brno3-cerit/home/opplatek/biocore/tichy_rna-seq_workshop
OUTPUT_DIR=${OUTPUT_DIR}/results/qc

APPENDIX="_001.fastq.gz"

ADAPTERS=/storage/brno2/home/opplatek/adapters_merge.txt # Adapters to scan for in fasta format

THREADS=$PBS_NUM_PPN

module add fastQC-0.11.5
FASTQC=$(which fastqc)
MINION=/storage/brno2/home/opplatek/tools/kraken/minion
SWAN=/storage/brno2/home/opplatek/tools/kraken/swan
module add python27-modules-gcc
#PYTHONPATH=$PYTHONPATH:/storage/brno2/home/opplatek/tools/MultiQC-1.0/lib/python2.7/site-packages/
#MULTIQC=/storage/brno2/home/opplatek/tools/MultiQC-1.0/bin/multiqc
MULTIQC=$(which multiqc)

# Check the tools versions
echo $FASTQC
echo $MINION
echo $SWAN
echo $MULTIQC
############################################################################################
### Copy inputs
cp $INPUT_DIR/*$APPENDIX $SCRATCH/
cp $ADAPTERS $SCRATCH/

cd $SCRATCH/

ADAPTERS=$(basename $ADAPTERS)

############################################################################################
### FastQC QC
### Minion adapter scan
mkdir -p $SCRATCH/minion

for i in *$APPENDIX
do
	$MINION search-adapter -i $i -show 3 -write-fasta $SCRATCH/minion/${i%.*}.minion.fasta
	$SWAN -r $SCRATCH/$ADAPTERS -q $SCRATCH/minion/${i%.*}.minion.fasta > $SCRATCH/minion/${i%.*}.minion.compare
done &

mkdir -p $SCRATCH/fastqc/raw

$FASTQC -t $THREADS -f fastq -o $SCRATCH/fastqc/raw *$APPENDIX
$MULTIQC -o $SCRATCH/fastqc/raw $SCRATCH/fastqc/raw/

cd $SCRATCH/fastqc/raw

echo "Number of reads before preprocessing is in fastqc/raw/total_sequences.txt"
for a in *zip
do
	unzip -q $a
	grep "Total Sequences" ${a%.zip*}/*.txt
done > total_sequences.txt
rm -r *fastqc

wait

# Get possible adapter sequences for all files - no checking, just taking 3rd row and/or match 
# HIGHLY EXPERIMENTAL AND NEEDS A MANUAL CHECK
for i in $SCRATCH/minion/*.compare
do
	echo $(basename $i .fastq.minion.compare) >> $SCRATCH/minion/adapters.txt
	grep -B1 -A1 "||||||||||||" $i | sed 's/^ //g' >> $SCRATCH/minion/adapters.txt
	echo -ne ">" >> $SCRATCH/minion/adapters_seq.txt
	echo $(basename $i .fastq.minion.compare) >> $SCRATCH/minion/adapters_seq.txt
	grep -A1 "||||||||||||" $i | sed 's/^ //g' | awk '{print $1}' | sed '/^||||||||||||/d'>> $SCRATCH/minion/adapters_seq.txt # Very rough extraction!
# Unsusable
#	echo -ne ">" >> $SCRATCH/minion/adapters.fa
#	echo ${i%.fastq*} >> $SCRATCH/minion/adapters.fa
#	cat $i | sed -n 3p | awk '{print $1}' >> $SCRATCH/minion/adapters.fa
done

############################################################################################
### Copy results
mkdir -p $OUTPUT_DIR

cp -r $SCRATCH/fastqc $OUTPUT_DIR/
cp -r $SCRATCH/minion $OUTPUT_DIR/

rm -r $SCRATCH/*
