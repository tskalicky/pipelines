#!/bin/bash
#PBS -l walltime=96:0:0
#PBS -q default@wagap-pro.cerit-sc.cz
#PBS -l select=1:ncpus=12:mem=50gb:scratch_local=250gb:os=debian9
#PBS -j oe
#PBS -N 02v4_fastQC_PE_FTO
#
# Preprocessing - RNA-Seq FastQC
#
# Requires FastQC
module add fastQC-0.11.5
### Variables
set -e
INPUT_DIR="/storage/brno9-ceitec/home/tskalicky/FTO_project/preprocessed/job_477994_wagap/data/trimmed"
OUTPUT_DIR="/storage/brno9-ceitec/home/tskalicky/FTO_project/preprocessed"
#
ADAPTERS="/storage/brno9-ceitec/home/tskalicky/tools/adapters_v3.fasta"
APPENDIX=".fastq.gz"
THREADS=$PBS_NUM_PPN # Number of threads to use
#
FASTQC=$(which fastqc)
# Check the tools versions
echo "INPUT_DIR path is" $INPUT_DIR > $INPUT_DIR/Pipeline_debuging.txt # debugging will be removed
echo "OUTPUT_DIR path is" $OUTPUT_DIR >> $INPUT_DIR/Pipeline_debuging.txt # debugging will be removed
echo "fastQC are located in" $FASTQC >> $INPUT_DIR/Pipeline_debuging.txt # debugging will be removed
cd $INPUT_DIR
echo "Entering INPUT_DIR" >> $INPUT_DIR/Pipeline_debuging.txt # debugging will be removed
############################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data are successfully copied out) or is killed
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
echo "SCRATCHDIR path is" $SCRATCHDIR >> $SCRATCHDIR/Pipeline_debuging.txt # debugging will be removed
echo "Start copying data to scratch" >> $SCRATCHDIR/Pipeline_debuging.txt # debugging will be removed
cp -v $INPUT_DIR/*.gz $SCRATCHDIR
cd $SCRATCHDIR
echo "Data for fastqc were copied to scratch" >> $SCRATCHDIR/Pipeline_debuging.txt # debugging will be removed

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created

############################################################################################
# commands
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME" >> $SCRATCHDIR/Pipeline_debuging.txt # debugging will be removed
mkdir $INPUT_DIR/secondary_fastqc  #creates whole subdirectory tree using -p and {}
QCDIR=$INPUT_DIR/secondary_fastqc/
#
echo "QCDIR path is:" $QCDIR
echo "QCDIR path is:" $QCDIR >> $INPUT_DIR/Pipeline_debuging.txt # debugging will be removed
#
echo "Now I am processing PE reads - Fastqc check" >> $INPUT_DIR/Pipeline_debuging.txt # debugging will be removed
fastqc --outdir $QCDIR --format fastq --threads $THREADS SRR3290147-TREX_WT_2_pass_LEFT_flexbar_2.fastq.gz \
SRR3290147-TREX_WT_2_pass_LEFT_flexbar_1.fastq.gz \
SRR3290145-FTO_KO_3_pass_LEFT_flexbar_2.fastq.gz \
SRR3290144-FTO_KO_2_pass_LEFT_flexbar_2.fastq.gz \
SRR3290143-FTO_KO_1_pass_LEFT_flexbar_2.fastq.gz \
SRR3290145-FTO_KO_3_pass_LEFT_flexbar_1.fastq.gz \
SRR3290144-FTO_KO_2_pass_LEFT_flexbar_1.fastq.gz \
SRR3290143-FTO_KO_1_pass_LEFT_flexbar_1.fastq.gz \
SRR3290148-TREX_WT_3_pass_LEFT_flexbar_2.fastq.gz \
SRR3290146-TREX_WT_1_pass_LEFT_flexbar_2.fastq.gz \
SRR3290148-TREX_WT_3_pass_LEFT_flexbar_1.fastq.gz \
SRR3290146-TREX_WT_1_pass_LEFT_flexbar_1.fastq.gz
echo "Done processing PE reads - Fastqc check"

cd $QCDIR
echo "entering" $QCDIR
echo "entering" $QCDIR >> $INPUT_DIR/Pipeline_debuging.txt # debugging will be removed

echo "Number of reads after preprocessing is in fastqc/trimmed/total_sequences.txt"
echo "Number of reads after preprocessing is in fastqc/trimmed/total_sequences.txt" >> $INPUT_DIR/Pipeline_debuging.txt # debugging will be removed
for a in *zip
do
	unzip -q $a
	grep "Total Sequences" ${a%.zip*}/*.txt
done > total_sequences.txt
rm -r ./*fastqc/

############################################################################################
#### Copy data from scratch back to home dir and clean scratch
#cp -avr $SCRATCHDIR $OUTPUT_DIR || export CLEAN_SCRATCH=false
#echo "Script finished on:"
#date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Script finished on:" >> $INPUT_DIR/Pipeline_debuging.txt # debugging will be removed
date +"%d/%m/%Y %H:%M:%S $HOSTNAME" >> $INPUT_DIR/Pipeline_debuging.txt # debugging will be removed