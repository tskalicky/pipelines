#!/bin/bash
#PBS -l walltime=96:0:0
#PBS -q uv@wagap-pro.cerit-sc.cz
#PBS -l select=1:ncpus=6:mem=40gb:scratch_local=50gb
#PBS -j oe
#PBS -N 02v43_fastQC_SE_FTO_Input+CLIP
#
# Preprocessing - RNA-Seq FastQC
#
# Requires FastQC
module add fastQC-0.11.5
### Variables
INPUT_DIR="/storage/brno3-cerit/home/tskalicky/FTO/preprocessed/job_537816_wagap/data/trimmed"
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/FTO/preprocessed/job_537816_wagap/fastqc/trimmed"
#
ADAPTERS="/storage/brno9-ceitec/home/tskalicky/tools/adapters_v4.fa"
APPENDIX=".fastq.gz"
THREADS=$PBS_NUM_PPN # Number of threads to use
#
FASTQC=$(which fastqc)
# Check the tools versions
echo "INPUT_DIR path is" $INPUT_DIR 
echo "OUTPUT_DIR path is" $OUTPUT_DIR
echo "fastQC is located in" $FASTQC
cd $INPUT_DIR
echo "Entering INPUT_DIR"
############################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data are successfully copied out) or is killed
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
echo "SCRATCHDIR path is" $SCRATCHDIR
echo "Start copying data to scratch"
cp -v $INPUT_DIR/*.gz $SCRATCHDIR
cd $SCRATCHDIR
ls -c1 $SCRATCHDIR #print filenames of copied files
echo "Data for fastqc were copied to scratch"

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created

############################################################################################
# commands
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
mkdir $SCRATCHDIR/secondary_fastqc  #creates whole subdirectory tree using -p and {}
QCDIR=$SCRATCHDIR/secondary_fastqc/
#
echo "QCDIR path is:" $QCDIR

#
echo "Now I am processing SE reads - Fastqc check"
fastqc --outdir $QCDIR --format fastq --threads $THREADS SRR3290136-FTO_CLIP1_pass_1_flexbar.fastq.gz \
SRR3290137-FTO_CLIP2_pass_1_flexbar.fastq.gz \
SRR3290138-FTO_CLIP3_pass_1_flexbar.fastq.gz \
SRR3290139-FTO_Input1_pass_1_flexbar.fastq.gz \
SRR3290140-FTO_Input2_pass_1_flexbar.fastq.gz \
SRR3290142-FTO_Input3_pass_1_flexbar.fastq.gz
echo "Done processing PE reads - Fastqc check"

cd $QCDIR
echo "entering" $QCDIR
echo "Number of reads after preprocessing is in fastqc/trimmed/total_sequences.txt"
#
for a in *zip
do
	unzip -q $a
	grep "Total Sequences" ${a%.zip*}/*.txt
done > total_sequences.txt
rm -r ./*fastqc/

############################################################################################
#### Copy data from scratch back to home dir and clean scratch
cp -avr $SCRATCHDIR $OUTPUT_DIR || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
