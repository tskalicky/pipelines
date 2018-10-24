#!/bin/bash
#PBS -l select=1:ncpus=30:mem=35gb:scratch_local=150gb
#PBS -l walltime=96:00:00
#PBS -q uv@wagap-pro.cerit-sc.cz
#PBS -j oe
#PBS -N 02_preprocess_SE_METTL16
#
# Preprocessing - RNA-Seq SE
#
# Requires FastQC, Trimmomatic (java), flexbar
#
## initialize the required application
module add fastQC-0.11.5
module add flexbar-3.0.3
module add trimmomatic-0.36
### Variables
INPUT_DIR=/storage/brno9-ceitec/home/tskalicky/FTO_project/raw_reads/paired_libs
OUTPUT_DIR=/storage/brno9-ceitec/home/tskalicky/FTO_project/preprocessed

ADAPTERS=/storage/brno9-ceitec/home/tskalicky/tools/adapters_v2.fasta
APPENDIX=".fastq.gz"
THREADS=$PBS_NUM_PPN # Number of threads to use

FASTQC=$(which fastqc)
FLEXBAR=$(which flexbar)
TRIMMOMATIC="java -jar /software/trimmomatic/0.36/dist/jar/trimmomatic-0.36.jar"
# Check the tools versions
echo $FASTQC
echo $FLEXBAR
echo $TRIMMOMATIC
############################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
cp -avr $INPUT_DIR $SCRATCHDIR
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR

############################################################################################
# commands
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
mkdir -p $SCRATCHDIR/{data/trimmed,fastqc/trimmed}  #creates whole subdirectory tree using -p and {}
QCDIR=$SCRATCHDIR/fastqc/trimmed/
TRIMDIR=$SCRATCHDIR/data/trimmed/
echo $QCDIR
echo $TRIMDIR
#trim reads witch flexbar
#usage of $FLEXBAR and $FASTQC variables in function is not working in PBS "Illegal instruction" error
#loop over every file and process them in background simultaniously
for i in *FTO_KO*$APPENDIX
do
	READ_FOR=$i
	SAMPLENAME=${i%.*.*}
	echo "Now I am processing SE reads $READ_FOR - Flexbar Trimming"
	
	flexbar -t $TRIMDIR$SAMPLENAME"_flexbar" -r $READ_FOR -a $ADAPTERS -ac -ao 3 -u 2 -m 17 -n 6 -z GZ -l MOD & 
	
	echo "Done processing SE reads $READ_FOR - Flexbar Trimming"

done

# Wait for all jobs to finish before exiting the job submission script
wait

cd $TRIMDIR

for j in *flexbar.fastq.gz 
do
	READ_FOR2=$j
	SAMPLENAME2=${j%.*.*}
	echo "Now I am processing SE reads $READ_FOR2 - Flexbar LEFT Trimming"

	flexbar -t $SAMPLENAME2"2" -r READ_FOR2 -a $ADAPTERS -ac -ae LEFT -ao 3 -u 2 -m 17 -n 6 -z GZ -l MOD &

	echo "Done processing SE reads $READ_FOR2 - Flexbar LEFT Trimming"

done

# Wait for all jobs to finish before exiting the job submission script
wait

fastqc --outdir $QCDIR --format fastq --threads $THREADS $SAMPLENAME_flexbar*

cd $QCDIR

echo "Number of reads after preprocessing is in fastqc/trimmed/total_sequences.txt"
for a in *zip
do
	unzip -q $a
	grep "Total Sequences" ${a%.zip*}/*.txt
done > total_sequences.txt
rm -r *fastqc

############################################################################################
### Copy data from scratch back to home dir and clean scratch
cp -avr $SCRATCHDIR $OUTPUT_DIR || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"