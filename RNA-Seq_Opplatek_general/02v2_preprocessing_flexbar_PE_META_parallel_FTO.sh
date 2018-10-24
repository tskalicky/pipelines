#!/bin/bash
#PBS -l select=1:ncpus=30:mem=30gb:scratch_local=250gb
#PBS -l walltime=96:00:00
#PBS -q uv@wagap-pro.cerit-sc.cz
#PBS -j oe
#PBS -N 02.1_trim_PE_FTO
#
# Preprocessing - RNA-Seq PE
#
# Requires FastQC, Trimmomatic (java), flexbar
#
## initialize the required application
module add fastQC-0.11.5
module add flexbar-3.0.3
module add trimmomatic-0.36
### Variables
INPUT_DIR=/mnt/storage-brno9-ceitec/home/tskalicky/FTO_project/preprocessed/job_471098.wagap-pro.cerit-sc.cz/data/trimmed
#INPUT_DIR=/storage/brno9-ceitec/home/tskalicky/FTO_project/raw_reads/paired_libs
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
cp -avr $INPUT_DIR/*.fastq.gz $SCRATCHDIR
cd $SCRATCHDIR
echo "Data were copied to scratch" > Pipeline_debuging.txt # debugging will be removed

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR

############################################################################################
# commands
echo "Starting trimming script on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
mkdir -p $SCRATCHDIR/{data/trimmed,fastqc/trimmed}  #creates whole subdirectory tree using -p and {}
QCDIR=$SCRATCHDIR/fastqc/trimmed/
TRIMDIR=$SCRATCHDIR/data/trimmed/
echo "fastqc directory is:" $QCDIR
echo "Trimmed reads dir is:" $TRIMDIR
#trim reads witch flexbar
#usage of $FLEXBAR and $FASTQC variables in function is not working in PBS "Illegal instruction" error
#How to do loops with paired reads? How to select correct pairs?
##First round of RIGHT trimming done - skipping
##echo "Now I am processing FTO PE reads - Flexbar Trimming"
#
#flexbar -t $TRIMDIR"SRR3290148-TREX_WT_3_pass_flexbar" -r SRR3290148-TREX_WT_3_pass_1.fastq.gz -p SRR3290148-TREX_WT_3_pass_2.fastq.gz -a $ADAPTERS -ac -ao 3 -u 2 -m 17 -n 10 -z GZ -l MOD & 
#flexbar -t $TRIMDIR"SRR3290147-TREX_WT_2_pass_flexbar" -r SRR3290147-TREX_WT_2_pass_1.fastq.gz -p SRR3290147-TREX_WT_2_pass_2.fastq.gz -a $ADAPTERS -ac -ao 3 -u 2 -m 17 -n 10 -z GZ -l MOD & 
#flexbar -t $TRIMDIR"SRR3290146-TREX_WT_1_pass_flexbar" -r SRR3290146-TREX_WT_1_pass_1.fastq.gz -p SRR3290146-TREX_WT_1_pass_2.fastq.gz -a $ADAPTERS -ac -ao 3 -u 2 -m 17 -n 10 -z GZ -l MOD & 
#flexbar -t $TRIMDIR"SRR3290145-FTO_KO_3_pass_flexbar" -r SRR3290145-FTO_KO_3_pass_1.fastq.gz -p SRR3290145-FTO_KO_3_pass_2.fastq.gz -a $ADAPTERS -ac -ao 3 -u 2 -m 17 -n 10 -z GZ -l MOD & 
#flexbar -t $TRIMDIR"SRR3290144-FTO_KO_2_pass_flexbar" -r SRR3290144-FTO_KO_2_pass_1.fastq.gz -p SRR3290144-FTO_KO_2_pass_2.fastq.gz -a $ADAPTERS -ac -ao 3 -u 2 -m 17 -n 10 -z GZ -l MOD & 
#flexbar -t $TRIMDIR"SRR3290143-FTO_KO_1_pass_flexbar" -r SRR3290143-FTO_KO_1_pass_1.fastq.gz -p SRR3290143-FTO_KO_1_pass_2.fastq.gz -a $ADAPTERS -ac -ao 3 -u 2 -m 17 -n 10 -z GZ -l MOD
#
## Wait for all jobs to finish before continuing the script
#wait
#echo "Done processing FTO PE reads - Flexbar Trimming"
#date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
#date +"%d/%m/%Y %H:%M:%S $HOSTNAME" >> Pipeline_debuging.txt # debugging will be removed
#
#cd $TRIMDIR

echo "Now I am processing FTO PE reads - Flexbar LEFT Trimming"

flexbar -t $TRIMDIR"SRR3290148-TREX_WT_3_pass_LEFT_flexbar" -r SRR3290148-TREX_WT_3_pass_flexbar_1.fastq.gz -p SRR3290148-TREX_WT_3_pass_flexbar_2.fastq.gz -a $ADAPTERS -ac -ae LEFT -ao 3 -u 2 -m 17 -n 10 -z GZ -l MOD & 
flexbar -t $TRIMDIR"SRR3290147-TREX_WT_2_pass_LEFT_flexbar" -r SRR3290147-TREX_WT_2_pass_flexbar_1.fastq.gz -p SRR3290147-TREX_WT_2_pass_flexbar_2.fastq.gz -a $ADAPTERS -ac -ae LEFT -ao 3 -u 2 -m 17 -n 10 -z GZ -l MOD & 
flexbar -t $TRIMDIR"SRR3290146-TREX_WT_1_pass_LEFT_flexbar" -r SRR3290146-TREX_WT_1_pass_flexbar_1.fastq.gz -p SRR3290146-TREX_WT_1_pass_flexbar_2.fastq.gz -a $ADAPTERS -ac -ae LEFT -ao 3 -u 2 -m 17 -n 10 -z GZ -l MOD & 
flexbar -t $TRIMDIR"SRR3290145-FTO_KO_3_pass_LEFT_flexbar" -r SRR3290145-FTO_KO_3_pass_flexbar_1.fastq.gz -p SRR3290145-FTO_KO_3_pass_flexbar_2.fastq.gz -a $ADAPTERS -ac -ae LEFT -ao 3 -u 2 -m 17 -n 10 -z GZ -l MOD & 
flexbar -t $TRIMDIR"SRR3290144-FTO_KO_2_pass_LEFT_flexbar" -r SRR3290144-FTO_KO_2_pass_flexbar_1.fastq.gz -p SRR3290144-FTO_KO_2_pass_flexbar_2.fastq.gz -a $ADAPTERS -ac -ae LEFT -ao 3 -u 2 -m 17 -n 10 -z GZ -l MOD & 
flexbar -t $TRIMDIR"SRR3290143-FTO_KO_1_pass_LEFT_flexbar" -r SRR3290143-FTO_KO_1_pass_flexbar_1.fastq.gz -p SRR3290143-FTO_KO_1_pass_flexbar_2.fastq.gz -a $ADAPTERS -ac -ae LEFT -ao 3 -u 2 -m 17 -n 10 -z GZ -l MOD

# Wait for all jobs to finish before continuing the script
wait
echo "Done processing FTO PE reads - Flexbar LEFT Trimming"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
#bylo nutne upravit volani souboru pro fastqc
fastqc --outdir $QCDIR --format fastq --threads $THREADS *pass_LEFT_flexbar*.fasq.gz

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
cd $SCRATCHDIR
tar -cvf 02_trim_PE_FTO.tar $SCRATCHDIR --remove-files
cp -avr $SCRATCHDIR $OUTPUT_DIR || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
