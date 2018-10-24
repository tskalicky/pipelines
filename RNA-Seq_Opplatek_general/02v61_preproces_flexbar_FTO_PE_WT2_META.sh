#!/bin/bash
#PBS -l walltime=336:0:0 
#PBS -q default 
#PBS -l select=1:ncpus=5:mem=50gb:scratch_local=200gb
#PBS -j oe
#PBS -N 02v61_PE_FTO_WT2
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
INPUT_DIR="storage-brno9-ceitec.metacentrum.cz:~/FTO_project/raw_reads/paired_libs"
OUTPUT_DIR="storage-brno9-ceitec.metacentrum.cz:~/FTO_project/preprocessed"

ADAPTERS="storage-brno9-ceitec.metacentrum.cz:~/tools/adapters_v4.fa"
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
#using SCRATCHDIR storage using scp (disks not connected via NFS)
# clean the SCRATCH when job finishes (and data are successfully copied out) or is killed
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
# use "scp -r ..." in case of copying directories
scp $INPUT_DIR/SRR3290147-TREX_WT_2_pass_1.fastq.gz $INPUT_DIR/SRR3290147-TREX_WT_2_pass_2.fastq.gz $ADAPTERS $SCRATCHDIR
cd $SCRATCHDIR

echo "Following reads were copyed to scratch:"
ls -c1
echo "Following reads were copyed to scratch:" >> $SCRATCHDIR/Pipeline_debuging.txt # debugging will be removed
ls -c1 >> $SCRATCHDIR/Pipeline_debuging.txt # debugging will be removed

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
echo "QCDIR path is" $QCDIR >> $SCRATCHDIR/Pipeline_debuging.txt # debugging will be removed
echo "TRIMDIR path is" $TRIMDIR >> $SCRATCHDIR/Pipeline_debuging.txt # debugging will be removed

#usage of $FLEXBAR and $FASTQC variables in function is not working in PBS "Illegal instruction" error
#loop over every file and process them in background simultaniously
#FLEXBAR is using only like 1/2 the cores assigned
echo "Now I am processing FTO PE reads - Flexbar ANY Trimming"

#flexbar -t $TRIMDIR"SRR3290148-TREX_WT_3_flexbar_qual" -r SRR3290148-TREX_WT_3_pass_1.fastq.gz -p SRR3290148-TREX_WT_3_pass_2.fastq.gz -q TAIL -qf i1.8 -qt 30 -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 5 -z GZ
flexbar -t $TRIMDIR"SRR3290147-TREX_WT_2_flexbar_qual" -r SRR3290147-TREX_WT_2_pass_1.fastq.gz -p SRR3290147-TREX_WT_2_pass_2.fastq.gz -q TAIL -qf i1.8 -qt 30 -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 5 -z GZ -l MOD
#flexbar -t $TRIMDIR"SRR3290146-TREX_WT_1_flexbar_qual" -r SRR3290146-TREX_WT_1_pass_1.fastq.gz -p SRR3290146-TREX_WT_1_pass_2.fastq.gz -q TAIL -qf i1.8 -qt 30 -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 5 -z GZ -l MOD
#flexbar -t $TRIMDIR"SRR3290145-FTO_KO_3_flexbar_qual" -r SRR3290145-FTO_KO_3_pass_1.fastq.gz -p SRR3290145-FTO_KO_3_pass_2.fastq.gz -q TAIL -qf i1.8 -qt 30 -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 5 -z GZ -l MOD
#flexbar -t $TRIMDIR"SRR3290144-FTO_KO_2_flexbar_qual" -r SRR3290144-FTO_KO_2_pass_1.fastq.gz -p SRR3290144-FTO_KO_2_pass_2.fastq.gz -q TAIL -qf i1.8 -qt 30 -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 5 -z GZ -l MOD
#flexbar -t $TRIMDIR"SRR3290143-FTO_KO_1_flexbar_qual" -r SRR3290143-FTO_KO_1_pass_1.fastq.gz -p SRR3290143-FTO_KO_1_pass_2.fastq.gz -q TAIL -qf i1.8 -qt 30 -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 5 -z GZ -l MOD

# Wait for all jobs to finish before continuing the script
#wait
echo "Done processing FTO PE reads - Flexbar ANY Trimming"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME" >> Pipeline_debuging.txt # debugging will be removed

cd $TRIMDIR

#bylo nutne upravit volani souboru pro fastqc
fastqc --outdir $QCDIR --format fastq --threads $THREADS *.fasq.gz
#
cd $QCDIR
echo "Number of reads after preprocessing is in fastqc/trimmed/total_sequences.txt"
for a in *zip
do
	unzip -q $a
	grep "Total Sequences" ${a%.zip*}/*.txt
done > total_sequences.txt
rm -r ./*fastqc/

############################################################################################
### Copy data from scratch back to home dir and clean scratch
cd $SCRATCHDIR
rm $SCRATCHDIR/*$APPENDIX
scp -r $SCRATCHDIR $OUTPUT_DIR || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
