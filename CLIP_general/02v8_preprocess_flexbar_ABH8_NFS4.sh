#!/bin/bash
#PBS -l walltime=96:0:0 
#PBS -q uv@wagap-pro.cerit-sc.cz
#PBS -l select=1:ncpus=18:mem=50gb:scratch_local=100gb
#PBS -j oe
#PBS -N 02v8_flexbar_fastqc_ABH8
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


INPUT_ABH8_CLIP="/storage/brno3-cerit/home/tskalicky/ABH8/CLIP/preprocessed/ABH8_*_flexbar3.fastq.gz"
INPUT_ABH8_RNA="/storage/brno3-cerit/home/tskalicky/ABH8/RNA-seq/ABH8_*_clipped_qf.fastq.gz"
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/ABH8/preprocessed"

ADAPTERS="/storage/brno3-cerit/home/tskalicky/tools/adapters_v4.fa"
APPENDIX=".fastq.gz"
THREADS=$PBS_NUM_PPN # Number of threads to use

FASTQC=$(which fastqc)
FLEXBAR=$(which flexbar)
TRIMMOMATIC="java -jar /software/trimmomatic/0.36/dist/jar/trimmomatic-0.36.jar"
# Check the tools versions
echo $FASTQC
echo $FLEXBAR
echo $TRIMMOMATIC

####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
cd /storage/brno3-cerit/home/tskalicky/ABH8
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
cp -av $INPUT_ABH8_CLIP $INPUT_ABH8_RNA $ADAPTERS $SCRATCHDIR

cd $SCRATCHDIR

echo "Following files were copied to scratch:"
ls -c1

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR

############################################################################################
# commands
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
mkdir -p $SCRATCHDIR/{data/trimmed,fastqc/trimmed}  #creates whole subdirectory tree using -p and {}
#
QCDIR=$SCRATCHDIR/fastqc/trimmed/
TRIMDIR=$SCRATCHDIR/data/trimmed/
echo "QCdir is:" $QCDIR
echo "TRIMdir is:" $TRIMDIR

#usage of $FLEXBAR and $FASTQC variables in function is not working in PBS "Illegal instruction" error
#loop over every file and process them in background simultaniously
#FLEXBAR is using only like 1/2 the cores assigned
# POZOR! Flexbar nebezi na Metacentru pri pouziti uzlu s os=debian9 !
echo "Now I am processing ABH8 CLIP and RNAseq SE reads - Flexbar ANY and quality Trimming"
#
flexbar -t $TRIMDIR"ABH8_1_RNAseq_flexbar" -r ABH8_1_clipped_qf.fastq.gz -q TAIL -qf i1.8 -qt 30 -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 5 -z GZ &
flexbar -t $TRIMDIR"ABH8_2_RNAseq_flexbar" -r ABH8_2_clipped_qf.fastq.gz -q TAIL -qf i1.8 -qt 30 -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 5 -z GZ &
flexbar -t $TRIMDIR"ABH8_3_RNAseq_flexbar" -r ABH8_3_clipped_qf.fastq.gz -q WIN -q TAIL -qf i1.8 -qt 30 -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 5 -z GZ &
flexbar -t $TRIMDIR"ABH8_1_CLIPseq_flexbar" -r ABH8_1_flexbar3.fastq.gz -q WIN -q TAIL -qf i1.8 -qt 30 -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 5 -z GZ &
flexbar -t $TRIMDIR"ABH8_2_CLIPseq_flexbar" -r ABH8_2_flexbar3.fastq.gz -q WIN -q TAIL -qf i1.8 -qt 30 -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 5 -z GZ &
flexbar -t $TRIMDIR"ABH8_3_CLIPseq_flexbar" -r ABH8_3_flexbar3.fastq.gz -q WIN -q TAIL -qf i1.8 -qt 30 -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 5 -z GZ
#Wait for all jobs to finish before continuing the script
wait
echo "Done processing ABH8 CLIP and RNAseq SE reads - Flexbar ANY and quality Trimming"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"

cd $TRIMDIR

#bylo nutne upravit volani souboru pro fastqc
fastqc --outdir $QCDIR --format fastq --threads $THREADS *.gz
#
cd $QCDIR
echo "Number of reads after preprocessing is in fastqc/trimmed/total_sequences.txt"
for i in *zip
do
	unzip -q $i
	grep -hF "Total Sequences" ${i%.zip*}/*.txt
done > total_sequences.txt
rm -r ./*fastqc/

############################################################################################
### Copy data from scratch back to home dir and clean scratch
cd $SCRATCHDIR
rm $SCRATCHDIR/*$APPENDIX
cp -avr $SCRATCHDIR $OUTPUT_DIR || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
