#!/bin/bash
#PBS -l walltime=96:0:0 
#PBS -q uv@wagap-pro.cerit-sc.cz
#PBS -l select=1:ncpus=15:mem=150gb:scratch_local=150gb
#PBS -j oe
#PBS -N 02v09_flexbar_fastqc_Dis3L2
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
INPUT_D3L2_WT="/storage/brno3-cerit/home/tskalicky/D3L2/RAW/D3L2-WT*fastq.gz"
INPUT_D3L2_KO="/storage/brno3-cerit/home/tskalicky/D3L2/RAW/D3L2-KO*fastq.gz"
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/ABH8/preprocessed"
#
ADAPTERS="/storage/brno3-cerit/home/tskalicky/tools/adapters_v4.fa"
APPENDIX=".fastq.gz"
THREADS=$PBS_NUM_PPN # Number of threads to use
#
FASTQC=$(which fastqc)
FLEXBAR=$(which flexbar)
TRIMMOMATIC="java -jar /software/trimmomatic/0.36/dist/jar/trimmomatic-0.36.jar"
# Check the tools versions
echo $FASTQC
echo $FLEXBAR
echo $TRIMMOMATIC
#
####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
cd /storage/brno3-cerit/home/tskalicky/D3L2/RAW
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
cp -av $INPUT_D3L2_WT $INPUT_D3L2_KO $ADAPTERS $SCRATCHDIR
#
cd $SCRATCHDIR
#
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
#FLEXBAR is using only like 1/2 the cores assigned (no more like 3-4 cores)
# POZOR! Flexbar nebezi na Metacentru pri pouziti uzlu s os=debian9 !
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Now I am processing Dis3L2 WT and KO PE reads - Flexbar ANY and quality Trimming"
# Zvoleno trimovani na kvalitu az po trimovani adapteru!
flexbar -t $TRIMDIR"Dis3L2_1_WT_flexbar" -r D3L2-KO1_FW.fastq.gz -p D3L2-KO1_RV.fastq.gz -q TAIL -qf i1.8 -qt 30 --qtrim-post-removal -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 4 -z GZ &
flexbar -t $TRIMDIR"Dis3L2_2_WT_flexbar" -r D3L2-KO2_FW.fastq.gz -p D3L2-KO2_RV.fastq.gz -q TAIL -qf i1.8 -qt 30 --qtrim-post-removal -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 4 -z GZ &
flexbar -t $TRIMDIR"Dis3L2_3_WT_flexbar" -r D3L2-KO3_FW.fastq.gz -p D3L2-KO3_RV.fastq.gz -q WIN -q TAIL -qf i1.8 -qt 30 --qtrim-post-removal -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 4 -z GZ &
flexbar -t $TRIMDIR"Dis3L2_1_KO_flexbar" -r D3L2-WT1_FW.fastq.gz -p D3L2-WT1_RV.fastq.gz -q WIN -q TAIL -qf i1.8 -qt 30 --qtrim-post-removal -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 4 -z GZ &
flexbar -t $TRIMDIR"Dis3L2_2_KO_flexbar" -r D3L2-WT2_FW.fastq.gz -p D3L2-WT2_RV.fastq.gz -q WIN -q TAIL -qf i1.8 -qt 30 --qtrim-post-removal -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 4 -z GZ &
flexbar -t $TRIMDIR"Dis3L2_3_KO_flexbar" -r D3L2-WT3_FW.fastq.gz -p D3L2-WT3_RV.fastq.gz -q WIN -q TAIL -qf i1.8 -qt 30 --qtrim-post-removal -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 4 -z GZ
#Wait for all jobs to finish before continuing the script
wait
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Done processing Dis3L2 WT and KO PE reads - Flexbar ANY and quality Trimming"

cd $TRIMDIR
#bylo nutne upravit volani souboru pro fastqc
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Started FastQC evaluation for all trimmed libraries"
fastqc --outdir $QCDIR --format fastq --threads $THREADS *.gz
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Finnished FastQC evaluation for all trimmed libraries"
#
cd $QCDIR
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Creating library summary info for all trimmed libraries"
echo "Number of reads after preprocessing is in fastqc/trimmed/total_sequences.txt"
for i in *zip
do
	unzip -q $i
	grep -hF "Total Sequences" ${i%.zip*}/*.txt
done > total_sequences.txt
rm -r ./*fastqc/
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Creating library summary info for all trimmed libraries"
echo "Number of reads in libraries are:"
cat total_sequences.txt

############################################################################################
### Copy data from scratch back to home dir and clean scratch
cd $SCRATCHDIR
rm $SCRATCHDIR/*$APPENDIX
cp -avr $SCRATCHDIR $OUTPUT_DIR || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
