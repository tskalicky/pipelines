#!/bin/bash
#PBS -l walltime=96:0:0 
#PBS -q default@wagap-pro.cerit-sc.cz 
#PBS -l select=1:ncpus=19:mem=40gb:scratch_local=200gb:os=debian9
#PBS -j oe
#PBS -N 02v7_flexbar_left_cut_METTL16
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


INPUT_FlagPAR_UV="/storage/brno3-cerit/home/tskalicky/METTL16/preprocessed/data/FlagPAR_FlagUV_jobs/FlagPAR_FlagUV_3_job_564502_wagap"
INPUT_METTL16_PAR1="/storage/brno3-cerit/home/tskalicky/METTL16/preprocessed/data/METTL16_PAR1_job_4260487_arien/METTL16PAR-1_SRR6048658_flexbar_qual.fastq.gz"
INPUT_METTL16_UV1="/storage/brno3-cerit/home/tskalicky/METTL16/preprocessed/data/METTL16_UV1_job_4260489_arien/METTL16UV-1_SRR6048655_flexbar_qual.fastq.gz"
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/METTL16/preprocessed"

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
cd /storage/brno3-cerit/home/tskalicky/METTL16/preprocessed/data
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
cp -av $INPUT_FlagPAR_UV/*$APPENDIX $INPUT_METTL16_PAR1 $INPUT_METTL16_UV1 $ADAPTERS $SCRATCHDIR

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

#usage of $FLEXBAR and $FASTQC variables in function is not working in PBS "Illegal instruction" error
#loop over every file and process them in background simultaniously
#FLEXBAR is using only like 1/2 the cores assigned
echo "Now I am processing METTL16 SE reads - Flexbar pre-trim-left Trimming"
flexbar -t $TRIMDIR"FlagPAR_SRR6048657_left_cut_flexbar" -r FlagPAR_SRR6048657_flexbar_qual3.fastq.gz --pre-trim-left 10 --min-read-length 17 -n 5 -z GZ &
flexbar -t $TRIMDIR"FlagUV_SRR6048654_left_cut_flexbar" -r FlagUV_SRR6048654_flexbar_qual3.fastq.gz --pre-trim-left 10 --min-read-length 17 -n 5 -z GZ &
flexbar -t $TRIMDIR"METTL16PAR-1_SRR6048658_left_cut_flexbar" -r METTL16PAR-1_SRR6048658_flexbar_qual.fastq.gz --pre-trim-left 10 --min-read-length 17 -n 5 -z GZ &
flexbar -t $TRIMDIR"METTL16UV-1_SRR6048655_left_cut_flexbar" -r METTL16UV-1_SRR6048655_flexbar_qual.fastq.gz --pre-trim-left 10 --min-read-length 17 -n 5 -z GZ
# Wait for all jobs to finish before continuing the script
wait
echo "Done processing METTL16 SE reads - Flexbar pre-trim-left Trimming"
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
