#!/bin/bash
#PBS -l walltime=24:0:0 
#PBS -q default@wagap-pro.cerit-sc.cz 
#PBS -l select=1:ncpus=5:mem=50gb:scratch_local=100gb
#PBS -j oe
#PBS -N 01_BAM_to_FASTQ_ALKBH8_RIPseq
#export PBS_SERVER=wagap-pro.cerit-sc.cz # needed only when executing from arien frontends
#
## initialize the required application
module add bamtools
module add picard-2.9.0 # will also initialize system variable $PICARD pointing into Picard Tools install dir.
module add samtools-1.4
############################################################################################
# CLIPseq analysis
# Will keep only unique mappings (with MAPQ >=1) and a minimal mapping size of 18 nt.
# Afterwards collapses PCR duplicates
############################################################################################
### Variables
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/ABH8/RIP-seq/RAW"
INPUT="/storage/brno3-cerit/home/tskalicky/ABH8/RIP-seq/RAW"

THREADS=$PBS_NUM_PPN
APPENDIX=".bam"
#APPENDIX2=".fq"
#APPENDIX3=".sam"
####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
cp -av $INPUT/*.bam $SCRATCHDIR
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR
echo "Following files/folders were copied to scratch:"
ls -c1

# Binaries
SAMTOOLS=$(which samtools)

# Check if tools are installed
which $PICARD
which $SAMTOOLS
####################################################################################################
### commands

for a in *.bam
do
	FILE=$a
	FILENAME=${a%.*}
	date +"%d/%m/%Y %H:%M:%S"
	echo "Converting file $FILE to fastq"
	java -XX:ParallelGCThreads=$THREADS -Xmx50g -jar $PICARD SamToFastq INPUT=$FILE FASTQ="$FILENAME"_R1.fastq SECOND_END_FASTQ="$FILENAME"_R2.fastq
	date +"%d/%m/%Y %H:%M:%S"
	echo "Finnished converting file $FILE to fastq"
done &
#wait
wait
############################################################################################
### Copy data from scratch back to home dir and clean scratch
mkdir -p $OUTPUT_DIR
rm *.bam
cp -avr $SCRATCHDIR $OUTPUT_DIR  || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"