#!/bin/bash
#PBS -l walltime=24:0:0 
#PBS -q default@wagap-pro.cerit-sc.cz 
#PBS -l select=1:ncpus=10:mem=200gb:scratch_local=50gb:os=debian9
#PBS -j oe
#PBS -N SAM_to_BAM
#export PBS_SERVER=wagap-pro.cerit-sc.cz # needed only when executing from arien frontends
#
## initialize the required application
module add bamtools
module add samtools-1.4
############################################################################################
# SAM to BAM conversion
############################################################################################
### Variables
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/FTO/mapping/PCR_collapsed_Picard"
INPUT_DIR="/storage/brno3-cerit/home/tskalicky/FTO/mapping/PCR_collapsed_Picard"

THREADS=$PBS_NUM_PPN
APPENDIX=".pcr_dedupl.md.sam"
#APPENDIX2=".fq"
#APPENDIX3=".sam"
####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
cp -av $INPUT_DIR/*.pcr_dedupl.md.sam $SCRATCHDIR
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR
echo "Following files/folders were copied to scratch:"
ls -c1

# Binaries
SAMTOOLS=$(which samtools)

# Check if tools are installed
which $SAMTOOLS
####################################################################################################
### commands

for a in *.sam
do
	FILE=$a
	FILENAME=${a%.*.*.*}
	$SAMTOOLS view -@ $THREADS -bS $FILENAME.pcr_dedupl.md.sam | $SAMTOOLS sort -@ $THREADS --output-fmt BAM - > $FILENAME.pcr_dedupl.sorted.md.bam
    $SAMTOOLS index -@ $THREADS $FILENAME.pcr_dedupl.sorted.md.bam $FILENAME.pcr_dedupl.sorted.md.bai
done
#wait
wait
############################################################################################
### Copy data from scratch back to home dir and clean scratch
#mkdir -p $OUTPUT_DIR
rm *.sam
cp -avr $SCRATCHDIR $OUTPUT_DIR  || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"