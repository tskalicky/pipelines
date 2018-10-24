#!/bin/bash
#PBS -l walltime=24:0:0 
#PBS -q default@wagap-pro.cerit-sc.cz 
#PBS -l select=1:ncpus=20:mem=300gb:scratch_local=50gb:os=debian9
#PBS -j oe
#PBS -N 05_CLIPseq_Picard_FTO_PCR_duplicate_collapse
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
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/FTO/mapping/PCR_collapsed_Picard"
FTO_INPUT="/storage/brno3-cerit/home/tskalicky/FTO/mapping/PCR_collapsed-genome/FTO_Input/SAM"
FTO_CLIP="/storage/brno3-cerit/home/tskalicky/FTO/mapping/PCR_collapsed-genome/FTO_CLIP/SAM"
FTO_KO="/storage/brno3-cerit/home/tskalicky/FTO/mapping/PCR_collapsed-genome/FTO_KO1-3/SAM"
FTO_WT="/storage/brno3-cerit/home/tskalicky/FTO/mapping/PCR_collapsed-genome/FTO_WT1-3/SAM"

THREADS=$PBS_NUM_PPN
APPENDIX=".sam"
#APPENDIX2=".fq"
#APPENDIX3=".sam"
####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
cp -av $FTO_INPUT/*.sam $FTO_CLIP/*.sam $FTO_KO/*.sam $FTO_WT/*.sam $SCRATCHDIR
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

for a in *.sam
do
	FILE=$a
	FILENAME=${a%.*.*}
	java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true \
    I=$FILE \
    O=$FILENAME.pcr_dedupl.md.sam \
    M=$FILENAME.pcr_dedupl_metrics.txt
    $SAMTOOLS view -bS $FILENAME.pcr_dedupl.md.sam | $SAMTOOLS sort - $FILENAME.pcr_dedupl.sorted.md.bam
    $SAMTOOLS index $FILENAME.pcr_dedupl.sorted.md.bam $FILENAME.pcr_dedupl.sorted.md.bai
done
#wait
wait
############################################################################################
### Copy data from scratch back to home dir and clean scratch
mkdir -p $OUTPUT_DIR
cp -avr $SCRATCHDIR $OUTPUT_DIR  || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"