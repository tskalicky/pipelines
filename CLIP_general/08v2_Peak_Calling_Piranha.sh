#!/bin/bash
#PBS -l walltime=24:0:0 
#PBS -q default
#PBS -l select=1:ncpus=16:mem=150gb:scratch_local=150gb
#PBS -j oe
#PBS -N 07v3_peak_calling_FTO
#
## initialize the required application
# module add gsl-1.16-intel #required by Piranha
# module add bamtools #required by Piranha
module add piranha-1.2.1
module add bedtools-2.26.0 
############################################################################################
#	HITS-CLIP Peak Calling using Piranha and reads after PCR collapsing using Picard tools
#	Requires: Piranha, bedtools
#
############################################################################################
### Variables
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/FTO/peak_calling"
INPUT_DIR="/storage/brno3-cerit/home/tskalicky/FTO/mapping/PCR_collapsed_Picard/pooled_BAM"
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
cp -av $INPUT_DIR/*.bam $SCRATCHDIR
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR
echo "Following files/folders were copied to scratch:"
ls -c1

# Binaries
PIRANHA=$(which Piranha)
BEDTOOLS=$(which bedtools)
# Check if tools are installed
echo "Used tools are installed here:"
which $PIRANHA
which $BEDTOOLS
####################################################################################################
### Commands
## converting BAM to BED
date +"%d/%m/%Y %H:%M:%S" 
echo "Converting all BAM files to BED."
for a in *bam
do
	FILE=$a
	FILENAME=${a%.*}
	$BEDTOOLS bamtobed -i $FILE > $FILENAME.bed &
done
date +"%d/%m/%Y %H:%M:%S" 
echo "Finnished conversion of all BAM files to BED." 
#
date +"%d/%m/%Y %H:%M:%S" 
echo "Started Peak Calling"
$PIRANHA -v -s -z 20 -p 0.05 FTO_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.bam FTO_WT1-3_GRCh38.pcr_dedupl.rg.merged.bed -o FTO_CLIP1-3_GRCh38_piranha_n20_vs_WT_pool_peaks &
$PIRANHA -v -s -z 20 -p 0.05 FTO_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.bam FTO_KO1-3_GRCh38.pcr_dedupl.rg.merged.bed -o FTO_CLIP1-3_GRCh38_piranha_n20_vs_KO_pool_peaks &
$PIRANHA -v -s -z 20 -p 0.05 FTO_INPUT1-3_GRCh38.pcr_dedupl.rg.merged.bam FTO_KO1-3_GRCh38.pcr_dedupl.rg.merged.bed -o FTO_INPUT1-3_GRCh38_piranha_n20_vs_KO_pool_peaks &
$PIRANHA -v -s -z 20 -p 0.05 FTO_INPUT1-3_GRCh38.pcr_dedupl.rg.merged.bam FTO_WT1-3_GRCh38.pcr_dedupl.rg.merged.bed -o FTO_INPUT1-3_GRCh38_piranha_n20_vs_WT_pool_peaks &
#
$PIRANHA -v -s -z 50 -p 0.05 FTO_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.bam FTO_WT1-3_GRCh38.pcr_dedupl.rg.merged.bed -o FTO_CLIP1-3_GRCh38_piranha_n50_vs_WT_pool_peaks &
$PIRANHA -v -s -z 50 -p 0.05 FTO_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.bam FTO_KO1-3_GRCh38.pcr_dedupl.rg.merged.bed -o FTO_CLIP1-3_GRCh38_piranha_n50_vs_KO_pool_peaks &
$PIRANHA -v -s -z 50 -p 0.05 FTO_INPUT1-3_GRCh38.pcr_dedupl.rg.merged.bam FTO_KO1-3_GRCh38.pcr_dedupl.rg.merged.bed -o FTO_INPUT1-3_GRCh38_piranha_n50_vs_KO_pool_peaks &
$PIRANHA -v -s -z 50 -p 0.05 FTO_INPUT1-3_GRCh38.pcr_dedupl.rg.merged.bam FTO_WT1-3_GRCh38.pcr_dedupl.rg.merged.bed -o FTO_INPUT1-3_GRCh38_piranha_n50_vs_WT_pool_peaks &
#
$PIRANHA -v -s -z 20 -p 0.05 FTO_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.bam FTO_WT1_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed \
FTO_WT2_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed FTO_WT3_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed \
-o FTO_CLIP1-3_GRCh38_piranha_n20_pool_vs_WT_separate_peaks &
$PIRANHA -v -s -z 20 -p 0.05 FTO_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.bam FTO_KO1_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed \
FTO_KO2_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed FTO_KO3_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed \
-o FTO_CLIP1-3_GRCh38_piranha_n20_pool_vs_KO_separate_peaks &
#
$PIRANHA -v -s -z 20 -p 0.05 FTO_INPUT1-3_GRCh38.pcr_dedupl.rg.merged.bam FTO_WT1_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed \
FTO_WT2_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed FTO_WT3_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed \
-o FTO_CLIP1-3_GRCh38_piranha_n20_pool_vs_WT_separate_peaks &
$PIRANHA -v -s -z 20 -p 0.05 FTO_INPUT1-3_GRCh38.pcr_dedupl.rg.merged.bam FTO_KO1_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed \
FTO_KO2_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed FTO_KO3_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed \
-o FTO_CLIP1-3_GRCh38_piranha_n20_pool_vs_KO_separate_peaks &
#
$PIRANHA -v -s -z 50 -p 0.05 FTO_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.bam FTO_WT1_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed \
FTO_WT2_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed FTO_WT3_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed \
-o FTO_CLIP1-3_GRCh38_piranha_n50_pool_vs_WT_separate_peaks &
$PIRANHA -v -s -z 50 -p 0.05 FTO_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.bam FTO_KO1_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed \
FTO_KO2_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed FTO_KO3_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed \
-o FTO_CLIP1-3_GRCh38_piranha_n50_pool_vs_KO_separate_peaks &
#
$PIRANHA -v -s -z 50 -p 0.05 FTO_INPUT1-3_GRCh38.pcr_dedupl.rg.merged.bam FTO_WT1_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed \
FTO_WT2_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed FTO_WT3_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed \
-o FTO_CLIP1-3_GRCh38_piranha_n50_pool_vs_WT_separate_peaks &
$PIRANHA -v -s -z 50 -p 0.05 FTO_INPUT1-3_GRCh38.pcr_dedupl.rg.merged.bam FTO_KO1_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed \
FTO_KO2_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed FTO_KO3_GRCh38.91_Aligned.pcr_dedupl.sorted.md.rg.bed \
-o FTO_CLIP1-3_GRCh38_piranha_n50_pool_vs_KO_separate_peaks &
#wait
wait
date +"%d/%m/%Y %H:%M:%S" 
echo "Finnished Peak Calling"

############################################################################################
### Copy data from scratch back to home dir and clean scratch
mkdir -p $OUTPUT_DIR
cp -avr $SCRATCHDIR $OUTPUT_DIR || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"

