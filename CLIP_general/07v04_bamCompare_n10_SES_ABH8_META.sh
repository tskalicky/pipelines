#!/bin/bash
#PBS -l walltime=4:0:0 
#PBS -q  uv@wagap-pro.cerit-sc.cz
#PBS -l select=1:ncpus=10:mem=100gb:scratch_local=50gb
#PBS -j oe
#PBS -N 07v04_bamCompare_ABH8
#
## initialize the required application
module add samtools-1.4
module add bedtools-2.25.0 
############################################################################################
# CLIP Seq pipeline
# 1) Sorting and index creation for BED files after PCR colapsing and bedgraph creation
# 2) Comparing CLIPseq vs RNAseq using Deeptools bamCompare
#
############################################################################################
### Variables
# 
# Requires manually installed Deeptools through Anaconda!!!
############################################################################################
### Variables
THREADS=$PBS_NUM_PPN
# Need to export manually installed binaries
export PATH="/storage/brno3-cerit/home/tskalicky/anaconda2/bin:$PATH"
#
INPUT="/storage/brno3-cerit/home/tskalicky/ABH8/mapping/PCR_collapsed_Picard/pooled"
CLIP="ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.sorted.bam"
RNASEQ="ABH8_RNAseq1-3_GRCh38.pcr_dedupl.rg.merged.sorted.bam"
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/ABH8/bamCompare"
############################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
cp -av $INPUT/*.sorted.bam $INPUT/*.sorted.bam.bai $SCRATCHDIR
cd $SCRATCHDIR
#
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR
echo "Following files/folders were copied to scratch:"
ls -c1
# Binaries
BAMCOMPARE=$(which bamCompare)
# Check if tools are installed
echo "Required software is installed here:"
which $BAMCOMPARE
# Commands
echo "Counting number of uniquely mapped reads in each sequencing library"
for a in *.bam
do
	SAMPLE=$a
 	SAMPLENAME=${a%.*}
 	date +"%d/%m/%Y %H:%M:%S" 
 	echo "Number of uniquely mapped reads in $SAMPLE is:"
	# For Single-end libraries
	# Excluding those marked as secondary or supplementary (i.e., multimappers or chimeric entries)
	samtools view -F 0x904 -c $SAMPLE
	# For Pair-end libraries
	# Excluding those marked as secondary or supplementary (i.e., multimappers or chimeric entries)
	# samtools view -F 0x4 $SAMPLE | cut -f 1 | sort | uniq | wc -l
done
#wait
wait
# cd $INPUT
# for a in *.bam.gz
# do
# 	FILE=$a
# 	FILENAME=${a%.*.*}
# 	date +"%d/%m/%Y %H:%M:%S" 
# 	echo "Now I am decompressing library $FILE"
# 	unpigz -v -p $THREADS $FILE
# 	echo "Done decompressing PE reads $FILE"
# done
# # wait
# wait
# for b in *.sorted.bam
# do
# 	FILE=$b
# 	FILENAME=${b%.*}
# 	date +"%d/%m/%Y %H:%M:%S" 
# 	#echo "Now I am sorting file $FILE"
# 	#igvtools "sort" $FILE $FILENAME".sorted.bam"
# 	#echo "Finnished sorting file $FILE"
# 	date +"%d/%m/%Y %H:%M:%S" 
# 	echo "Started indexing file $FILE"
# 	igvtools "index" $FILE
# done
# wait
## Compare 2 datasets
# mkdir /home/tomas/CEITEC_lab/ABH8/bamCompare
date +"%d/%m/%Y %H:%M:%S" 
echo "Now I am comparing files using bamCompare SES enhance scalling subtract normalisation"
bamCompare \
-b1 ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.sorted.bam \
-b2 ABH8_RNAseq1-3_GRCh38.pcr_dedupl.rg.merged.sorted.bam \
--numberOfProcessors $THREADS \
--binSize 10 \
--scaleFactorsMethod 'SES' \
--sampleLength 10000 \
--numberOfSamples 100000 \
--operation 'subtract' \
--skipNonCoveredRegions \
-of bedgraph \
-o ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_subtract_nosmooth.bed
date +"%d/%m/%Y %H:%M:%S" 
echo "bamCompare SES enhance scalling subtract normalisation is finnished"
#
date +"%d/%m/%Y %H:%M:%S" 
echo "Now I am comparing files using bamCompare SES enhance scalling subtract smooth normalisation"
bamCompare \
-b1 ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.sorted.bam \
-b2 ABH8_RNAseq1-3_GRCh38.pcr_dedupl.rg.merged.sorted.bam \
--numberOfProcessors $THREADS \
--binSize 10 \
--scaleFactorsMethod 'SES' \
--sampleLength 10000 \
--numberOfSamples 100000 \
--operation 'subtract' \
--smoothLength 30 \
--skipNonCoveredRegions \
-of bedgraph \
-o ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_subtract_smooth.bed
date +"%d/%m/%Y %H:%M:%S" 
echo "bamCompare SES enhance scalling subtract smooth normalisation is finnished"
#
date +"%d/%m/%Y %H:%M:%S" 
echo "Now I am comparing files using bamCompare RPKM normalisation ONLY"
bamCompare \
-b1 ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.sorted.bam \
-b2 ABH8_RNAseq1-3_GRCh38.pcr_dedupl.rg.merged.sorted.bam \
--numberOfProcessors $THREADS \
--binSize 10 \
--scaleFactorsMethod 'None' \
--operation 'subtract' \
--normalizeUsing 'RPKM' \
--skipNonCoveredRegions \
-of bedgraph \
-o ABH8_CLIP_vs_RNAseq_bamCompare_bin10_RPKM_nosmooth.bed
date +"%d/%m/%Y %H:%M:%S" 
echo "bamCompare SES enhance subtract RPKM smooth normalisation ONLY is finnished"
#
date +"%d/%m/%Y %H:%M:%S" 
echo "Now I am comparing files using bamCompare RPKM smooth normalisation ONLY smooth n=30"
bamCompare \
-b1 ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.sorted.bam \
-b2 ABH8_RNAseq1-3_GRCh38.pcr_dedupl.rg.merged.sorted.bam \
--numberOfProcessors $THREADS \
--binSize 10 \
--scaleFactorsMethod 'None' \
--operation 'subtract' \
--normalizeUsing 'RPKM' \
--smoothLength 30 \
--skipNonCoveredRegions \
-of bedgraph \
-o ABH8_CLIP_vs_RNAseq_bamCompare_bin10_RPKM_smooth.bed
date +"%d/%m/%Y %H:%M:%S" 
echo "bamCompare RPKM smooth normalisation is finnished ONLY smooth n=30"

############################################################################################
### Copy data from scratch back to home dir and clean scratch
mkdir -p $OUTPUT_DIR
rm $SCRATCHDIR/*.sorted.bam $SCRATCHDIR/*.bai
cp -avr $SCRATCHDIR $OUTPUT_DIR || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"