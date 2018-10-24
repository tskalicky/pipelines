#!/bin/bash
#PBS -l walltime=96:0:0 
#PBS -q default@wagap-pro.cerit-sc.cz 
#PBS -l select=1:ncpus=6:mem=200gb:scratch_local=250gb:os=debian9
#PBS -j oe
#PBS -N 05_FTO_CLIP_Collaps_PCR_dupl
#
## initialize the required application
#module add star-2.5.2b
#module add samtools-1.6
#module add python27-modules-intel #required by multiQC
#module add bedtools-2.26.0
#module add ctk-1.0.7
# module add gsl-1.16-intel #required by Piranha
# module add bamtools #required by Piranha
#
# Alignment - SE RNA-Seq
# Pipeline to align RNA-Seq experiment as recommended in GENCODE project https://github.com/ENCODE-DCC/long-rna-seq-pipeline/tree/master/dnanexus/align-star-pe; https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/dnanexus/Readme.md with few moditifications: --twopassMode Basic; --outSAMattributes All --outFilterMismatchNoverReadLmax 0.05
# To increase sensitivity you might try to add --seedSearchStartLmax 30
#
# Requires STAR, pigz, samtools, multiQC, bedGraphToBigWig
#
# Note: do not use modified GTF (added features) to alignment as it can cause issues later on
############################################################################################
### Variables
INPUT_DIR=/Users/tomas/Data/FTO/mapping/FTO_CLIP_INPUT_results/alignment/genome/SAM
OUTPUT_DIR=/Users/tomas/Data/FTO/mapping/FTO_CLIP_INPUT_results/alignment/genome/BED

MY_RAM=50 # Max RAM memory for Samtools sort
THREADS=4
APPENDIX=".out.bam"
APPENDIX2=".sam"
APPENDIX3=".tag.bed"

# Binaries
#STAR=$(which STAR)
#SAMTOOLS=$(which samtools)
#MULTIQC=$(which multiqc)
#BEDTOOLS=$(which bedtools)
TAG2COLLAPSE=$(which tag2collapse.pl)
PARSEALIGN=$(which parseAlignment.pl)

# Check the tools versions
echo "Used software is located here:"
#which $STAR
#which $SAMTOOLS
#which $MULTIQC
#which $BEDTOOLS
which $TAG2COLLAPSE
which $PARSEALIGN

### Commands
# mkdir -p $SCRATCHDIR/{SAM,BED}
#for i in *$APPENDIX2
#do
#	READ_FOR=$i
#	SAMPLENAME=${i%.*}
#	echo "Now I am converting BAM file $READ_FOR to SAM"
#	$SAMTOOLS view -h $READ_FOR > $SCRATCHDIR/SAM/$SAMPLENAME".sam" &
#done
## Wait for all jobs to finish before exiting the job submission script
#wait
#echo "Done convering ALL BAM to SAM files"
mkdir $OUTPUT_DIR
cd $INPUT_DIR
for a in *$APPENDIX2
do
	READ_FOR2=$a
	SAMPLENAME2=${a%.*}
	echo "Now I am parsing SAM file $READ_FOR2 for unique mapppings only"
	perl $PARSEALIGN -v --map-qual 1 --min-len 18 --mutation-file $INPUT_DIR/BED/$SAMPLENAME2".mutation.txt" $READ_FOR2 $INPUT_DIR/BED/$SAMPLENAME2".tag.bed"
done
# Wait for all jobs to finish before exiting the job submission script
wait
echo "Done parsing ALL SAM files for unique mappings only"
cd $INPUT_DIR/BED
wc -l $OUTPUT_DIR/*.tag.bed > Uniquely_mapped_read_proportions.txt # Keep track what proportion of reads can be mapped uniquely.
#
for b in *APPENDIX3
do
	READ_FOR3=$b
	SAMPLENAME3=${a%.*.*}
	echo "Now I am collapsing PCR duplicates for $READ_FOR3 file"
	perl $TAG2COLLAPSE -v -big -EM 30 --seq-error-model alignment \
	-weight --weight-in-name --keep-max-score --keep-tag-name $b $SAMPLENAME3".tag.uniq.bed"
done
# Wait for all jobs to finish before exiting the job submission script
wait
echo "Done collapsing PCR duplicates for ALL files"
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"







