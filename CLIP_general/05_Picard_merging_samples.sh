#!/bin/bash
#PBS -l walltime=4:0:0 
#PBS -q default@wagap-pro.cerit-sc.cz 
#PBS -l select=1:ncpus=2:mem=100gb:scratch_local=100gb:os=debian9
#PBS -j oe
#PBS -N 05_CLIPseq_Picard_pooling_samples
#export PBS_SERVER=wagap-pro.cerit-sc.cz # needed only when executing from arien frontends
#
## initialize the required application
module add bamtools
module add picard-2.9.0 # will also initialize system variable $PICARD pointing into Picard Tools install dir.
module add samtools-1.4
############################################################################################
# CLIPseq analysis
# Will add/modify @RG header in SAM/BAM files that they can bee pooled for loading in IGV
############################################################################################
### Variables
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/FTO/mapping/PCR_collapsed_Picard/pooled_BAM"
INPUT_DIR="/storage/brno3-cerit/home/tskalicky/FTO/mapping/PCR_collapsed_Picard/BAM"

THREADS=$PBS_NUM_PPN
APPENDIX=".bam"
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
SAMTOOLS=$(which samtools)

# Check if tools are installed
which $PICARD
which $SAMTOOLS
####################################################################################################
### commands
counter=0
for a in *.bam
do
	let counter++
	FILE=$a
	FILENAME=${a%.*.*}
	SAMPLE=$(echo $FILE | cut -d"_" -f1,2)
	echo "Adding Read Groups to file $FILE" 
	java -jar $PICARD AddOrReplaceReadGroups \
    I=$FILE \
    O=$FILENAME.md.rg.bam \
    RGID=$SAMPLE \
    RGLB=lib"$counter" \
    RGPL=illumina \
    RGPU=unit"$counter" \
    RGSM=$SAMPLE
    echo "Done adding Read Groups to file $FILE"
done
#wait
wait
# declare an array variables
declare -a FTO_CLIP=(FTO_CLIP*.md.rg.bam)
declare -a FTO_INPUT=(FTO_Input*.md.rg.bam)
declare -a FTO_KO=(FTO_KO*.md.rg.bam)
declare -a FTO_WT=(FTO_WT*.md.rg.bam)
#
$SAMTOOLS merge FTO_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.bam ${FTO_CLIP[@]}
$SAMTOOLS merge FTO_INPUT1-3_GRCh38.pcr_dedupl.rg.merged.bam ${FTO_INPUT[@]}
$SAMTOOLS merge FTO_KO1-3_GRCh38.pcr_dedupl.rg.merged.bam ${FTO_KO[@]}
$SAMTOOLS merge FTO_WT1-3_GRCh38.pcr_dedupl.rg.merged.bam ${FTO_WT[@]}

for b in *.rg.merged.bam
do
	FILE2=$b
	FILENAME2=${b%.*}
	$SAMTOOLS index -@ $THREADS $FILE2 $FILENAME.bai
done
#wait
wait

############################################################################################
### Copy data from scratch back to home dir and clean scratch
mkdir -p $OUTPUT_DIR
rm *.md.bam
cp -avr $SCRATCHDIR $OUTPUT_DIR  || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"