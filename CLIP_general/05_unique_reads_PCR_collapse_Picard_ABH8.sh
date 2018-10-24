#!/bin/bash
#PBS -l walltime=24:0:0 
#PBS -q default@wagap-pro.cerit-sc.cz 
#PBS -l select=1:ncpus=4:mem=60gb:scratch_local=50gb:os=debian9
#PBS -j oe
#PBS -N 05_CLIPseq_Picard_ABH8_PCR_duplicate_collapse
#export PBS_SERVER=wagap-pro.cerit-sc.cz # needed only when executing from arien frontends
#
## initialize the required application
module add bamtools
module add picard-2.9.0 # will also initialize system variable $PICARD pointing into Picard Tools install dir.
module add samtools-1.4
############################################################################################
# CLIPseq analysis
# Will keep only unique mappings (with MAPQ >=1) and a minimal mapping size of 17 nt.
# Afterwards collapses PCR duplicates
# PICARD tools accept both SAM or BAM input
############################################################################################
### Variables
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/ABH8/mapping/PCR_collapsed_Picard"
ABH8_CLIP="/storage/brno3-cerit/home/tskalicky/ABH8/mapping/deduplicated/CLIPseq/alignment/genome"
ABH8_RNAseq="/storage/brno3-cerit/home/tskalicky/ABH8/mapping/deduplicated/RNAseq/alignment/genome"
#
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
cp -av $ABH8_CLIP/*.bam $ABH8_RNAseq/*.bam $SCRATCHDIR
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
# Removing PCR duplicates based on their alignment
date +"%d/%m/%Y %H:%M:%S" 
echo "Removing PCR duplicates based on their alignment."
for a in *.bam
do
	FILE=$a
	FILENAME=${a%.*.*}
	java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true \
    I=$FILE \
    O=$FILENAME.pcr_dedupl.bam \
    M=$FILENAME.pcr_dedupl_metrics.txt
    $SAMTOOLS view -bS $FILENAME.pcr_dedupl.bam | $SAMTOOLS sort -@ $THREADS - -o $FILENAME.pcr_dedupl.sorted.bam
    $SAMTOOLS index $FILENAME.pcr_dedupl.sorted.bam $FILENAME.pcr_dedupl.sorted.bai
done
#wait
wait
date +"%d/%m/%Y %H:%M:%S" 
echo "Finnished removing PCR duplicates based on their alignment."
#
# Merging several coresponding libraries from biological replicates together in one file
#
# Adding RG Group sign to each file for merging later
date +"%d/%m/%Y %H:%M:%S" 
echo "Adding RG Group signs to each file."
#
counter=0
for b in *.pcr_dedupl.sorted.bam
do
	let counter++
	FILE2=$b
	FILENAME2=${b%.*}
	SAMPLE=$(echo $FILE2 | cut -d"_" -f1,2)
	echo "Adding Read Groups to file $FILE2" 
	java -jar $PICARD AddOrReplaceReadGroups \
    I=$FILE2 \
    O=$FILENAME2.rg.bam \
    RGID=$SAMPLE \
    RGLB=lib"$counter" \
    RGPL=illumina \
    RGPU=unit"$counter" \
    RGSM=$SAMPLE
    echo "Done adding Read Groups to file $FILE2"
done
#wait
wait
date +"%d/%m/%Y %H:%M:%S" 
echo "Finnished adding RG Group signs to each file."
# declare an array variables
declare -a ABH8_CLIP=(ABH8-*_CLIPseq*.rg.bam)
declare -a ABH8_RNAseq=(ABH8-*_RNAseq*.rg.bam)
#
# Merging biological replicates
date +"%d/%m/%Y %H:%M:%S" 
echo "Merging biological replicates."
$SAMTOOLS merge ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.bam ${ABH8_CLIP[@]} &
$SAMTOOLS merge ABH8_RNAseq1-3_GRCh38.pcr_dedupl.rg.merged.bam ${ABH8_RNAseq[@]}
#wait
wait
# Indexing files
date +"%d/%m/%Y %H:%M:%S" 
echo "Sorting and indexing merged BAM files."
for c in *.rg.merged.bam
do
	FILE2=$c
	FILENAME2=${c%.*}
    $SAMTOOLS view -bS $FILE2 | $SAMTOOLS sort -@ $THREADS - -o $FILENAME2.sorted.bam
    $SAMTOOLS index -@ $THREADS $FILENAME2.sorted.bam $FILENAME2.sorted.bai
done
#wait
wait
date +"%d/%m/%Y %H:%M:%S" 
echo "Finnished sorting and indexing merged BAM files."

############################################################################################
### Copy data from scratch back to home dir and clean scratch
mkdir -p $OUTPUT_DIR
cp -avr $SCRATCHDIR $OUTPUT_DIR  || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"