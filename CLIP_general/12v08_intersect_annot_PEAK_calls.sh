#!/bin/bash
#PBS -l walltime=24:0:0 
#PBS -q uv@wagap-pro.cerit-sc.cz 
#PBS -l select=1:ncpus=2:mem=50gb:scratch_local=50gb
#PBS -j oe
#PBS -N 12_08_intersect_annot_all
#export PBS_SERVER=wagap-pro.cerit-sc.cz # needed only when executing from arien frontends
#
## initialize the required application
module add bedtools-2.25.0 # version 2.26.0 is BROKEN!
# module add bamtools
# module add picard-2.9.0 # will also initialize system variable $PICARD pointing into Picard Tools install dir.
# module add samtools-1.4
############################################################################################
# CLIPseq analysis
# Will take identified PEAKS from peak calling and intersect them with already annotated results from BamCompare 
# from CLIPseq vs. RNAseq
############################################################################################
# Check if tools are installed
BEDTOOLS=$(which bedtools)
SORTBED=$(which sortBed)
MERGEBED=$(which mergeBed)
COMPLEMENTBED=$(which complementBed)
INTERSECTBED=$(which intersectBed)
BAMTOBED=$(which bamToBed)
#
which $BEDTOOLS
which $SORTBED
which $MERGEBED
which $COMPLEMENTBED
which $INTERSECTBED
### Variables
# AWK renaming on already produced renamed and sorted files because Strand information in bed format has to be in 6th collumn!!!
# Else it will not be processed by BedIntersect. 
# After bedCompare and AWK field switch perform bedIntersect and only after that mergeBed or you will get wrong possitions,
# spanning over several tRNA genes! -> This is because bedCompare even with bins=10bp will produce many adjacent bins that
# don't correspond to gene intervals (consider even smaller bins? like n=5 ?)
SES_LOG2="/mnt/storage-brno3-cerit/nfs4/home/tskalicky/ABH8/bamCompare/Galaxy_bamCompare_ABH8_CLIP_vs_RNAseq_bin10_SES_log2_nosmooth.bed"
SES_SUBTRACT="/mnt/storage-brno3-cerit/nfs4/home/tskalicky/ABH8/bamCompare/Galaxy_ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_subtract_nosmooth.bed"
PEAKS="/storage/brno3-cerit/home/tskalicky/ABH8/peak_calling/PIPE-CLIP_ABH8_enrichedClusters_n50.bed"
ANNOTATED="/storage/brno3-cerit/home/tskalicky/ABH8/anotace/bamcompare_annot/Galaxy_ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_subtract_nosmooth/renamed_merged/temp_2_tRNA.rsm.bed"
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/ABH8/anotace/bamcompare_annot"
####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
cp -av $PEAKS $SES_LOG2 $SES_SUBTRACT $SCRATCHDIR
# cp -av $ANNOTATED $SCRATCHDIR
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR
echo "Following files/folders were copied to scratch:"
ls -c1

### Commands
####################################################################################################
shopt nullglob # If set, Bash allows filename patterns which match no files to expand to a null string, rather than themselves.
# declare an array variables
declare -a FW_BED=(*.FW.bed)
declare -a RV_BED=(*.RV.bed)
#I know, it is not ideal solution and ugly code but it works with minimum effort. :-D
if [[ -f ${FW_BED[0]} ]]; then
	for a in *.FW.bed; do
		FILE=$a
		FILENAME=${a%.*}
    	date +"%d/%m/%Y %H:%M:%S"
		echo "Start filling strand info to file $FILE"
		awk -F'[\t]' -v OFS='\t' '{print $1, $2, $3, $4, $4, "+"}' $FILE > "$FILENAME".fill.bed
		date +"%d/%m/%Y %H:%M:%S"
		echo "Finnished filling strand info to file $FILE"
	done
else
	date +"%d/%m/%Y %H:%M:%S"
	echo "There is no *.FW.bed file for strand info filling!" && exit 1
fi
if [[ -f ${RV_BED[0]} ]]; then
	for b in *.RV.bed; do
		FILE2=$b
		FILENAME2=${b%.*}
    	date +"%d/%m/%Y %H:%M:%S"
		echo "Start filling strand info to file $FILE2"
		awk -F'[\t]' -v OFS='\t' '{print $1, $2, $3, $4, $4, "-"}' $FILE2 > "$FILENAME2".fill.bed
		date +"%d/%m/%Y %H:%M:%S"
		echo "Finnished filling strand info to file $FILE2"
	done
else
	date +"%d/%m/%Y %H:%M:%S"
	echo "There is no *.RV.bed file for strand info filling!" && exit 1
fi
#
wait
#
for c in *.FW.fill.bed; do
		FILE3=$c
		FILENAME3=${c%.*.*.*}
		if [[ -a "$FILENAME3".RV.fill.bed ]]; then
			date +"%d/%m/%Y %H:%M:%S"
			echo "Start combining $FILE3 and "$FILENAME3".RV.fill.bed files with filled strand info"
			cat "$FILE3" "$FILENAME3".RV.fill.bed |
			sortBed -i stdin > "$FILENAME3".filled.bed
			echo "Finnished combining $FILE3 and "$FILENAME3".RV.fill.bed files with filled strand info"
		else
			echo "There are no matching FW and RV files for merging!" && exit 1
		fi
done
#wait
wait
#
declare -a SAMPLES=(*.filled.bed)
# get length of an array
smp_num="${#SAMPLES[@]}"
#
echo "There are $smp_num samples that will be cross-intersected."
echo "Sample names are: ${SAMPLES[@]}"

   SAMPLE="${SAMPLES[0]}"
SAMPLENAME="${SAMPLE%.*.*}"
SAMPLE2="${SAMPLES[1]}"
SAMPLENAME2="${SAMPLE%.*.*}"
   # echo "Renaming reads in file $FILE"
   # Strand information in bed format has to be in 6th collumn!!!
   # awk -F'[\t]' -v OFS='\t' '{print $1, $2, $3, $4, $6, $5}' $a | mergeBed -s -c 4,5,6 -o distinct,median,distinct -i stdin > $FILENAME".merged.bed"
date +"%d/%m/%Y %H:%M:%S"
echo "Start intersecting of sample"
intersectBed -wo -s -f 0.5 -split \
-a "${SAMPLES[0]}" \
-b "${SAMPLES[1]}" \
> ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_log2_and_subtract_nosmooth.bed
#
intersectBed -v -split \
-a "${SAMPLES[0]}" \
-b "${SAMPLES[1]}" \
> NOT_intersect_ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_log2_and_subtract_nosmooth.bed
#
date +"%d/%m/%Y %H:%M:%S"
echo "Reshuffling fields in bed file ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_log2_and_subtract_nosmooth.bed"
# awk -F'[\t]' -v OFS='\t' '{print $1, $2, $3, $8, $4}' ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_log2_and_subtract_nosmooth.bed \
# > ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_log2_and_subtract_nosmooth.ren.bed
awk -F'[\t]' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $11 }' ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_log2_and_subtract_nosmooth.bed \
> ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_log2_and_subtract_nosmooth.shuf.bed
date +"%d/%m/%Y %H:%M:%S"
echo "Finnished Reshuffling fields in bed file ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_log2_and_subtract_nosmooth.bed"
#
intersectBed -wo -s -f 0.5 -split \
-a "PIPE-CLIP_ABH8_enrichedClusters_n50.bed" \
-b "ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_log2_and_subtract_nosmooth.shuf.bed" \
> Intersect_ABH8_bamcompared_and_PIPE-CLIP_enrichedClusters_n50.bed
#
intersectBed -v -split \
-a "PIPE-CLIP_ABH8_enrichedClusters_n50.bed" \
-b "ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_log2_and_subtract_nosmooth.shuf.bed" \
> NOT_Intersect_ABH8_bamcompared_and_PIPE-CLIP_enrichedClusters_n50.bed
#
date +"%d/%m/%Y %H:%M:%S"
echo "Sorting and merging shuffled file Intersect_ABH8_bamcompared_and_PIPE-CLIP_enrichedClusters_n50.bed"
sortBed -i Intersect_ABH8_bamcompared_and_PIPE-CLIP_enrichedClusters_n50.bed | \
mergeBed -s -c 4,5,6,7 -o median,median,distinct,median -i stdin > Intersect_ABH8_bamcompared_and_PIPE-CLIP_enrichedClusters_n50.sm.bed
date +"%d/%m/%Y %H:%M:%S"
echo "Finnished Sorting and merging shuffled file Intersect_ABH8_bamcompared_and_PIPE-CLIP_enrichedClusters_n50.bed"
# # Need to finish this!
# intersectBed -wo -s -f 0.5 -split \
# -a "temp_2_tRNA.rsm.bed" \
# -b "Intersect_ABH8_bamcompared_and_PIPE-CLIP_enrichedClusters_n50.sm.bed" \
# > Intersect_ABH8_bamcomp_PIPE-CLIP_and_tRNA.bed
# #
# intersectBed -v -split \
# -a "temp_2_tRNA.rsm.bed" \
# -b "Intersect_ABH8_bamcompared_and_PIPE-CLIP_enrichedClusters_n50.sm.bed" \
# > NOT_intersect_ABH8_bamcomp_PIPE-CLIP_and_tRNA.bed
# #
# echo "Reshuffling fields in bed file Intersect_ABH8_bamcompared_and_PIPE-CLIP_enrichedClusters_n50.bed"
# # awk -F'[\t]' -v OFS='\t' '{print $1, $2, $3, $8, $4}' ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_log2_and_subtract_nosmooth.bed \
# # > ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_log2_and_subtract_nosmooth.ren.bed
# awk -F'[\t]' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "," $11 "\t" $6 "\t" $13}' Intersect_ABH8_bamcompared_and_PIPE-CLIP_enrichedClusters_n50.bed \
# > Intersect_ABH8_bamcompared_and_PIPE-CLIP_enrichedClusters_n50.shuf.bed
# date +"%d/%m/%Y %H:%M:%S"
# echo "Finnished Reshuffling fields in bed file Intersect_ABH8_bamcompared_and_PIPE-CLIP_enrichedClusters_n50.bed"
# wait
wait
############################################################################################
### Copy data from scratch back to home dir and clean scratch
mkdir -p $OUTPUT_DIR
cp -avr $SCRATCHDIR $OUTPUT_DIR || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"