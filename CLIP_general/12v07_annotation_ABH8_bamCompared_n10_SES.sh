#!/bin/bash
#PBS -l walltime=24:0:0 
#PBS -q uv@wagap-pro.cerit-sc.cz
#PBS -l select=1:ncpus=2:mem=50gb:scratch_local=70gb
#PBS -j oe
#PBS -N 12v06_annotation_ABH8_bamcompared
#
## initialize the required application
# module add gcc-4.9.2 # required libraries for bedtools on Zewura
module add bedtools-2.25.0 # version 2.26.0 is BROKEN!
# module add bamtools
############################################################################################
# CLIP Seq annotation pipeline
# Usage of previous bed intersect results as a database to probe results from Deeptools bamCompare
# comparison od CLIPseq and RNAseq to get annotation and enrichment of observed CLIPseq mapped reads.
# !!! If you did not used --samFlagExclude 16 (excludes RV strand reads) and --samFlagInclude 16 (work only on RV strand reads)
# in bamCompare than you have to work WITHOUT checking the STRAND !!!
#
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
#
# variables
# ABH8_CLIP_COMPARED="/storage/brno3-cerit/home/tskalicky/ABH8/bamCompare/Galaxy_ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_subtract_nosmooth.bed"
ABH8_CLIP_COMPARED="/storage/brno3-cerit/home/tskalicky/ABH8/bamCompare/Intersect_ABH8_bamcompared_and_PIPE-CLIP_enrichedClusters_n50.shuf2.sm.bed"
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/ABH8/anotace/bamcompare_annot"
ANNOTATION_DB="/storage/brno3-cerit/home/tskalicky/ABH8/anotace/final/ABH8_CLIP1-3/renamed_merged"

####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
# cp -av $ABH8/*.bam.gz $ABH8_POOLED/*.bam.gz $ABH8_PEAKS/*.bed $SCRATCHDIR
cp -av $ABH8_CLIP_COMPARED $SCRATCHDIR
mkdir $SCRATCHDIR/annotation_db
cp -av $ANNOTATION_DB/temp*.bed $SCRATCHDIR/annotation_db
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR
echo "Following files/folders were copied to scratch:"
ls -c1

### Commands
####################################################################################################
# BAM to BED conversion
# for a in *.bam
# do
# 	FILE=$a
#   	FILENAME=${a%.*}
#   	# For spliced alignments (like STAR RNA) need to create BED12 format
#   	# otherwise not only spliced reads on exons but the whole range will be 
#   	# reported as interval !!!
#   	date +"%d/%m/%Y %H:%M:%S" 
#   	echo "Converting BAM file $FILE to BED format."
#   	bamToBed -bed12 -i $FILE | sortBed -i stdin > $FILENAME.bed
#   	date +"%d/%m/%Y %H:%M:%S" 
#   	echo "Finnished converting BAM file $FILE to BED format."
# done
# wait

# for a in *.bam.gz
# do
#   FILE=$a
#   FILENAME=${a%.*}
#   # For spliced alignments (like STAR RNA) need to create BED12 format
#   # otherwise not only spliced reads on exons but the whole range will be 
#   # reported as interval !!!
#   date +"%d/%m/%Y %H:%M:%S" 
#   echo "Converting BAM file $FILE to BED format."
#   gzip -dc $FILE | bamToBed -bed12 -i stdin | sortBed -i stdin > $FILENAME.sorted.bed
#   date +"%d/%m/%Y %H:%M:%S" 
#   echo "Finnished converting BAM file $FILE to BED format."
# done
# wait
shopt nullglob # If set, Bash allows filename patterns which match no files to expand to a null string, rather than themselves.
# declare an array variables
declare -a SAMPLES=(*.bed)
declare -a ANNOT_DB=("temp_1_rRNA.ABH8_CLIP1-3.rsm.bed" "temp_2_tRNA.ABH8_CLIP1-3.rsm.bed" "temp_3_snoRNA.ABH8_CLIP1-3.rsm.bed" \
"temp_4_snRNA.ABH8_CLIP1-3.rsm.bed" "temp_5_miRNA.ABH8_CLIP1-3.rsm.bed" "temp_6_mRNA_exons.ABH8_CLIP1-3.rsm.bed" \
"temp_7_mRNA_introns.ABH8_CLIP1-3.rsm.bed" "temp_8_lincRNA.ABH8_CLIP1-3.rsm.bed" "temp_9_repetitive_elements.ABH8_CLIP1-3.rsm.bed" \
"temp_10_antisense_RNA.ABH8_CLIP1-3.rsm.bed" "temp_11_intergenic.ABH8_CLIP1-3.rsm.bed" "temp_11_NOT_intergenic.ABH8_CLIP1-3.rsm.bed")

# get length of an array
smp_num="${#SAMPLES[@]}"
annot_num="${#ANNOT_DB[@]}"
#
echo "There are $smp_num samples that will be processed."
echo "Sample names are: ${SAMPLES[@]}"
#
# use for loop to read all values and indexes
for (( a=0;a<${smp_num};a++ )); # I know, I know programmers COUNT from zero :-D
do
	SAMPLE="${SAMPLES[$a]}"
	SAMPLENAME="${SAMPLE%.*}"
	mkdir $SCRATCHDIR/$SAMPLENAME
	#
	date +"%d/%m/%Y %H:%M:%S"
	echo "Start intersecting of sample $SAMPLE with databases."
	## rRNA ##
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Start intersecting of sample $SAMPLE with database ${ANNOT_DB[0]}"
	intersectBed -wo -s -f 0.5 -split \
	-a "$SCRATCHDIR/$SAMPLE" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[0]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_1_rRNA.bed"
	#
	intersectBed -v -split \
	-a "$SCRATCHDIR/$SAMPLE" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[0]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_1_NOT.bed"
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[0]}"
	## tRNA ##
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Start intersecting of sample $SAMPLE with database ${ANNOT_DB[1]}"
	intersectBed -wo -s -f 0.5 -split \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_1_NOT.bed" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[1]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_2_tRNA.bed"
	#
	intersectBed -v -split \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_1_NOT.bed" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[1]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_2_NOT.bed"
	rm $SCRATCHDIR/$SAMPLENAME/temp_1_NOT.bed
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[1]}"
	## snoRNA
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Start intersecting of sample $SAMPLE with database ${ANNOT_DB[2]}"
	intersectBed -wo -s -f 0.5 -split \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_2_NOT.bed" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[2]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_3_snoRNA.bed"
	#
	intersectBed -v -split \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_2_NOT.bed" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[2]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_3_NOT.bed"
	rm $SCRATCHDIR/$SAMPLENAME/temp_2_NOT.bed
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[2]}"
	## snRNA
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Start intersecting of sample $SAMPLE with database ${ANNOT_DB[3]}"
	intersectBed -wo -s -f 0.5 -split \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_3_NOT.bed" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[3]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_4_snRNA.bed"
	#
	intersectBed -v -split \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_3_NOT.bed" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[3]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_4_NOT.bed"
	rm $SCRATCHDIR/$SAMPLENAME/temp_3_NOT.bed
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[3]}"
	## miRNA
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Start intersecting of sample $SAMPLE with database ${ANNOT_DB[4]}"
	intersectBed -wo -s -f 0.5 -split \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_4_NOT.bed" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[4]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_5_miRNA.bed"
	#
	intersectBed -v -split \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_4_NOT.bed" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[4]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_5_NOT.bed"
	rm $SCRATCHDIR/$SAMPLENAME/temp_4_NOT.bed
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[4]}"
	## exons mRNA
	# In this version I used MERGED EXON Database
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Start intersecting of sample $SAMPLE with database ${ANNOT_DB[5]}"
	intersectBed -wo -s -f 0.5 -split \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_5_NOT.bed" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[5]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_6_mRNA_exons.bed"
	#
	intersectBed -v -split \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_5_NOT.bed" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[5]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_6_NOT.bed"
	rm $SCRATCHDIR/$SAMPLENAME/temp_5_NOT.bed
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[5]}"
	## introns mRNA
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Start intersecting of sample $SAMPLE with database ${ANNOT_DB[6]}"
	intersectBed -wo -s -f 0.5 -split \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_6_NOT.bed" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[6]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_7_mRNA_introns.bed"
	#
	intersectBed -v -split \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_6_NOT.bed" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[6]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_7_NOT.bed"
	rm $SCRATCHDIR/$SAMPLENAME/temp_6_NOT.bed
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[6]}"
	## lincRNA
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Start intersecting of sample $SAMPLE with database ${ANNOT_DB[7]}"
	intersectBed -wo -s -f 0.5 -split \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_7_NOT.bed" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[7]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_8_lincRNA.bed"
	#
	intersectBed -v -split \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_7_NOT.bed" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[7]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_8_NOT.bed"
	rm $SCRATCHDIR/$SAMPLENAME/temp_7_NOT.bed
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[7]}"
	## repetitive elements
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Start intersecting of sample $SAMPLE with database ${ANNOT_DB[8]}"
	intersectBed -wo -s -f 0.5 -split \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_8_NOT.bed" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[8]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_9_repetitive_elements.bed"
	#
	intersectBed -v -split \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_8_NOT.bed" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[8]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_9_NOT.bed"
	rm $SCRATCHDIR/$SAMPLENAME/temp_8_NOT.bed
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[8]}"
	## antisense RNA
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Start intersecting of sample $SAMPLE with database ${ANNOT_DB[9]}"
	intersectBed -wo -s -f 0.5 -split \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_9_NOT.bed" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[9]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_10_antisense_RNA.bed"
	#
	intersectBed -v -split \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_9_NOT.bed" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[9]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_10_NOT.bed"
	rm $SCRATCHDIR/$SAMPLENAME/temp_9_NOT.bed
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[9]}"
	## intergenic region
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Start intersecting of sample $SAMPLE with database ${ANNOT_DB[10]}"
	intersectBed -wo -s -f 0.5 -split \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_10_NOT.bed" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[10]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_11_intergenic.bed"
	#
	intersectBed -v -split \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_10_NOT.bed" \
	-b "$SCRATCHDIR/annotation_db/${ANNOT_DB[10]}" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_11_NOT_intergenic.bed"
	rm $SCRATCHDIR/$SAMPLENAME/temp_10_NOT.bed
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[10]}"
	## Merge all annotated exone genome elements.
	mkdir $SCRATCHDIR/$SAMPLENAME/renamed_merged
	REN_MERGED="$SCRATCHDIR/$SAMPLENAME/renamed_merged"
	cd $SCRATCHDIR/$SAMPLENAME
	for c in temp*.bed
	do
		SAMPLE3=$c
		SAMPLENAME3=${c%.*}
		date +"%d/%m/%Y %H:%M:%S"
		echo "Reshuffling fields in bed file $SAMPLE3"
		awk -F'[\t]' -v OFS='\t' '{print $1, $2, $3, $8, $4, $10}' "$SCRATCHDIR/$SAMPLENAME/$SAMPLE3" > "$REN_MERGED/$SAMPLENAME3".ren.bed
		date +"%d/%m/%Y %H:%M:%S"
		echo "Finnished Reshuffling fields in bed file $SAMPLE3"
		#
		date +"%d/%m/%Y %H:%M:%S"
		echo "Sorting and merging renamed file $SAMPLE3"
		sortBed -i "$REN_MERGED/$SAMPLENAME3".ren.bed | \
		mergeBed -c 4,5,6 -o distinct,median,distinct -i stdin > "$REN_MERGED/$SAMPLENAME3".rsm.bed
		date +"%d/%m/%Y %H:%M:%S"
		echo "Finnished Sorting and merging renamed file $SAMPLE3."
	done &
	#
	wait
	# 
	# Count number of lines for each file and create summary report
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Creating summary report."
	for b in "$REN_MERGED/"*.rsm.bed
	do
		SAMPLE2=$b
		SAMPLENAME2=${b%.*}
		wc -l $SAMPLE2
	done | sort -k1,1 -rn >> "$REN_MERGED/Summary.rsm.txt"
	## Manual count one-liner
	# for b in temp*.bed; do wc -l $b; done  | sort -k1,1 -rn >> Summary_METTL16-PAR.txt
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Summary report created."
	## Merge files
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Merging all annotated files."
	cat "$REN_MERGED/*.rsm.bed" \
	> "$REN_MERGED/$SAMPLENAME_all_annot.rsm.bed"
done
# Wait for all background jobs to finish before the script continues
wait
# gzip -f "$SCRATCHDIR/$SAMPLE/$SAMPLEname.ann.bed"
## Remove temp
# rm "$SCRATCHDIR/$smp"/temp_?_NOT.$SAMPLENAME.bed

############################################################################################
### Copy data from scratch back to home dir and clean scratch
mkdir -p $OUTPUT_DIR
rm -r $SCRATCHDIR/annotation_db $SCRATCHDIR/*.bam
cp -avr $SCRATCHDIR $OUTPUT_DIR || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
