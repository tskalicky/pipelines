#!/bin/bash
#PBS -l walltime=24:0:0 
#PBS -q uv@wagap-pro.cerit-sc.cz 
#PBS -l select=1:ncpus=2:mem=200gb:scratch_local=200gb
#PBS -j oe
#PBS -N 11v03_annotation_ABH8
#export PBS_SERVER=wagap-pro.cerit-sc.cz # needed only when executing from arien frontends#
#
## initialize the required application
module add gcc-4.9.2 # required libraries for bedtools
module add bedtools-2.26.0
module add bamtools
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
ABH8_BED="/storage/brno3-cerit/home/tskalicky/ABH8/anotace/job_768846"
ABH8="/storage/brno3-cerit/home/tskalicky/ABH8/mapping/PCR_collapsed_Picard"
ABH8_POOLED="/storage/brno3-cerit/home/tskalicky/ABH8/mapping/PCR_collapsed_Picard/pooled"
ABH8_PEAKS="/storage/brno3-cerit/home/tskalicky/ABH8/peak_calling"
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/ABH8/anotace"
ANNOTATION_DB="/storage/brno3-cerit/home/tskalicky/genomes/human/annotation/my_annot_DB/ensembl/bed"
####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
cp -av $ABH8_BED/*.bed $ABH8_PEAKS/*.bed $SCRATCHDIR
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR
echo "Following files/folders were copied to scratch:"
ls -c1

### Commands
####################################################################################################
# BAM to BED conversion
for a in *.bam
do
	FILE=$a
  	FILENAME=${a%.*}
  	# For spliced alignments (like STAR RNA) need to create BED12 format
  	# otherwise not only spliced reads on exons but the whole range will be 
  	# reported as interval !!!
  	date +"%d/%m/%Y %H:%M:%S" 
  	echo "Converting BAM file $FILE to BED format."
  	bamToBed -i $FILE | sortBed -bed12 -i stdin > $FILENAME.bed
  	date +"%d/%m/%Y %H:%M:%S" 
  	echo "Finnished converting BAM file $FILE to BED format."
done
wait

# for a in *.bam.gz
# do
#   FILE=$a
#   FILENAME=${a%.*.*}
#   gzip -dc $FILE | bamToBed -i stdin > $FILENAME.bed
# done
# wait
#
shopt -s nullglob # If set, Bash allows filename patterns which match no files to expand to a null string, rather than themselves.
# declare an array variables
declare -a SAMPLES=(*.bed)
declare -a ANNOT_DB=("GRCh38.91.rRNA.bed" "GRCh38.91-tRNAs.ren.bed" "GRCh38.91.snoRNA.bed" \
  "GRCh38.91.snRNA.bed" "GRCh38.91.miRNA.bed" "GRCh38.91.mRNA_exons.merged.bed" 
  "GRCh38.91.mRNA_introns.bed" "GRCh38.91.lincRNA.bed" "ucsc_hg38_repetitive_elements.ren.bed" \
  "GRCh38.91.antisense_mRNA.bed" "GRCh38.91.intergenic.bed")
# get length of an array
smp_num="${#SAMPLES[@]}"
annot_num="${#ANNOT_DB[@]}"
#
echo "There are $smp_num samples that will be processed."
echo "Sample names are: ${SAMPLES[@]}"
#
# use for loop to read all values and indexes
for (( a=0;a<=${smp_num};a++ )); # I know, I know programmers COUNT from zero :-D
do
	SAMPLE="${SAMPLES[$a]}"
	SAMPLENAME="${SAMPLE%.*}"
	mkdir $SCRATCHDIR/$SAMPLENAME
	#
	date +"%d/%m/%Y %H:%M:%S"
	echo "Start intersecting of sample $SAMPLE with databases."
	## rRNA ##
	intersectBed -wo -s \
	-a "$INPUT/$SAMPLE" \
	-b "$ANNOTATION_DB/GRCh38.91.rRNA.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_1_rRNA.$SAMPLENAME.bed"
	#
	intersectBed -v -s \
	-a "$INPUT/$SAMPLE" \
	-b "$ANNOTATION_DB/GRCh38.91.rRNA.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_1_NOT.$SAMPLENAME.bed"
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[0]}"
	## tRNA ##
	intersectBed -wo -s \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_1_NOT.$SAMPLENAME.bed" \
	-b "$ANNOTATION_DB/GRCh38.91-tRNAs.ren.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_2_tRNA.$SAMPLENAME.bed"
	#
	intersectBed -v -s \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_1_NOT.$SAMPLENAME.bed" \
	-b "$ANNOTATION_DB/GRCh38.91-tRNAs.ren.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_2_NOT.$SAMPLENAME.bed"
	rm $SCRATCHDIR/$SAMPLENAME/temp_1_NOT.$SAMPLENAME.bed
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[1]}"
	## snoRNA
	intersectBed -wo -s \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_2_NOT.$SAMPLENAME.bed" \
	-b "$ANNOTATION_DB/GRCh38.91.snoRNA.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_3_snoRNA.$SAMPLENAME.bed"
	#
	intersectBed -v -s \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_2_NOT.$SAMPLENAME.bed" \
	-b "$ANNOTATION_DB/GRCh38.91.snoRNA.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_3_NOT.$SAMPLENAME.bed"
	rm $SCRATCHDIR/$SAMPLENAME/temp_2_NOT.$SAMPLENAME.bed
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[2]}"
	## snRNA
	intersectBed -wo -s \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_3_NOT.$SAMPLENAME.bed" \
	-b "$ANNOTATION_DB/GRCh38.91.snRNA.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_4_snRNA.$SAMPLENAME.bed"
	#
	intersectBed -v -s \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_3_NOT.$SAMPLENAME.bed" \
	-b "$ANNOTATION_DB/GRCh38.91.snRNA.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_4_NOT.$SAMPLENAME.bed"
	rm $SCRATCHDIR/$SAMPLENAME/temp_3_NOT.$SAMPLENAME.bed
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[3]}"
	## miRNA
	intersectBed -wo -s \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_4_NOT.$SAMPLENAME.bed" \
	-b "$ANNOTATION_DB/GRCh38.91.miRNA.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_5_miRNA.$SAMPLENAME.bed"
	#
	intersectBed -v -s \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_4_NOT.$SAMPLENAME.bed" \
	-b "$ANNOTATION_DB/GRCh38.91.miRNA.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_5_NOT.$SAMPLENAME.bed"
	rm $SCRATCHDIR/$SAMPLENAME/temp_4_NOT.$SAMPLENAME.bed
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[4]}"
	## exons mRNA
	# In this version I used MERGED EXON Database
	intersectBed -wo -s \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_5_NOT.$SAMPLENAME.bed" \
	-b "$ANNOTATION_DB/GRCh38.91.mRNA_exons.merged.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_6_mRNA_exons.$SAMPLENAME.bed"
	#
	intersectBed -v -s \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_5_NOT.$SAMPLENAME.bed" \
	-b "$ANNOTATION_DB/GRCh38.91.mRNA_exons.merged.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_6_NOT.$SAMPLENAME.bed"
	rm $SCRATCHDIR/$SAMPLENAME/temp_5_NOT.$SAMPLENAME.bed
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[5]}"
	## introns mRNA
	intersectBed -wo -s \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_6_NOT.$SAMPLENAME.bed" \
	-b "$ANNOTATION_DB/GRCh38.91.mRNA_introns.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_7_mRNA_introns.$SAMPLENAME.bed"
	#
	intersectBed -v -s \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_6_NOT.$SAMPLENAME.bed" \
	-b "$ANNOTATION_DB/GRCh38.91.mRNA_introns.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_7_NOT.$SAMPLENAME.bed"
	rm $SCRATCHDIR/$SAMPLENAME/temp_6_NOT.$SAMPLENAME.bed
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[6]}"
	## lincRNA
	intersectBed -wo -s \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_7_NOT.$SAMPLENAME.bed" \
	-b "$ANNOTATION_DB/GRCh38.91.lincRNA.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_8_lincRNA.$SAMPLENAME.bed"
	#
	intersectBed -v -s \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_7_NOT.$SAMPLENAME.bed" \
	-b "$ANNOTATION_DB/GRCh38.91.lincRNA.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_8_NOT.$SAMPLENAME.bed"
	rm $SCRATCHDIR/$SAMPLENAME/temp_7_NOT.$SAMPLENAME.bed
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[7]}"
	## repetitive elements
	intersectBed -wo -s \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_8_NOT.$SAMPLENAME.bed" \
	-b "$ANNOTATION_DB/ucsc_hg38_repetitive_elements.ren.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_9_repetitive_elements.$SAMPLENAME.bed"
	#
	intersectBed -v -s \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_8_NOT.$SAMPLENAME.bed" \
	-b "$ANNOTATION_DB/ucsc_hg38_repetitive_elements.ren.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_9_NOT.$SAMPLENAME.bed"
	rm $SCRATCHDIR/$SAMPLENAME/temp_8_NOT.$SAMPLENAME.bed
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[8]}"
	## antisense RNA
	intersectBed -wo -s \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_9_NOT.$SAMPLENAME.bed" \
	-b "$ANNOTATION_DB/GRCh38.91.antisense_mRNA.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_10_antisense_RNA.$SAMPLENAME.bed"
	#
	intersectBed -v -s \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_9_NOT.$SAMPLENAME.bed" \
	-b "$ANNOTATION_DB/GRCh38.91.antisense_mRNA.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_10_NOT.$SAMPLENAME.bed"
	rm $SCRATCHDIR/$SAMPLENAME/temp_9_NOT.$SAMPLENAME.bed
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[9]}"
	## intergenic region
	intersectBed -wo \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_10_NOT.$SAMPLENAME.bed" \
	-b "$ANNOTATION_DB/GRCh38.91.intergenic.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_11_intergenic.$SAMPLENAME.bed"
	#
	intersectBed -v \
	-a "$SCRATCHDIR/$SAMPLENAME/temp_10_NOT.$SAMPLENAME.bed" \
	-b "$ANNOTATION_DB/GRCh38.91.intergenic.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/temp_11_NOT.$SAMPLENAME.bed"
	rm $SCRATCHDIR/$SAMPLENAME/temp_10_NOT.$SAMPLENAME.bed
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Finnished intersecting of sample $SAMPLE with database ${ANNOT_DB[10]}"
	## Merge all annotated exone genome elements.
	# date +"%d/%m/%Y %H:%M:%S"
	# echo "Merging overlapping exons."
	# sortBed -i $SCRATCHDIR/$SAMPLENAME/temp_6_mRNA_exons.$SAMPLENAME.bed | \
	# mergeBed -c 10,12 -o distinct,distinct -i stdin > $SCRATCHDIR/$SAMPLENAME/temp_6_mRNA_exons.merged.$SAMPLENAME.bed
	# date +"%d/%m/%Y %H:%M:%S"
	# echo "Finnished merging overlaping exons."
	# echo "Finnished intersecting all samples."
	# manual merging
	# sortBed -i temp_6_mRNA_exons.METTL16UV1-2_GRCh38.pcr_dedupl.rg.bed > temp_6_mRNA_exons.METTL16UV1-2_GRCh38.pcr_dedupl.rg.bed
	# mergeBed -c 10,12 -o distinct,distinct -i temp_6_mRNA_exons.METTL16UV1-2_GRCh38.pcr_dedupl.rg.bed > temp_6_mRNA_exons.merged.METTL16UV1-2_GRCh38.pcr_dedupl.rg.bed
	#
	# Manual Intersecting of not assigned leftover temp_11_NOT.$SAMPLENAME.bed again against exons
	# intersectBed -wo \
	# -a "/home/tomas/CEITEC_lab/ABH8/anotace/ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged/temp_11_NOT.ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.bed" \
	# -b "/home/tomas/ownCloud/CEITEC_lab/genomes/annotation/my_annot_DB/ensembl/bed/GRCh38.91.mRNA_exons.merged.bed" \
	# > "/home/tomas/CEITEC_lab/ABH8/anotace/ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged/test/temp_12_exons_ABH8_CLIP1-3_GRCh38.bed"
	# #
	# intersectBed -v \
	# -a "/home/tomas/CEITEC_lab/ABH8/anotace/ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged/test/temp_12_exons_ABH8_CLIP1-3_GRCh38.bed" \
	# -b "/home/tomas/ownCloud/CEITEC_lab/genomes/annotation/my_annot_DB/ensembl/bed/GRCh38.91.mRNA_exons.merged.bed" \
	# > "/home/tomas/CEITEC_lab/ABH8/anotace/ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged/test/temp_12_NOT_exons.merged_ABH8_CLIP1-3_GRCh38.bed"
	# 
	# Count number of lines for each file
	for b in $SCRATCHDIR/$SAMPLENAME/temp*.bed
	do
		SAMPLE2=$b
		SAMPLENAME2=${b%.*}
		wc -l $SAMPLE2
	done | sort -k1,1 -rn >> $SCRATCHDIR/$SAMPLENAME/Summary_"$SAMPLENAME".txt
	## Manual count one-liner
	# for b in temp*.bed; do wc -l $b; done  | sort -s -k1,1 -rn >> Summary_METTL16-PAR.txt
	## Merge files
	date +"%d/%m/%Y %H:%M:%S" 
	echo "Merging all annotated files."
	cat "$SCRATCHDIR/$SAMPLENAME/temp_1_rRNA.$SAMPLENAME.bed" \
	"$SCRATCHDIR/$SAMPLENAME/temp_2_tRNA.$SAMPLENAME.bed" \
	"$SCRATCHDIR/$SAMPLENAME/temp_3_snoRNA.$SAMPLENAME.bed" \
	"$SCRATCHDIR/$SAMPLENAME/temp_4_snRNA.$SAMPLENAME.bed" \
	"$SCRATCHDIR/$SAMPLENAME/temp_5_miRNA.$SAMPLENAME.bed" \
	"$SCRATCHDIR/$SAMPLENAME/temp_6_mRNA_exons.$SAMPLENAME.bed" \
	"$SCRATCHDIR/$SAMPLENAME/temp_7_mRNA_introns.$SAMPLENAME.bed" \
	"$SCRATCHDIR/$SAMPLENAME/temp_8_lincRNA.$SAMPLENAME.bed" \
	"$SCRATCHDIR/$SAMPLENAME/temp_9_repetitive_elements.$SAMPLENAME.bed" \
	"$SCRATCHDIR/$SAMPLENAME/temp_10_antisense_RNA.$SAMPLENAME.bed" \
	"$SCRATCHDIR/$SAMPLENAME/temp_11_intergenic.$SAMPLENAME.bed" \
	> "$SCRATCHDIR/$SAMPLENAME/$SAMPLENAME.all_annot.bed"
done
# Wait for all background jobs to finish before the script continues
wait
# gzip -f "$INPUT/$SAMPLE/$SAMPLEname.ann.bed"
## Remove temp
# rm "$INPUT/$smp"/temp_?_NOT.$SAMPLENAME.bed

############################################################################################
### Copy data from scratch back to home dir and clean scratch
mkdir -p $OUTPUT_DIR
rm $SCRATCHDIR/*.bam
cp -avr $SCRATCHDIR $OUTPUT_DIR || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
