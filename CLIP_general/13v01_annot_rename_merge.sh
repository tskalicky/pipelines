#!/bin/bash
#PBS -l walltime=24:0:0 
#PBS -q default@wagap-pro.cerit-sc.cz
#PBS -l select=1:ncpus=2:mem=150gb:scratch_local=150gb
#PBS -j oe
#PBS -N 12v01_ren_merge_ABH8
#
## initialize the required application
# module add gcc-4.9.2 # required libraries for bedtools on Zewura
module add bedtools-2.25.0 # version 2.26.0 is BROKEN!
# module add bamtools
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
ABH8_CLIP_ANNOT="/storage/brno3-cerit/home/tskalicky/ABH8/anotace/job_774301/ABH8_CLIP1-3_merged_sorted"
ABH8_Piranha_ANNOT1="/storage/brno3-cerit/home/tskalicky/ABH8/anotace/job_774301/Piranha_ABH8_CLIP_covariety_RNAseq_bin-20_p-value-0.05"
ABH8_Piranha_ANNOT2="/storage/brno3-cerit/home/tskalicky/ABH8/anotace/job_774301/Piranha_ABH8_CLIP_covariety_RNAseq_bin-50_p-value-0.05"
ABH8_PIPECLIP_ANNOT="/storage/brno3-cerit/home/tskalicky/ABH8/anotace/job_774301/PIPE-CLIP_ABH8_CLIP_enrichedClusters_n20"
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/ABH8/anotace/job_774301/renamed_merged"

####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
# cp -av $ABH8/*.bam.gz $ABH8_POOLED/*.bam.gz $ABH8_PEAKS/*.bed $SCRATCHDIR
cp -avr $ABH8_CLIP_ANNOT $ABH8_Piranha_ANNOT1 $ABH8_Piranha_ANNOT2 $ABH8_PIPECLIP_ANNOT $SCRATCHDIR
mkdir -p $SCRATCHDIR/{annotation_db,renamed_merged/{ABH8_CLIP1-3_merged_sorted,Piranha_ABH8_CLIP_covariety_RNAseq_bin-20,Piranha_ABH8_CLIP_covariety_RNAseq_bin-50,PIPE-CLIP_ABH8_CLIP_enrichedClusters_n20}}
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR
echo "Following files/folders were copied to scratch:"
ls -c1

### Commands
####################################################################################################

## rename all reads to contain name of the annotated gene
## Peak calling softwares have different structure of BED files!!!
# to execute on all copied folders use:
# dirs=($(find . -type d))
# for dir in "${dirs[@]}"
# do
#	cd "$dir"
#   echo "Now processing files in folder $dir"
# done
cd $SCRATCHDIR/ABH8_CLIP1-3_merged_sorted
for a in *.bed
do
 	FILE=$a
 	FILENAME=${a%.*}
 	date +"%d/%m/%Y %H:%M:%S"
 	echo "Renaming reads in file $FILE"
 	awk -F'[\t]' -v OFS='\t' '{print $1, $2, $3, $16, $5, $6}' $FILE > "$SCRATCHDIR"/renamed_merged/ABH8_CLIP1-3_merged_sorted/"$FILENAME"_ren.bed
 	date +"%d/%m/%Y %H:%M:%S"
 	echo "Finnished renaming reads in file $FILE"
	## Merge all annotated exone genome elements in order to get intronic regions only
	# keep gene names and their strand orientation
	# -s = force merging only on the same strand
	date +"%d/%m/%Y %H:%M:%S"
 	echo "Sorting and merging reads in file $FILE"
	sortBed -i $SCRATCHDIR/renamed_merged/ABH8_CLIP1-3_merged_sorted/$FILENAME_ren.bed | \
	mergeBed -s -c 4,5,6 -o distinct,distinct,distinct -i stdin > $SCRATCHDIR/renamed_merged/ABH8_CLIP1-3_merged_sorted/$FILENAME_ren.merged.sorted.bed
	date +"%d/%m/%Y %H:%M:%S"
 	echo "Finnished sorting and merging reads in file $FILE"
done
##
cd $SCRATCHDIR/Piranha_ABH8_CLIP_covariety_RNAseq_bin-20_p-value-0.05
for a in *.bed
do
 	FILE=$a
 	FILENAME=${a%.*}
 	date +"%d/%m/%Y %H:%M:%S"
 	echo "Renaming reads in file $FILE"
 	awk -F'[\t]' '{print $1, $2, $3, $12, $5, $6}' $FILE > $SCRATCHDIR/renamed_merged/Piranha_ABH8_CLIP_covariety_RNAseq_bin-20/$FILENAME_ren.bed
 	date +"%d/%m/%Y %H:%M:%S"
 	echo "Finnished renaming reads in file $FILE"
	## Merge all annotated exone genome elements in order to get intronic regions only
	# keep gene names and their strand orientation
	# -s = force merging only on the same strand
	date +"%d/%m/%Y %H:%M:%S"
 	echo "Sorting and merging reads in file $FILE"
	sortBed -i $SCRATCHDIR/renamed_merged/Piranha_ABH8_CLIP_covariety_RNAseq_bin-20/$FILENAME_ren.bed | \
	mergeBed -s -c 4,5,6 -o distinct,distinct,distinct -i stdin > $SCRATCHDIR/renamed_merged/Piranha_ABH8_CLIP_covariety_RNAseq_bin-20/$FILENAME_ren.merged.sorted.bed
	date +"%d/%m/%Y %H:%M:%S"
 	echo "Finnished sorting and merging reads in file $FILE"
done
wait
##
cd $SCRATCHDIR/Piranha_ABH8_CLIP_covariety_RNAseq_bin-50_p-value-0.05
for a in *.bed
do
 	FILE=$a
 	FILENAME=${a%.*}
 	date +"%d/%m/%Y %H:%M:%S"
 	echo "Renaming reads in file $FILE"
 	awk -F'[\t]' '{print $1, $2, $3, $12, $5, $6}' $FILE > $SCRATCHDIR/renamed_merged/Piranha_ABH8_CLIP_covariety_RNAseq_bin-50/$FILENAME_ren.bed
 	date +"%d/%m/%Y %H:%M:%S"
 	echo "Finnished renaming reads in file $FILE"
	## Merge all annotated exone genome elements in order to get intronic regions only
	# keep gene names and their strand orientation
	# -s = force merging only on the same strand
	date +"%d/%m/%Y %H:%M:%S"
 	echo "Sorting and merging reads in file $FILE"
	sortBed -i $SCRATCHDIR/renamed_merged/Piranha_ABH8_CLIP_covariety_RNAseq_bin-50/$FILENAME_ren.bed | \
	mergeBed -s -c 4,5,6 -o distinct,distinct,distinct -i stdin > $SCRATCHDIR/renamed_merged/Piranha_ABH8_CLIP_covariety_RNAseq_bin-50/$FILENAME_ren.merged.sorted.bed
	date +"%d/%m/%Y %H:%M:%S"
 	echo "Finnished sorting and merging reads in file $FILE"
done
wait
##
cd $SCRATCHDIR/PIPE-CLIP_ABH8_CLIP_enrichedClusters_n20
for a in *.bed
do
 	FILE=$a
 	FILENAME=${a%.*}
 	date +"%d/%m/%Y %H:%M:%S"
 	echo "Renaming reads in file $FILE"
 	awk -F'[\t]' '{print $1, $2, $3, $12, $5, $6}' $FILE > $SCRATCHDIR/renamed_merged/PIPE-CLIP_ABH8_CLIP_enrichedClusters_n20/$FILENAME_ren.bed
 	date +"%d/%m/%Y %H:%M:%S"
 	echo "Finnished renaming reads in file $FILE"
	## Merge all annotated exone genome elements in order to get intronic regions only
	# keep gene names and their strand orientation
	# -s = force merging only on the same strand
	date +"%d/%m/%Y %H:%M:%S"
 	echo "Sorting and merging reads in file $FILE"
	sortBed -i $SCRATCHDIR/renamed_merged/PIPE-CLIP_ABH8_CLIP_enrichedClusters_n20/$FILENAME_ren.bed | \
	mergeBed -s -c 4,5,6 -o distinct,distinct,distinct -i stdin > $SCRATCHDIR/renamed_merged/PIPE-CLIP_ABH8_CLIP_enrichedClusters_n20/$FILENAME_ren.merged.sorted.bed
	date +"%d/%m/%Y %H:%M:%S"
 	echo "Finnished sorting and merging reads in file $FILE"
done
wait
#
############################################################################################
### Copy data from scratch back to home dir and clean scratch
mkdir -p $OUTPUT_DIR
# rm -r $SCRATCHDIR/ABH8_CLIP1-3_merged_sorted $SCRATCHDIR/Piranha_ABH8_CLIP_covariety_RNAseq_bin-20_p-value-0.05 \
# $SCRATCHDIR/Piranha_ABH8_CLIP_covariety_RNAseq_bin-50_p-value-0.05 $SCRATCHDIR/PIPE-CLIP_ABH8_CLIP_enrichedClusters_n20
cp -avr $SCRATCHDIR $OUTPUT_DIR || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
