#!/bin/bash
#PBS -l walltime=96:0:0 
#PBS -q uv@wagap-pro.cerit-sc.cz 
#PBS -l select=1:ncpus=1:mem=100gb:scratch_local=100gb
#PBS -j oe
#PBS -N 05v3_METTL16_UV1_Collapsing_PCR_duplicates
#
## initialize the required application
module add star-2.5.2b
module add samtools-1.6
module add python27-modules-intel #required by multiQC
module add bedtools-2.26.0
module add ctk-1.0.7
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
INPUT_DIR="/storage/brno3-cerit/home/tskalicky/METTL16/mapping/METTL16_UV1/results/alignment/genome"
MY_GENOME="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
MY_GTF="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/Homo_sapiens.GRCh38.91.gtf.gz"
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/METTL16/mapping/METTL16_UV1/results/alignment/genome/PCR_collapsed"

MY_RAM=25 # Max RAM memory for Samtools sort
THREADS=$PBS_NUM_PPN
APPENDIX=".out.bam"
APPENDIX2=".md.sam"
APPENDIX3=".tag.bed"

# Binaries
STAR=$(which STAR)
SAMTOOLS=$(which samtools)
MULTIQC=$(which multiqc)
BEDTOOLS=$(which bedtools)
TAG2COLLAPSE=$(which tag2collapse.pl)
PARSEALIGN=$(which parseAlignment.pl)

# Check the tools versions
which $STAR
which $SAMTOOLS
which $MULTIQC
which $BEDTOOLS
which $TAG2COLLAPSE
which $PARSEALIGN

####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
cd $INPUT_DIR
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
cp -av $INPUT_DIR/*$APPENDIX $MY_GENOME $SCRATCHDIR
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR
echo "Following files were copied to scratch:"
ls -c1

### Commands
## Unpack genomes for samtools
GENOME_NAME=$(basename $MY_GENOME)
#
echo "Unpacking genome $GENOME_NAME"
unpigz -p $THREADS $GENOME_NAME

GENOME=${GENOME_NAME%.*}

# Set MAX ram
RAM=$[$MY_RAM-1] # lower it by one to leave some space
MEM_LIMIT=$RAM # $[$RAM/$THREADS] # In case we need GB per thread
## Converting BAM to SAM files that are needed for PCR collapsing of duplicates using CTK tools
mkdir -p $SCRATCHDIR/{SAM,BED}
#for i in *$APPENDIX
#do
#	READ_FOR=$i
#	SAMPLE=${i%.*.*}
#	echo "Now I am converting BAM file $READ_FOR to SAM"
#	$SAMTOOLS view -h $READ_FOR > $SCRATCHDIR/SAM/$SAMPLE".sam" &
#done
## Wait for all jobs to finish before exiting the job submission script
#wait
#echo "Done convering ALL BAM to SAM files"
#
## The parsing script relies on MD tags, which is an optional field without strict definition in SAM file format specification. 
## Some aligners might have slightly different format how they report mismatches. If other aligners than bwa is used,
## one should run the samtools view, sort and fillmd before running the parsealign.pl script
#cd $SCRATCHDIR
for c in *$APPENDIX
do
	READ_BAM=$c
	SAMPLENAME=${c%.*}
	echo "Now I am filling MD tags to SAM file $READ_BAM"
	## next line needed only if we have SAM files processed with CTK pipeline. samtools view will convert SAM to BAM and sort will sort them.
	#samtools view -bS $READ_SAM | samtools sort - $SAMPLENAME.bam.sorted 
	## fillmd will fill MD tag for visualising mismatches and insertions in an alignment of a read to a reference genome and convert file to SAM
	samtools fillmd $READ_BAM $SCRATCHDIR/$GENOME > $SCRATCHDIR/SAM/$SAMPLENAME.md.sam &
done
# Wait for all jobs to finish before exiting the job submission script
wait
echo "Done filling MD tags to ALL SAM files."
#
cd $SCRATCHDIR/SAM
for a in *$APPENDIX2
do
	READ_FOR2=$a
	SAMPLENAME2=${a%.*}
	echo "Now I am parsing SAM file $READ_FOR2 for unique mapppings only"
	# This will keep only unique mappings (with MAPQ >=1) and a minimal mapping size of 18 nt.
	perl $PARSEALIGN -v --map-qual 1 --min-len 18 --mutation-file $SCRATCHDIR/BED/$SAMPLENAME2".mutation.txt" \
	$READ_FOR2 $SCRATCHDIR/BED/$SAMPLENAME2".tag.bed" &
done
# Wait for all jobs to finish before exiting the job submission script
wait
echo "Done parsing ALL SAM files for unique mappings only"
cd $SCRATCHDIR/BED
# Keep track what proportion of reads can be mapped uniquely.
wc -l $SCRATCHDIR/BED/*.tag.bed > Uniquely_mapped_read_proportions.txt
#
for b in *$APPENDIX3
do
	READ_FOR3=$b
	SAMPLENAME3=${b%.*.*}
	echo "Now I am collapsing PCR duplicates for $READ_FOR3 file"
	# It is critical to collapse PCR duplicates, not only for the exact duplicates collapsed above, 
	# but also for those with slight differences due to sequencing errors.
	# A model-based algorithm called EM is used to identify "sufficiently distinct" barcodes. 
	# Details of the algorithm was described in the following paper: https://doi.org/10.1016/j.cell.2011.06.013
	# EM options applicable only for libraries generated with random BARCODES!
	perl $TAG2COLLAPSE -v -big --random-barcode -EM 30 --seq-error-model alignment \
	-weight --weight-in-name --keep-max-score --keep-tag-name $READ_FOR3 $SAMPLENAME3".tag.uniq.bed" &
done
echo "Done collapsing PCR duplicates for ALL files"
############################################################################################
### Copy data from scratch back to home dir and clean scratch
mkdir -p $OUTPUT_DIR
rm $SCRATCHDIR/*$APPENDIX $SCRATCHDIR/$GENOME
cp -avr $SCRATCHDIR $OUTPUT_DIR || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"







