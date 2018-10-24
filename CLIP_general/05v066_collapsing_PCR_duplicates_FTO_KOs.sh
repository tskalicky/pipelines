#!/bin/bash
#PBS -l walltime=48:0:0 
#PBS -q uv@wagap-pro.cerit-sc.cz 
#PBS -l select=1:ncpus=2:mem=50gb:scratch_local=150gb
#PBS -j oe
#PBS -N 05v066_FTO_KOs_Collapsing_PCR_duplicates
#
## initialize the required application
module add samtools-1.6
module add python27-modules-gcc #required by multiQC, import gcc not intel version!
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
INPUT_DIR="/storage/brno3-cerit/home/tskalicky/FTO/mapping/deduplicated/KO1/results/alignment/genome"
INPUT_DIR2="/storage/brno3-cerit/home/tskalicky/FTO/mapping/deduplicated/KO2/results/alignment/genome"
INPUT_DIR3="/storage/brno3-cerit/home/tskalicky/FTO/mapping/deduplicated/KO3/results/alignment/genome"
MY_GENOME="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
MY_GTF="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/Homo_sapiens.GRCh38.91.gtf.gz"
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/FTO/mapping/deduplicated/PCR_collapsed-genome/FTO_KO1-3"

MY_RAM=50 # Max RAM memory for Samtools sort
THREADS=$PBS_NUM_PPN
SAMPLES="FTO_KO1-3"
APPENDIX=".sortedByCoord.out.bam"
APPENDIX2=".md.sam"
APPENDIX3=".tag.bed"
APPENDIX4=".tag.uniq.bed"

# Binaries
SAMTOOLS=$(which samtools)
MULTIQC=$(which multiqc)
BEDTOOLS=$(which bedtools)
TAG2COLLAPSE=$(which tag2collapse.pl)
PARSEALIGN=$(which parseAlignment.pl)
SELECTROW=$(which selectRow.pl)
BED2RGB=$(which bed2rgb.pl)
TAG2PROFILE=$(which tag2profile.pl)
BED2ANNOTATION=$(which bed2annotation.pl)

# Check the tools versions
which $STAR
which $SAMTOOLS
which $MULTIQC
which $BEDTOOLS
which $TAG2COLLAPSE
which $PARSEALIGN
which $SELECTROW
which $BED2RGB
which $TAG2PROFILE
which $BED2ANNOTATION
####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
cd $INPUT_DIR
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
cp -av $INPUT_DIR/*$APPENDIX $INPUT_DIR2/*$APPENDIX $INPUT_DIR3/*$APPENDIX $MY_GENOME $SCRATCHDIR
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR
echo "Following files/folders were copied to scratch:"
ls -c1

### Commands
####################################################################################################
## redirecting temp files from system /var/temp which is restricted to 1GB
## until now needed only for preprocessing
# mkdir $SCRATCHDIR/temp
# export TMPDIR=$SCRATCHDIR/temp
####################################################################################################
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
# Wait for all background jobs to finish before the script continues
wait
echo "Done filling MD tags to ALL SAM files."
#
cd $SCRATCHDIR/SAM
for a in *$APPENDIX2
do
	READ_FOR2=$a
	SAMPLENAME2=${a%.*}
	echo "Now I am parsing SAM file $READ_FOR2 for unique mapppings only"
	# This will keep only unique mappings (with MAPQ >=1) and a minimal mapping size of 17 nt.
	perl $PARSEALIGN -v --map-qual 1 --min-len 17 --mutation-file $SCRATCHDIR/BED/$SAMPLENAME2".mutation.txt" \
	$READ_FOR2 $SCRATCHDIR/BED/$SAMPLENAME2".tag.bed" &
done
# Wait for all background jobs to finish before the script continues
wait
echo "Done parsing ALL SAM files for unique mappings only"
#
cd $SCRATCHDIR/BED
# Keep track what proportion of reads can be mapped uniquely.
wc -l $SCRATCHDIR/BED/*.tag.bed > Uniquely_mapped_read_proportions.txt
echo "Statistics of proportions of uniquely mapped reads stored in file:"
ls -c1 *proportions.txt
#
####################################################################################################
## Collapsing PCR duplicates
mkdir $SCRATCHDIR/BED/bedgraph
for b in *$APPENDIX3 #APPENDIX3=".tag.bed"
do
	READ_FOR3=$b
	SAMPLENAME3=${b%.*.*}
	echo "Now I am collapsing PCR duplicates for file $READ_FOR3"
	# It is critical to collapse PCR duplicates, not only for the exact duplicates collapsed above, 
	# but also for those with slight differences due to sequencing errors.
	# A model-based algorithm called EM is used to identify "sufficiently distinct" barcodes. 
	# Details of the algorithm was described in the following paper: https://doi.org/10.1016/j.cell.2011.06.013
	# EM options applicable only for libraries generated with random BARCODES!
	# If randome barcodes are present, use this EM settings:
	# perl $TAG2COLLAPSE -v -big --random-barcode -EM 30 --seq-error-model alignment \
	# -weight --weight-in-name --keep-max-score --keep-tag-name $READ_FOR3 $SAMPLENAME3".tag.uniq.bed"
	# If barcodes were not used in library prep, use this settings:
	perl $TAG2COLLAPSE -v -big -weight --weight-in-name --keep-max-score --keep-tag-name \
	$READ_FOR3 $SAMPLENAME3".tag.uniq.bed"
	echo "Done collapsing PCR duplicates for file $READ_FOR3"
	echo "Created PCR collapsed file:"
	ls -c1 *.tag.uniq.bed
	## Get the mutations in unique tags
	echo "Start calculation of mutations in unique tags for file $READ_FOR3"
	perl $SELECTROW -q 3 -f 3 $SAMPLENAME3".mutation.txt" $SAMPLENAME3".tag.uniq.bed" > $SAMPLENAME3".tag.uniq.mutation.txt"
	echo "Done calculation of mutations in unique tags for file $READ_FOR3"
	echo "Created unique tags matutation file:"
	ls -c1 *.tag.uniq.mutation.txt
	# Prepare a bedGraph file of unique tags for visualization in genome browser.
	# Visualize each individual experiment separately to make sure all experiments
	# work as expected, and to evaluate the reproducibility of biological replicates.
	echo "Preparing bedgraph visualisation of unique tags for sample $SAMPLENAME3.tag.uniq.bed"
	perl $TAG2PROFILE -v -ss -exact -of bedgraph -n Unique_Tag_Profile $SAMPLENAME3".tag.uniq.bed" $SAMPLENAME3".tag.uniq.bedgraph"
	echo "Finnished bedgraph visualisation of unique tags for sample $SAMPLENAME3.tag.uniq.bed"
	echo "Created bedgraph visualisation file of unique tags:"
	ls -c1 *.tag.uniq.bedgraph
done
# Wait for all background jobs to finish before the script continues
wait

## After getting the unique tags of each library, one might concatenate biological replicates, 
## which are distinguished by different colors. Perform ONLY if you have biological REPLICATES
echo "Checking if we have biological replicates"
shopt -s nullglob # If set, Bash allows filename patterns which match no files to expand to a null string, rather than themselves.
# declare an array variables
declare -a BEDNAMES=(*.tag.uniq.bed)
# get length of an array
quantity=${#BEDNAMES[@]}
#
if [[ $quantity -gt 1 ]]; then
	echo "Number of biological replicates is:" $quantity
	echo "Samples that will be processed as biological preplicates are: ${BEDNAMES[@]}" 
	# declare an array variables
	declare -a COLORarray=("blue" "red" "green" "brown" "cyan" "darkblue" "darkgreen" "darkred" "lightblue" "lightgreen" "pink" "purple" "skyblue" "yellow")
	# get length of an array
	length2=${#COLORarray[@]}
	# use for loop to read all values and indexes
	for (( c=0;c<${quantity};c++ )); # I know, I know programmers COUNT from zero :-D
	do
		SAMPLE=${BEDNAMES[$c]}
		COLOR=${COLORarray[$c]}
		SAMPLENAME4=${SAMPLE%.tag*}
		#
		echo "Start $COLOR color taging of biological replicates for sample: $SAMPLE"
		perl $BED2RGB -v -col $COLOR $SAMPLE $SAMPLENAME4".tag.uniq.rgb.bed" #as one example
		echo "Finnished $COLOR color taging of biological replicates for sample: $SAMPLE"
		#
		## repeat the above step for all other uniq.bed files to generate rgb.bed files, 
		## but use different rgb colors "x,x,x" (see http://www.rapidtables.com/web/color/RGB_Color.htm or other charts)
	done
	# Wait for all background jobs to finish before the script continues
	wait
	echo "Start concatenation of biological replicates"
	cat *tag.uniq.rgb.bed > $SAMPLES".pool.tag.uniq.rgb.bed"
	cat *tag.uniq.mutation.txt > $SAMPLES".pool.tag.uniq.mutation.txt" 
	echo "Created concatenated files:"
	ls -c1 *.pool.*
	echo "Preparing bedgraph visualisation of pooled tags from samples $SAMPLES"
	perl $TAG2PROFILE -v -ss -exact -of bedgraph -n Unique_Tag_Profile $SAMPLES".pool.tag.uniq.rgb.bed" $SAMPLES".pool.tag.uniq.bedgraph"
	echo "Finnished bedgraph visualisation of pooled tags from samples $SAMPLES"
	mv -v *.tag.uniq.bedgraph $SCRATCHDIR/BED/bedgraph/
	echo "All Bedgraph files from separate and pooled samples were moved to a folder $SCRATCHDIR/BED/bedgraph/"
	# get genomic distribution of CLIP tags
	# Check the summary file (<tag.uniq.annot.summary.txt>) for the percentage of tags mapped to CDS, 3'UTR, introns, etc.
	echo "Now I am getting genomic distribution of CLIP tags for pooled samples $SAMPLES"
	perl $BED2ANNOTATION -dbkey hg19 -ss -big -region -v -summary $SAMPLES".pool.tag.uniq.annot.summary.txt" \
	$SAMPLES".pool.tag.uniq.rgb.bed" $SAMPLES".pool.tag.uniq.annot.txt"
	echo "Done getting genomic distribution of CLIP tags for pooled samples $SAMPLES"
	echo "Created summary files:"
	ls -c1 *.annot.*
	#
	for d in *$APPENDIX4
	do
		READ_FOR5=$d
		SAMPLENAME5=${d%.tag.*.*}
		# As a diagnostic step, get the length distribution of unique tags, 
		# which should be a more faithful representation of the library:
		echo "Now I am getting the length distribution of unique tags for sample:" $READ_FOR5
		awk '{print $3-$2}' $READ_FOR5 | sort -n | uniq -c | awk '{print $2"\t"$1}' > $SAMPLENAME5".uniq.len.dist.txt"
		echo "Done getting the length distribution of unique tags for sample:" $READ_FOR5
		echo "Created length distribution of unique tags files:"
		ls -c1 *.uniq.len.dist.txt
	done
	# Wait for all background jobs to finish before the script continues
	wait

else
    echo "There is $quantity biological replicate! Skipping concatenation and color taging"
    #
	# get genomic distribution of CLIP tags
	# Check the summary file (<tag.uniq.annot.summary.txt>) for the percentage of tags mapped to CDS, 3'UTR, introns, etc.
	for e in *$APPENDIX3 #APPENDIX3=".tag.bed"
	do
		READ_FOR6=$e
		SAMPLENAME6=${e%.*.*}
		echo "Now I am getting genomic distribution of CLIP tags for sample $READ_FOR6"
		perl $BED2ANNOTATION -dbkey hg19 -ss -big -region -v -summary $SAMPLENAME6".tag.uniq.annot.summary.txt" \
		$SAMPLENAME6".tag.uniq.bed" $SAMPLENAME6".tag.uniq.annot.txt"
		echo "Done getting genomic distribution of CLIP tags for sample $READ_FOR3"
		echo "Created summary files:"
		ls -c1 *.annot.*
	done
	# Wait for all background jobs to finish before the script continues
	wait
	for d in *$APPENDIX4
	do
		READ_FOR7=$d
		SAMPLENAME7=${d%.tag.*.*}
		# As a diagnostic step, get the length distribution of unique tags, 
		# which should be a more faithful representation of the library:
		echo "Now I am getting the length distribution of unique tags for sample:" $READ_FOR7
		awk '{print $3-$2}' $READ_FOR7 | sort -n | uniq -c | awk '{print $2"\t"$1}' > $SAMPLENAME7".uniq.len.dist.txt"
		echo "Done getting the length distribution of unique tags for sample:" $READ_FOR7
		echo "Created length distribution of unique tags files:"
		ls -c1 *.uniq.len.dist.txt
	done
	# Wait for all background jobs to finish before the script continues
	wait
fi
############################################################################################
### Copy data from scratch back to home dir and clean scratch
mkdir -p $OUTPUT_DIR
rm $SCRATCHDIR/*$APPENDIX $SCRATCHDIR/$GENOME
cp -avr $SCRATCHDIR $OUTPUT_DIR || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"







