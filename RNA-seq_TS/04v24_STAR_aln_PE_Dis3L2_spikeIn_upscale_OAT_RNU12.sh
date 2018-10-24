#!/bin/bash
#PBS -l walltime=96:0:0 
#PBS -q uv@wagap-pro.cerit-sc.cz
#PBS -l select=1:ncpus=30:mem=1500gb:scratch_local=450gb
#PBS -j oe
#PBS -N STAR_aln_Dis3L2_SpikeIn_upscape
#
# For 1 library you need cca 35GB RAM
#
## initialize the required application
module add samtools-1.4 # Samtools v1.6 are broken! 
module add python27-modules-gcc #required by multiQC
module add star-2.5.2b # METAcentrum has only old v2.5.2b which has problems with 2nd pass mapping
# module add ctk-1.0.7 # Reads were already deduplicated in failed run
#
# Alignment - PE RNA-Seq
# Pipeline to align RNA-Seq experiment as recommended in GENCODE project https://github.com/ENCODE-DCC/long-rna-seq-pipeline/tree/master/dnanexus/align-star-pe; https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/dnanexus/Readme.md with few moditifications: --twopassMode Basic; --outSAMattributes All --outFilterMismatchNoverReadLmax 0.05
# To increase sensitivity you might try to add --seedSearchStartLmax 30
#
# Requires STAR, pigz, samtools, multiQC, bedGraphToBigWig
#
# Note: do not use modified GTF (added features) to alignment as it can cause issues later on
############################################################################################
### Variables
INPUT_DIR="/storage/brno3-cerit/home/tskalicky/Dis3L2/Dasa_spikein/new_upscale/data/trimmed"
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/Dis3L2/Dasa_spikein/OAT_RNU12_upscale_NO_ANNOT_mapping"
OUTPUT_DIR=${OUTPUT_DIR}/results
OUTPUT_DIR_QC=${OUTPUT_DIR}/qc

STRANDED="false" # "true" or "false" - STAR parameters: strandedness, affects bedGraph (wiggle) files and XS tag in BAM 

APPENDIX1="_1.fastq"
APPENDIX2="_2.fastq"
APPENDIX=".fastq"
APPENDIX3=".fastq"
APPENDIX4=".dedupl.fastq"
APPENDIX5=".dedupl.fastq"
APPENDIX6="_1.dedupl.fastq"
APPENDIX7="_2.dedupl.fastq"

NO_MISMATCHES=999 # Number of mismatches per read - since we do not trim; 14 is recommended in QuantSeq FWD protocol for SE and 16 for PE; to turn this off set it to 999
PER_MISMATCHES=0.04 # Percent of mismatches per read; I usually keep this at 0.05 or 0.1; for untrimmed QuantSeq FWD it might be 0.3 -> for 50 bp single end it gives max. 15 mismatches
MAX_INTRON=1000000 # Default used by ENCODE is 1000000; to turn this off set 1
MAX_MATE_DIST=1000000 # Default used by ENCODE is 1000000; For "normal" fragments 1000 should be enough but for special cases, like chimeric we should increase this

# Read length for sjdbOverhang; --sjdbOverhang 100 should work fine for most of the data, but more specific setting based on the real read length should be more sensitive https://groups.google.com/forum/#!msg/rna-star/h9oh10UlvhI/BfSPGivUHmsJ
RD_LENGTH=100 # Read length from the sequencing. Illumina sometimes reports N+1 http://seqanswers.com/forums/archive/index.php/t-31154.html; in case you change this value uncomment next line as well
#RD_LENGTH=$[$RD_LENGTH-$TRIM_LEFT-$TRIM_RIGHT-1] # Modified read length for the index creation and mapping - should be read length -1

# Genome and annotation
GENOME_INDEX="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/STAR_index/OAT_RNU12_RD_LENGTH_100_STAR2.5.2" # Genome index with --sjdbOverhang 79 - Only for Dis3L2 RNA-seq. length=125
MY_GENOME="/storage/brno3-cerit/home/tskalicky/genomes/human/OAT_U12/OAT_U12.fa.gz"
MY_GTF="/storage/brno3-cerit/home/tskalicky/genomes/human/OAT_U12/OAT_U12_mapping.gtf.gz"


MY_RAM=1500 # Max RAM memory for Samtools sort
THREADS=$PBS_NUM_PPN
# Adding RSEM and STAR binaries into the PATH
# Will work ONLY on NFS4 connected servers
# export PATH="/storage/brno3-cerit/home/tskalicky/tools/RSEM-1.3.0/bin:$PATH"
# export PATH="/storage/brno3-cerit/home/tskalicky/anaconda2/bin:$PATH" # this version is broken

# Binaries
STAR=$(which STAR)
$STAR --version
SAMTOOLS=$(which samtools)
BEDGRAPH2BIGWIG=/storage/brno3-cerit/home/tskalicky/tools/ucsc/bedGraphToBigWig
MULTIQC=$(which multiqc)
# FASTQ2COLLAPSE=$(which fastq2collapse.pl)

# Check the tools versions
which $STAR
which $SAMTOOLS
which $BEDGRAPH2BIGWIG
which $MULTIQC
# which $FASTQ2COLLAPSE

####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
cd $INPUT_DIR
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
cp -av $MY_GENOME $MY_GTF $SCRATCHDIR
find "$INPUT_DIR" -maxdepth 2 -name "*.fastq.gz" -exec cp -vt "$SCRATCHDIR" {} + # this is not working at every occasion!!
cp -avr $GENOME_INDEX $SCRATCHDIR/ 
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR
echo "Following files were copied to scratch:"
ls -Rc1
####################################################################################################
# ## redirecting temp files from system /var/temp which is restricted to 1GB
# mkdir $SCRATCHDIR/temp
# export TMPDIR=$SCRATCHDIR/temp
####################################################################################################
### Collapsing exact duplicates
# WARNING !!!
# You can NOT deduplex SpikeIn samples!!!!! The LIBRARIES WERE CREATED FROM PCR AMPLICONS !!!
# for z in *$APPENDIX
# do
# 	READ_LIB=$z
# 	echo "Now I am decompressing PE reads $READ_LIB"
# 	unpigz -v -p $THREADS $READ_LIB
# 	echo "Done decompressing PE reads $READ_LIB"
# done
# # wait
# wait
#
# for y in *$APPENDIX3
# do
# 	READ_LIB2=$y
# 	SAMPLENAME=${y%.*}
# 	echo "Now I am collapsing exact duplicates in PE reads $READ_LIB2" 
# 	$FASTQ2COLLAPSE -v $READ_LIB2 $SAMPLENAME".dedupl.fastq" &
# 	# fastq2collapse.pl script uses awk, sort and other bash tools that will exceed /var/tmp disk quota on Metacentrum.
# 	# Setting biger scratch will not help you!
# 	# You need to redirect temp files from /var/temp which is restricted to 1GB!
# 	# Use command "export TMPDIR=$SCRATCHDIR" 
# done
# #wait
# wait
# echo "Done collapsing all PE libraries"
# echo "These deduplicated files were created"
# ls -c1 *.dedupl.fastq
####################################################################################################
## compresing deduplicated files
# for x in *$APPENDIX4
# do
# 	READ_LIB3=$x
# 	echo "Now I am compressing PE reads $READ_LIB3"
# 	pigz -v -p $THREADS $READ_LIB3
# 	echo "Done compressing PE reads $READ_LIB3"
# done
# # wait
# wait
# rm $SCRATCHDIR/*_1.fastq $SCRATCHDIR/*_2.fastq # deleting original data
####################################################################################################
# Set extra flags for STAR
module unload ctk-1.0.7 # Needed to unload ctk modules due to incompatibility with STAR module
msg=""
extra_flags_star_motif="" # Default for SE as in the manual
extra_flags_star_wig=""
if [ "$STRANDED" == "true" ]; then
	extra_flags_star_motif="" # For STAR outSAMstrandField
	extra_flags_star_wig="--outWigStrand Stranded" # For START bedGraph
	msg="Stranded experiment"
else
	extra_flags_star_motif="--outSAMstrandField intronMotif" # For STAR outSAMstrandField
	extra_flags_star_wig="--outWigStrand Unstranded" # For START bedGraph
	msg="Unstranded experiment"
fi

echo "Running as $msg"
####################################################################################################
### Genome and annotation preparation
# GENOME_NAME=$(basename $MY_GENOME) # for normal alignment, not for RSEM analysis prep
# GTF_NAME=$(basename $RSEM_GTF) # for normal alignment, not for RSEM analysis prep
GENOME_NAME=$(basename $MY_GENOME)
GTF_NAME=$(basename $MY_GTF)
GEN_DIR=$(basename $GENOME_INDEX)
# unpigz -p $THREADS $GENOME_NAME # for normal alignment, not for RSEM analysis prep
# unpigz -p $THREADS $GTF_NAME # for normal alignment, not for RSEM analysis prep
# GENOME=${GENOME_NAME%.*} # for normal alignment, not for RSEM analysis prep
# GTF=${GTF_NAME%.*} # for normal alignment, not for RSEM analysis prep
GENOME=${GENOME_NAME%.*}
GTF=${GTF_NAME%.*}
#
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Decompressing files:"
unpigz -v -p $THREADS -d $GENOME_NAME
unpigz -v -p $THREADS -d $GTF_NAME
find "$SCRATCHDIR" -maxdepth 1 -name "*.fastq.gz" -exec unpigz -v -p $THREADS -d {} + 
#
$SAMTOOLS faidx $SCRATCHDIR/$GENOME & # For bedGraphToBigWig
####################################################################################################
shopt nullglob # If set, Bash allows filename patterns which match no files to expand to a null string, rather than themselves.
IFS=$'\n' # split on newline only, needed for filling arrays with output from find
set -f    # disable globbing, needed for filling arrays with output from find
# declare an array variables
declare -a OAT_cyto=($(find . -maxdepth 1 -name 'DIS3l2_upscale_OAT_cyto_flexbar*.fastq' -exec basename {} \; | sort -n))
declare -a OAT_nonfrac=($(find . -maxdepth 1 -name 'DIS3l2_upscale_OAT_nonfrac_flexbar*.fastq' -exec basename {} \; | sort -n))
declare -a OAT_nucl=($(find . -maxdepth 1 -name 'DIS3l2_upscale_OAT_nucl_flexbar*.fastq' -exec basename {} \; | sort -n))
declare -a U12_cyto=($(find . -maxdepth 1 -name 'DIS3l2_upscale_U12_cyto_flexbar*.fastq' -exec basename {} \; | sort -n))
declare -a U12_nonfrac=($(find . -maxdepth 1 -name 'DIS3l2_upscale_U12_nonfrac_flexbar*.fastq' -exec basename {} \; | sort -n))
declare -a U12_nucl=($(find . -maxdepth 1 -name 'DIS3l2_upscale_U12_nucl_flexbar*.fastq' -exec basename {} \; | sort -n))
declare -a ARRAY_NAMES=("OAT_cyto" "OAT_nonfrac" "OAT_nucl" "U12_cyto" "U12_nonfrac" "U12_nucl")
declare -a ALL_ARRAYS=("${OAT_cyto[@]}" "${OAT_nonfrac[@]}" "${OAT_nucl[@]}" "${U12_cyto[@]}" "${U12_nonfrac[@]}" "${U12_nucl[@]}" )
set +f # enable globbing, because is needed for other parts of the script, like variable expansions !!!
# get length of an array
OAT_cyto_num="${#OAT_cyto[@]}"
OAT_nonfrac_num="${#OAT_nonfrac[@]}"
OAT_nucl_num="${#OAT_nucl[@]}"
U12_cyto_num="${#U12_cyto[@]}"
U12_nonfrac_num="${#U12_nonfrac[@]}"
U12_nucl_num="${#U12_nucl[@]}"
quantity="${#ALL_ARRAYS[@]}"
lib_count="${#ARRAY_NAMES[@]}"
#
echo "There are $OAT_cyto_num OAT_cyto samples that will be mapped."
echo "Sample names are: ${OAT_cyto[@]}"
echo "There are $OAT_nonfrac_num OAT_nonfrac_num samples that will be mapped."
echo "Sample names are: ${OAT_nonfrac[@]}"
echo "There are $OAT_nucl_num OAT_nucl_num samples that will be mapped."
echo "Sample names are: ${OAT_nucl[@]}"
echo "There are $U12_cyto_num U12_cyto_num samples that will be mapped."
echo "Sample names are: ${U12_cyto[@]}"
echo "There are $U12_nonfrac_num U12_nonfrac_num samples that will be mapped."
echo "Sample names are: ${U12_nonfrac[@]}"
echo "There are $U12_nucl_num U12_nucl_num samples that will be mapped."
echo "Sample names are: ${U12_nucl[@]}"
# Set MAX RAM and CPUs
RAM=$[$MY_RAM-1] # lower it by one to leave some space
MEM_LIMIT=$RAM # $[$RAM/$THREADS] # In case we need GB per thread
# In case we need to divide CPUs between several running jobs
# To round up (Function ceiling):
THREADS_LIMIT="$(echo "$THREADS $lib_count" | awk '{print int( ($1/$2) + 1 )}')"
# To round down (Function floor):
# THREADS_LIMIT="$(echo "$THREADS $quantity" | awk '{print int($1/$2)}')"
#
echo "There are $quantity LIBRARIES that will be mapped."
echo "Sample names are: ${LIBRARIES[@]}"
### Alignment
# --limitSjdbInsertNsj 4000000 = need to increase limit for juction detection if STAR ends up after 1st pass with error
for (( w=0;w<${quantity};w +=2 ));
do
	SAMPLE="${ALL_ARRAYS[$w]}"
	SAMPLE2="${ALL_ARRAYS[$w+1]}"
	if [[ -f $SAMPLE ]]; then
		date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
		echo "Now I am processing PE reads $SAMPLE and $SAMPLE2 - alignment"
		$STAR --runMode alignReads --runThreadN $THREADS --genomeDir $GEN_DIR \
		--genomeLoad LoadAndRemove --readFilesIn $SAMPLE $SAMPLE2 \
		--sjdbOverhang $RD_LENGTH --sjdbGTFfile $GTF \
		--outFileNamePrefix ${SAMPLE%"_flexbar_"*}.${GTF%.*} \
		--outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
		--outFilterMismatchNmax $NO_MISMATCHES --outFilterMismatchNoverReadLmax $PER_MISMATCHES \
		--alignIntronMin 20 --alignIntronMax $MAX_INTRON --alignMatesGapMax 1000000 \
		--outFilterMatchNmin 0 --outFilterMatchNminOverLread 0.66 \
		--outSAMheaderHD @HD VN:1.4 SO:coordinate --chimSegmentMin 30 --chimOutType SeparateSAMold \
		--outReadsUnmapped Fastx --outSAMunmapped Within --outFilterType Normal --outSAMattributes All \
		$extra_flags_star_motif --quantMode GeneCounts TranscriptomeSAM --sjdbScore 1 --twopassMode Basic \
		--limitSjdbInsertNsj 10000000 --limitIObufferSize  1500000000 --outMultimapperOrder Random \
		--outSAMtype BAM SortedByCoordinate & # --readFilesCommand zcat for reading .gz files; --outWigType bedGraph # --outSAMstrandField intronMotif - ENCODE uses this just for UNSTRANDED mapping
		date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
		echo "Done processing PE reads $SAMPLE and $SAMPLE2 - alignment"
	else
		date +"%d/%m/%Y %H:%M:%S"
		echo "There is no $SAMPLE file for mapping!" && exit 1
	fi
done
wait
####################################################################################################
### STAR mark duplicates
#for i in *Aligned.sortedByCoord.out.bam 
#do
#$STAR --inputBAMfile $i --bamRemoveDuplicatesType UniqueIdentical --runMode inputAlignmentsFromBAM --bamRemoveDuplicatesMate2basesN 15 --outFileNamePrefix markdup. --limitBAMsortRAM 30000000000
#done

####################################################################################################
### Make signal bedGraph->bigWig files
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Creating BedGraph SIGNAL files"
#
mkdir $SCRATCHDIR/star_signal

for i in *Aligned.sortedByCoord.out.bam 
do
	echo "Sorting and renaming BAM file for UCSC browser:"
	echo "$i"
#	We need to change chromosome names to visualize the output in UCSC Browser http://seqanswers.com/forums/archive/index.php/t-22504.html
	$SAMTOOLS view -@ $THREADS -h $i | \
	awk 'BEGIN{FS=OFS="\t"} (/^@/ && !/@SQ/){print $0} $2~/^SN:[1-9]|^SN:X|^SN:Y|^SN:MT/{print $0} $3~/^[1-9]|X|Y|MT/{$3="chr"$3; print $0} ' | \
	sed 's/SN:/SN:chr/g' | sed 's/chrMT/chrM/g' | samtools view -@ $THREADS -b - > "$i.tmp"
	# Preparing STAR signal for IGV
	echo "Preparing STAR signal for IGV"
	$STAR --runMode inputAlignmentsFromBAM --runThreadN $THREADS --inputBAMfile $i --outWigType bedGraph $extra_flags_star_wig \
	--outWigNorm RPM --outFileNamePrefix $SCRATCHDIR/star_signal/"${i%.out.bam}"_IGV_
	# Preparing STAR signal for UCSC
	echo "Preparing STAR signal for UCSC"
	$STAR --runMode inputAlignmentsFromBAM --runThreadN $THREADS --inputBAMfile "$i.tmp" --outWigType bedGraph $extra_flags_star_wig \
	--outWigNorm RPM --outWigReferencesPrefix chr --outFileNamePrefix $SCRATCHDIR/star_signal/"${i%.out.bam}"_UCSC_ # --outWigReferencesPrefix chr suitable for UCSC
#	rm $i.tmp
done
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Finnished creating BedGraph SIGNAL files"
echo "Created files:"
ls -Rc1 $SCRATCHDIR/star_signal/*

# Prepare .fai for conversion
#	We need to change chromosome names to visualize the output in UCSC Browser http://seqanswers.com/forums/archive/index.php/t-22504.html
echo " Prepare .fai for conversion"
cat $SCRATCHDIR/${GENOME}.fai | sed -e 's/^\([0-9XY]\)/chr\1/' -e 's/^MT/chrM/' | grep "^chr" > $SCRATCHDIR/${GENOME}.fai.tmp
#
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Sorting BedGraph SIGNAL files"
for i in $SCRATCHDIR/star_signal/*.bg # https://ycl6.gitbooks.io/rna-seq-data-analysis/visualization.html
do
	FILE=$i
	if [[ "$FILE" =~ "_UCSC_" ]]; then
		FILENAME=${i%.*}
		date +"%d/%m/%Y %H:%M:%S"
		echo "Sorting BedGraph SIGNAL file for UCSC genome browser"
		echo "Filename is: $i"
		cat "$i" | sed -e 's/^\([0-9XY]\)/chr\1/' -e 's/^MT/chrM/' | \
		grep "^chr" > "$FILENAME".bg #	We need to change chromosome names to visualize the output in UCSC Browser http://seqanswers.com/forums/archive/index.php/t-22504.html
		LC_COLLATE=C sort -k1,1 -k2,2n -o "$FILENAME".tmp "$FILENAME".bg # Sort to proper order
		mv -v "$FILENAME".tmp "$FILENAME".bg
		echo "Converting file to BigWig format"
		$BEDGRAPH2BIGWIG "$FILENAME".bg $SCRATCHDIR/${GENOME}.fai.tmp ${i%.*}.wg
		rm "$FILENAME".bg
	elif [[ "$FILE" =~ "_IGV_" ]]; then
		FILENAME2=${i%.*}
		date +"%d/%m/%Y %H:%M:%S"
		echo "Sorting BedGraph SIGNAL file for IGV genome browser"
		echo "Filename is: $i"
		$BEDGRAPH2BIGWIG $FILE $SCRATCHDIR/${GENOME}.fai.tmp $FILENAME2.wg
	else
		date +"%d/%m/%Y %H:%M:%S"
		echo "There are no *.bg files for sorting BedGraph SIGNAL!"
	fi
	rm $SCRATCHDIR/${GENOME}.fai.tmp
done
#
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Finnished sorting BedGraph SIGNAL files"
#
mkdir -p $SCRATCHDIR/star_signal/bg; mkdir $SCRATCHDIR/star_signal/wg; mkdir $SCRATCHDIR/star_signal/log
mv -v $SCRATCHDIR/star_signal/*.bg $SCRATCHDIR/star_signal/bg/; mv $SCRATCHDIR/star_signal/*.wg $SCRATCHDIR/star_signal/wg/; mv $SCRATCHDIR/star_signal/*.out $SCRATCHDIR/star_signal/log/

for i in $SCRATCHDIR/star_signal/bg/*.bg # Compress .bg files
do 
	pigz $i
done &

####################################################################################################
### Cleaning results and indexing

# Move transcriptome mapping
echo "Move transcriptome mapping"
mkdir $SCRATCHDIR/transcriptome
mv -v $SCRATCHDIR/*Transcriptome.out.bam $SCRATCHDIR/transcriptome/

# Sort transcriptome BAMs
# Prepare for RSEM: sort transcriptome BAM to ensure the order of the reads, to make RSEM output (not pme) deterministic

date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Prepare for RSEM - Sort transcriptome BAM"
#
cd $SCRATCHDIR/transcriptome
#
for i in *Transcriptome.out.bam # Untested!
do
	cat <( $SAMTOOLS view -H $i ) <( $SAMTOOLS view -@ $THREADS $i | \
		awk '{printf $0 " "; getline; print}' | \
		sort -S ${MEM_LIMIT}G -T ./ | tr ' ' '\n' ) | \
		$SAMTOOLS view -@ $THREADS -b - > $i.tmp
	mv $i.tmp $i
done
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Finnished preparation for RSEM - Sort transcriptome BAM"
#
# Move STAR logs
mkdir $SCRATCHDIR/star_log
mv -v $SCRATCHDIR/*.out $SCRATCHDIR/star_log/

$MULTIQC -o $SCRATCHDIR/star_log $SCRATCHDIR/star_log/ &

# Chimeric to BAM
cd $SCRATCHDIR/

mkdir $SCRATCHDIR/chimeric
#
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Chimeric from SAM to BAM conversion"
#
for i in *.sam
do
	$SAMTOOLS view -@ $THREADS -b $i | $SAMTOOLS sort -@ $THREADS -T $SCRATCHDIR/tmp.sort \
	-o $SCRATCHDIR/chimeric/${i%.*}.bam -

	rm $i

	$SAMTOOLS index -@ $THREADS $SCRATCHDIR/chimeric/${i%.*}.bam &
done
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Finnished conversion of Chimeric SAM to BAM."

# Move STAR gene counts - should be mainly used for RSEM as stradness determination help/confirmation
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo " Move STAR gene counts - mainly for RSEM usage"
mkdir -p $SCRATCHDIR/star_gc
mv -v $SCRATCHDIR/*ReadsPerGene.out.tab $SCRATCHDIR/star_gc
#
$MULTIQC -o $SCRATCHDIR/star_gc $SCRATCHDIR/star_gc/
#
wait
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Finnished moving STAR gene counts - mainly for RSEM usage"
####################################################################################################
### Finalize and copy results
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Finalize and copy results"
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR_QC
# rm -v $SCRATCHDIR/*.dedupl.fastq # only when the deduplicated reads were created in previous run

mkdir -p $SCRATCHDIR/alignment/genome
mv -v $SCRATCHDIR/*.bam $SCRATCHDIR/alignment/genome/
mv -v $SCRATCHDIR/*.bai $SCRATCHDIR/alignment/genome/
mv -v $SCRATCHDIR/transcriptome $SCRATCHDIR/alignment/
mkdir -p $SCRATCHDIR/other/junction
mv -v $SCRATCHDIR/*junction $SCRATCHDIR/chimeric/
mv -v $SCRATCHDIR/*SJ.out.tab $SCRATCHDIR/other/junction/
mv -v $SCRATCHDIR/chimeric $SCRATCHDIR/other/

rm $SCRATCHDIR/*.fa.gz
rm $SCRATCHDIR/$GENOME*
rm $SCRATCHDIR/$GTF*
rm -rv $SCRATCHDIR/$GEN_DIR
rm -rv $SCRATCHDIR/*"_STARpass1"
rm -rv $SCRATCHDIR/*"_STARgenome"
rm -rv $SCRATCHDIR/*"_STARtmp"

cp -avr $SCRATCHDIR/star_signal $OUTPUT_DIR_QC/
cp -avr $SCRATCHDIR/star_log $OUTPUT_DIR_QC/
cp -avr $SCRATCHDIR/star_gc $OUTPUT_DIR_QC/
rm -rv $SCRATCHDIR/star_log $SCRATCHDIR/star_gc $SCRATCHDIR/star_signal
cp -avr $SCRATCHDIR/* $OUTPUT_DIR/ || export CLEAN_SCRATCH=false
#
rm -rv $SCRATCHDIR/*
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"