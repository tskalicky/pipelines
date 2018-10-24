#!/bin/bash
#PBS -l nodes=1:ppn=30
#PBS -l walltime=96:00:00
#PBS -l mem=200gb
#PBS -k oe
#PBS -N STAR_Dis3L2_human
#PBS -M tomas.skalicky@seznam.cz
#PBS -m abe
#
#PARU KRTECEK server is using TORQUE scheduling system !!!
#
## initialize the required application
# For 1 library you need cca 35GB RAM
#
# Alignment - PE RNA-Seq
# Requires STAR, pigz, samtools, multiQC, bedGraphToBigWig
############################################################################################
# Adding RSEM and STAR binaries into the PATH
# Need to export system PATH on Krtecek server or bash commands wil NOT work111 :-/
# Will work ONLY on NFS4 connected servers
# export PATH="/storage/brno3-cerit/home/tskalicky/tools/RSEM-1.3.0/bin:$PATH"
export PATH="/bin/:$PATH"
export PATH="/usr/bin/:$PATH"
export PATH="/home/users/tskalicky/anaconda2/bin/:$PATH"
echo "PATH after modification is:"
echo "$PATH"
#
# Set number of CPU
THREADS=$PBS_NUM_PPN
# Set MAX ram
MY_RAM="200" # SET THE AMOUNT OF RAM HERE !!!
RAM=$[$MY_RAM-1] # lower it by one to leave some space
MEM_LIMIT=$RAM # $[$RAM/$THREADS] # In case we need GB per thread
#
# Note: do not use modified GTF (added features) to alignment as it can cause issues later on
############################################################################################
### Variables
INPUT_DIR="/home/users/tskalicky/CEITEC/Dis3L2/Dasa_SpikeIn/mapping/Scerevisiae_genome/unmapped_reads"
OUTPUT_DIR="/home/users/tskalicky/CEITEC/Dis3L2/Dasa_SpikeIn/mapping/human_genome/STAR/temp_spikein"
OUTPUT_DIR=${OUTPUT_DIR}/results
OUTPUT_DIR_QC=${OUTPUT_DIR}/qc
# Genome and annotation
MY_GENOME="/home/users/tskalicky/CEITEC/genomes/human/ensembl91/GRCh38.dna.primary_assembly.fa.gz"
MY_GTF="/home/users/tskalicky/CEITEC/genomes/human/ensembl91/GRCh38.91.gtf.gz"
GENOME_INDEX="/home/users/tskalicky/CEITEC/genomes/human/ensembl91/STAR_index/GRCh38_with_tRNA_RD_LENGTH_100_STAR2.6c" 
# Genome index with --sjdbOverhang 79 - Only for Dis3L2 RNA-seq. length=125
#
STRANDED="false" # "true" or "false" - STAR parameters: strandedness, affects bedGraph (wiggle) files and XS tag in BAM 
#
APPENDIX1="_1.fastq.gz"
APPENDIX2="_2.fastq.gz"
APPENDIX=".fastq.gz"
APPENDIX3=".fastq"
APPENDIX4=".dedupl.fastq"
APPENDIX5=".dedupl.fastq.gz"
APPENDIX6="_1.dedupl.fastq.gz"
APPENDIX7="_2.dedupl.fastq.gz"
#
NO_MISMATCHES=999 # Number of mismatches per read - since we do not trim; 14 is recommended in QuantSeq FWD protocol for SE and 16 for PE; to turn this off set it to 999
PER_MISMATCHES=0.04 # Percent of mismatches per read; I usually keep this at 0.05 or 0.1; for untrimmed QuantSeq FWD it might be 0.3 -> for 50 bp single end it gives max. 15 mismatches
MAX_INTRON=1000000 # Default used by ENCODE is 1000000; to turn this off set 1
MAX_MATE_DIST=1000000 # Default used by ENCODE is 1000000; For "normal" fragments 1000 should be enough but for special cases, like chimeric we should increase this
#
# Read length for sjdbOverhang; --sjdbOverhang 100 should work fine for most of the data, but more specific setting based on the real read length should be more sensitive https://groups.google.com/forum/#!msg/rna-star/h9oh10UlvhI/BfSPGivUHmsJ
RD_LENGTH=100 # Read length from the sequencing. Illumina sometimes reports N+1 http://seqanswers.com/forums/archive/index.php/t-31154.html; in case you change this value uncomment next line as well
#RD_LENGTH=$[$RD_LENGTH-$TRIM_LEFT-$TRIM_RIGHT-1] # Modified read length for the index creation and mapping - should be read length -1
#
# Binaries
UNPIGZ=$(which unpigz)
STAR=$(which STAR)
# BBMAP=$(which bbmap.sh)
# SAMTOOLS=$(which samtools)
SAMTOOLS="/home/users/tskalicky/anaconda2/bin/samtools"
BEDGRAPH2BIGWIG=$(which bedGraphToBigWig)
MULTIQC=$(which multiqc)
#
# FASTQ2COLLAPSE=$(which fastq2collapse.pl)
#
# Check the tools versions
echo "Checking required tools:"
# which $BBMAP
which $SAMTOOLS
$SAMTOOLS --version
#which $BEDGRAPH2BIGWIG
which $MULTIQC
which $UNPIGZ
which $STAR
$STAR --version
#
####################################################################################################
# copy input data
# use cp -avr when copying directories
mkdir -p $OUTPUT_DIR
#
cp -av $MY_GENOME $MY_GTF $OUTPUT_DIR
find "$INPUT_DIR" -maxdepth 1 -name *.fq.gz -exec cp -vt "$OUTPUT_DIR" {} +
cp -avr $GENOME_INDEX $OUTPUT_DIR/ 
cd $OUTPUT_DIR
#
if [ ! -d "$OUTPUT_DIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if output directory is created
echo "OUTPUT_DIR path is:" $OUTPUT_DIR
echo "Following files were copied to scratch:"
ls -Rc1
#
####################################################################################################
###
# Setting Variables
shopt nullglob # If set, Bash allows filename patterns which match no files to expand to a null string, rather than themselves.
IFS=$'\n' # split on newline only, needed for filling arrays with output from find
set -f    # disable globbing, ONLY needed for filling arrays with output from find
# declare an array variables
declare -a OAT_CYTO=($(find . -maxdepth 1 -name  '*OAT_cyto*.fq.gz' -exec basename {} \; | sort -n))
declare -a OAT_NONFRAC=($(find . -maxdepth 1 -name  '*OAT_nonfrac*.fq.gz' -exec basename {} \; | sort -n))
declare -a OAT_NUCL=($(find . -maxdepth 1 -name  '*OAT_nucl*.fq.gz' -exec basename {} \; | sort -n))
declare -a U12_CYTO=($(find . -maxdepth 1 -name  '*U12_cyto*.fq.gz' -exec basename {} \; | sort -n))
declare -a U12_NONFRAC=($(find . -maxdepth 1 -name  '*U12_nonfrac*.fq.gz' -exec basename {} \; | sort -n))
declare -a U12_NUCL=($(find . -maxdepth 1 -name  '*U12_nucl*.fq.gz' -exec basename {} \; | sort -n))
declare -a ARRAY_NAMES=("OAT_CYTO" "OAT_NONFRAC" "OAT_NUCL" "U12_CYTO" "U12_NONFRAC" "U12_NUCL")
declare -a ALL_ARRAYS=("${OAT_CYTO[@]}" "${OAT_NONFRAC[@]}" "${OAT_NUCL[@]}" "${U12_CYTO[@]}" "${U12_NONFRAC[@]}" "${U12_NUCL[@]}" )
# declare -a ARRAY_NAMES=("OAT_NONFRAC" "OAT_NUCL" "U12_CYTO" "U12_NONFRAC" "U12_NUCL")
# declare -a ALL_ARRAYS=("${OAT_NONFRAC[@]}" "${OAT_NUCL[@]}" "${U12_CYTO[@]}" "${U12_NONFRAC[@]}" "${U12_NUCL[@]}" )
set +f # enable globbing, because is needed for other parts of the script, like variable expansions !!!
# get length of an array
# OAT_cyto_num="${#OAT_CYTO[@]}"
OAT_nonfrac_num="${#OAT_NONFRAC[@]}"
OAT_nucl_num="${#OAT_NUCL[@]}"
U12_cyto_num="${#U12_CYTO[@]}"
U12_nonfrac_num="${#U12_NONFRAC[@]}"
U12_nucl_num="${#U12_NUCL[@]}"
quantity="${#ALL_ARRAYS[@]}"
#
####################################################################################################
### Genome and annotation preparation
GENOME_NAME=$(basename $MY_GENOME)
GTF_NAME=$(basename $MY_GTF)
GEN_DIR=$(basename $GENOME_INDEX)

GENOME=${GENOME_NAME%.*}
GTF=${GTF_NAME%.*}
#
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Decompressing files:"
# $UNPIGZ -v -p $THREADS -d $GENOME_NAME
# unpigz -v -p $THREADS -d $GTF_NAME
gzip -vd $GENOME_NAME
gzip -vd $GTF_NAME
#
### Commands 
## STAR - Alignment
$SAMTOOLS faidx $OUTPUT_DIR/$GENOME & # For bedGraphToBigWig
#
echo "There are $OAT_cyto_num OAT_CYTO samples that will be mapped."
echo "Sample names are: ${OAT_CYTO[@]}"
echo "There are $OAT_nonfrac_num OAT_nonfrac_num samples that will be mapped."
echo "Sample names are: ${OAT_NONFRAC[@]}"
echo "There are $OAT_nucl_num OAT_nucl_num samples that will be mapped."
echo "Sample names are: ${OAT_NUCL[@]}"
echo "There are $U12_cyto_num U12_cyto_num samples that will be mapped."
echo "Sample names are: ${U12_CYTO[@]}"
echo "There are $U12_nonfrac_num U12_nonfrac_num samples that will be mapped."
echo "Sample names are: ${U12_NONFRAC[@]}"
echo "There are $U12_nucl_num U12_nucl_num samples that will be mapped."
echo "Sample names are: ${U12_NUCL[@]}"
#
for (( w=0;w<${quantity};w +=2 )); do
	SAMPLE="${ALL_ARRAYS[$w]}"
	SAMPLE2="${ALL_ARRAYS[$w+1]}"
	SAMPLENAME=${SAMPLE%"_unmapped_only_"*}
	if [[ -f $SAMPLE ]]; then
		date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
		echo "Now I am mapping PE reads $SAMPLE and $SAMPLE2 - alignment"
		$STAR --runMode alignReads --runThreadN $THREADS --genomeDir $GEN_DIR \
		--readFilesIn $SAMPLE $SAMPLE2 \
		--readFilesCommand zcat --sjdbOverhang $RD_LENGTH --sjdbGTFfile $GTF \
		--outFileNamePrefix "$SAMPLENAME".${GTF%.*} \
		--outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
		--outFilterMismatchNmax $NO_MISMATCHES --outFilterMismatchNoverReadLmax $PER_MISMATCHES \
		--alignIntronMin 20 --alignIntronMax $MAX_INTRON --alignMatesGapMax 1000000 \
		--outFilterMatchNmin 0 --outFilterMatchNminOverLread 0.66 \
		--outSAMheaderHD @HD VN:1.4 SO:coordinate --chimSegmentMin 30 --chimOutType SeparateSAMold \
		--outReadsUnmapped Fastx --outSAMunmapped None --outFilterType Normal --outSAMattributes All \
		$extra_flags_star_motif --quantMode GeneCounts TranscriptomeSAM --sjdbScore 1 --twopassMode Basic \
		--limitSjdbInsertNsj 10000000 --limitIObufferSize  1500000000 --outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate # \ --outWigType bedGraph # --outSAMstrandField intronMotif - ENCODE uses this just for UNSTRANDED mapping
	else
		date +"%d/%m/%Y %H:%M:%S"
		echo "There is no $SAMPLE file for mapping!" && exit 1
	fi
done
wait
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Done mapping all PE reads in parralel."
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
mkdir $OUTPUT_DIR/star_signal

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
	--outWigNorm RPM --outFileNamePrefix $OUTPUT_DIR/star_signal/"${i%.out.bam}"_IGV_
	# Preparing STAR signal for UCSC
	echo "Preparing STAR signal for UCSC"
	$STAR --runMode inputAlignmentsFromBAM --runThreadN $THREADS --inputBAMfile "$i.tmp" --outWigType bedGraph $extra_flags_star_wig \
	--outWigNorm RPM --outWigReferencesPrefix chr --outFileNamePrefix $OUTPUT_DIR/star_signal/"${i%.out.bam}"_UCSC_ # --outWigReferencesPrefix chr suitable for UCSC
#	rm $i.tmp
done
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Finnished creating BedGraph SIGNAL files"
echo "Created files:"
ls -Rc1 $OUTPUT_DIR/star_signal/*

# Prepare .fai for conversion
#	We need to change chromosome names to visualize the output in UCSC Browser http://seqanswers.com/forums/archive/index.php/t-22504.html
echo " Prepare .fai for conversion"
cat $OUTPUT_DIR/${GENOME}.fai | sed -e 's/^\([0-9XY]\)/chr\1/' -e 's/^MT/chrM/' | grep "^chr" > $OUTPUT_DIR/${GENOME}.fai.tmp
#
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Sorting BedGraph SIGNAL files"
for i in $OUTPUT_DIR/star_signal/*.bg # https://ycl6.gitbooks.io/rna-seq-data-analysis/visualization.html
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
		$BEDGRAPH2BIGWIG "$FILENAME".bg $OUTPUT_DIR/${GENOME}.fai.tmp ${i%.*}.wg
		rm "$FILENAME".bg
	elif [[ "$FILE" =~ "_IGV_" ]]; then
		FILENAME2=${i%.*}
		date +"%d/%m/%Y %H:%M:%S"
		echo "Sorting BedGraph SIGNAL file for IGV genome browser"
		echo "Filename is: $i"
		$BEDGRAPH2BIGWIG $FILE $OUTPUT_DIR/${GENOME}.fai.tmp $FILENAME2.wg
	else
		date +"%d/%m/%Y %H:%M:%S"
		echo "There are no *.bg files for sorting BedGraph SIGNAL!"
	fi
	rm $OUTPUT_DIR/${GENOME}.fai.tmp
done
#
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Finnished sorting BedGraph SIGNAL files"
#
mkdir -p $OUTPUT_DIR/star_signal/bg; mkdir $OUTPUT_DIR/star_signal/wg; mkdir $OUTPUT_DIR/star_signal/log
mv -v $OUTPUT_DIR/star_signal/*.bg $OUTPUT_DIR/star_signal/bg/; mv $OUTPUT_DIR/star_signal/*.wg $OUTPUT_DIR/star_signal/wg/; mv $OUTPUT_DIR/star_signal/*.out $OUTPUT_DIR/star_signal/log/

for i in $OUTPUT_DIR/star_signal/bg/*.bg # Compress .bg files
do 
	gzip -v $i
done &

####################################################################################################
### Cleaning results and indexing

# Move transcriptome mapping
echo "Move transcriptome mapping"
mkdir $OUTPUT_DIR/transcriptome
mv -v $OUTPUT_DIR/*Transcriptome.out.bam $OUTPUT_DIR/transcriptome/

# Sort transcriptome BAMs
# Prepare for RSEM: sort transcriptome BAM to ensure the order of the reads, to make RSEM output (not pme) deterministic

date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Prepare for RSEM - Sort transcriptome BAM"
#
cd $OUTPUT_DIR/transcriptome
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
mkdir $OUTPUT_DIR/star_log
mv -v $OUTPUT_DIR/*.out $OUTPUT_DIR/star_log/

$MULTIQC -o $OUTPUT_DIR/star_log $OUTPUT_DIR/star_log/ &

# Chimeric to BAM
cd $OUTPUT_DIR/

mkdir $OUTPUT_DIR/chimeric
#
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Chimeric from SAM to BAM conversion"
#
for i in *.sam
do
	$SAMTOOLS view -@ $THREADS -b $i | $SAMTOOLS sort -@ $THREADS -T $OUTPUT_DIR/tmp.sort \
	-o $OUTPUT_DIR/chimeric/${i%.*}.bam -

	rm $i

	$SAMTOOLS index -@ $THREADS $OUTPUT_DIR/chimeric/${i%.*}.bam &
done
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Finnished conversion of Chimeric SAM to BAM."

# Move STAR gene counts - should be mainly used for RSEM as stradness determination help/confirmation
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo " Move STAR gene counts - mainly for RSEM usage"
mkdir -p $OUTPUT_DIR/star_gc
mv -v $OUTPUT_DIR/*ReadsPerGene.out.tab $OUTPUT_DIR/star_gc
#
$MULTIQC -o $OUTPUT_DIR/star_gc $OUTPUT_DIR/star_gc/
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
# rm $OUTPUT_DIR/*.dedupl.fastq.gz # only when the deduplicated reads were created in previous run

mkdir -p $OUTPUT_DIR/alignment/genome
mv -v $OUTPUT_DIR/*.bam $OUTPUT_DIR/alignment/genome/
mv -v $OUTPUT_DIR/*.bai $OUTPUT_DIR/alignment/genome/
mv -v $OUTPUT_DIR/transcriptome $OUTPUT_DIR/alignment/
mkdir -p $OUTPUT_DIR/other/junction
mv -v $OUTPUT_DIR/*junction $OUTPUT_DIR/chimeric/
mv -v $OUTPUT_DIR/*SJ.out.tab $OUTPUT_DIR/other/junction/
mv -v $OUTPUT_DIR/chimeric $OUTPUT_DIR/other/

rm $OUTPUT_DIR/*.fa.gz
rm $OUTPUT_DIR/$GENOME*
rm $OUTPUT_DIR/$GTF*
rm -r $OUTPUT_DIR/$GEN_DIR
rm -r $OUTPUT_DIR/*_STARpass1
rm -r $OUTPUT_DIR/*_STARgenome
rm -r $OUTPUT_DIR/*_STARtmp

cp -avr $OUTPUT_DIR/star_signal $OUTPUT_DIR_QC/
cp -avr $OUTPUT_DIR/star_log $OUTPUT_DIR_QC/
cp -avr $OUTPUT_DIR/star_gc $OUTPUT_DIR_QC/
rm -r $OUTPUT_DIR/star_log $OUTPUT_DIR/star_gc $OUTPUT_DIR/star_signal
# cp -avr $OUTPUT_DIR/* $OUTPUT_DIR/
#
# rm -r $OUTPUT_DIR/*
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"