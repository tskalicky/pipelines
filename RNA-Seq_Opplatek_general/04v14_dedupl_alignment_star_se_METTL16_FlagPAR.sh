#!/bin/bash
#PBS -l walltime=96:0:0 
#PBS -q default@wagap-pro.cerit-sc.cz 
#PBS -l select=1:ncpus=10:mem=50gb:scratch_local=250gb:os=debian9
#PBS -j oe
#PBS -N 04_dedupl_align_se_METL16_FlagPAT
#
## initialize the required application
module add star-2.5.2b
module add samtools-1.6
module add python27-modules-intel #required by multiQC
module add bedtools-2.26.0
module add ctk-1.0.7
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

INPUT_DIR=/storage/brno3-cerit/home/tskalicky/METTL16/preprocessed/data/final/FlagPAR
OUTPUT_DIR=/storage/brno3-cerit/home/tskalicky/FTO/mapping/deduplicated/FlagPAR
OUTPUT_DIR=${OUTPUT_DIR}/results
OUTPUT_DIR_QC=${OUTPUT_DIR}/qc

MY_RAM=50 # Max RAM memory for Samtools sort
THREADS=$PBS_NUM_PPN

STRANDED="false" # "true" or "false" - STAR parameters: strandedness, affects bedGraph (wiggle) files and XS tag in BAM 
# FTO libraries are not strand specific

#APPENDIX1="_R1_trim.fastq.gz"
#APPENDIX2="_R2_trim.fastq.gz"
APPENDIX=".fastq.gz"
APPENDIX2=".dedupl.fastq"

NO_MISMATCHES=999 # Number of mismatches per read - since we do not trim; 14 is recommended in QuantSeq FWD protocol for SE and 16 for PE; to turn this off set it to 999
PER_MISMATCHES=0.05 # Percent of mismatches per read; I usually keep this at 0.05 or 0.1; for untrimmed QuantSeq FWD it might be 0.3 -> for 50 bp single end it gives max. 15 mismatches
MAX_INTRON=1000000 # Default used by ENCODE is 1000000; to turn this off set 1
MAX_MATE_DIST=1000000 # Default used by ENCODE is 1000000; For "normal" fragments 1000 should be enough but for special cases, like chimeric we should increase this

# Read length for sjdbOverhang; --sjdbOverhang 100 should work fine for most of the data, but more specific setting based on the real read length should be more sensitive https://groups.google.com/forum/#!msg/rna-star/h9oh10UlvhI/BfSPGivUHmsJ
RD_LENGTH=49 # Read length from the sequencing. Illumina sometimes reports N+1 http://seqanswers.com/forums/archive/index.php/t-31154.html; in case you change this value uncomment next line as well
#RD_LENGTH=$[$RD_LENGTH-$TRIM_LEFT-$TRIM_RIGHT-1] # Modified read length for the index creation and mapping - should be read length -1

# Genome and annotation
GENOME_INDEX="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/STAR_index/Homo_sapiens.GRCh38.dna.primary_assembly_RD_LENGTH_49" # Genome index with --sjdbOverhang 124 - Only for FTO CLIP and INPUT with lib. length=50
MY_GENOME="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
MY_GTF="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/Homo_sapiens.GRCh38.91.gtf.gz"

# Binaries
STAR=$(which STAR)
SAMTOOLS=$(which samtools)
BEDGRAPH2BIGWIG=/storage/brno3-cerit/home/tskalicky/tools/ucsc/bedGraphToBigWig
MULTIQC=$(which multiqc)
FASTQ2COLLAPSE=$(which fastq2collapse.pl)

# Check the tools versions
which $STAR
which $SAMTOOLS
which $BEDGRAPH2BIGWIG
which $MULTIQC
which $FASTQ2COLLAPSE

####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
cd $INPUT_DIR
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
cp -av $MY_GENOME $MY_GTF $INPUT_DIR/*$APPENDIX $SCRATCHDIR
cp -avr $GENOME_INDEX $SCRATCHDIR/ 
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR
echo "Following reads were copied to scratch:"
ls -c1
####################################################################################################
### Collapsing exact duplicates
for z in *$APPENDIX
do
	READ_LIB=$z
	echo "Now I am decompressing SE reads $READ_LIB"
	unpigz -p $THREADS $READ_LIB
	echo "Done decompressing SE reads $READ_LIB"
done
# wait
wait
#
for y in *.fastq
do
	READ_LIB2=$y
	SAMPLENAME=${y%.*}
	echo "Now I am collapsing exact duplicates in SE reads $READ_LIB2" 
	$FASTQ2COLLAPSE -v $READ_LIB2 $SAMPLENAME".dedupl.fastq" &
done
#wait
wait
echo "Done collapsing all SE libraries"
#####################################################################################################
# Set extra flags for STAR
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

GENOME_NAME=$(basename $MY_GENOME)
GTF_NAME=$(basename $MY_GTF)
GEN_DIR=$(basename $GENOME_INDEX)

unpigz -p $THREADS $GENOME_NAME
unpigz -p $THREADS $GTF_NAME

GENOME=${GENOME_NAME%.*}
GTF=${GTF_NAME%.*}

# Set MAX ram
RAM=$[$MY_RAM-1] # lower it by one to leave some space
MEM_LIMIT=$RAM # $[$RAM/$THREADS] # In case we need GB per thread

$SAMTOOLS faidx $SCRATCHDIR/$GENOME & # For bedGraphToBigWig

####################################################################################################
### Alignment
for i in *$APPENDIX2
do
	READ_FOR=$i
	echo "Now I am processing SE reads $READ_FOR - alignment"

	$STAR --runMode alignReads --runThreadN $THREADS --genomeDir $GEN_DIR \
	--readFilesIn $READ_FOR \
	--readFilesCommand zcat --sjdbOverhang $RD_LENGTH --sjdbGTFfile $GTF \
	--outFileNamePrefix ${i%.fastq*}.${GTF%.*} \
	--outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
	--outFilterMismatchNmax $NO_MISMATCHES --outFilterMismatchNoverReadLmax $PER_MISMATCHES \
	--alignIntronMin 20 --alignIntronMax $MAX_INTRON --alignMatesGapMax 1000000 \
	--outFilterMatchNmin 0 --outFilterMatchNminOverLread 0.66 \
	--outSAMheaderHD @HD VN:1.6 SO:coordinate --chimSegmentMin 30 --chimOutType SeparateSAMold \
	--outSAMunmapped Within --outFilterType BySJout --outSAMattributes All \
	$extra_flags_star_motif --quantMode GeneCounts TranscriptomeSAM --sjdbScore 1 --twopassMode Basic \
	--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate # \ --outWigType bedGraph # --outSAMstrandField intronMotif - ENCODE uses this just for UNSTRANDED mapping

	echo "Done processing SE reads $READ_FOR - alignment"
done

####################################################################################################
### STAR mark duplicates
#for i in *Aligned.sortedByCoord.out.bam 
#do
#$STAR --inputBAMfile $i --bamRemoveDuplicatesType UniqueIdentical --runMode inputAlignmentsFromBAM --bamRemoveDuplicatesMate2basesN 15 --outFileNamePrefix markdup. --limitBAMsortRAM 30000000000
#done

####################################################################################################
### Make signal bedGraph->bigWig files
mkdir $SCRATCHDIR/star_signal

for i in *Aligned.sortedByCoord.out.bam 
do
	echo $i
#	We need to change chromosome names to visualize the output in UCSC Browser http://seqanswers.com/forums/archive/index.php/t-22504.html
	#$SAMTOOLS view -@ $THREADS -h $i | awk 'BEGIN{FS=OFS="\t"} (/^@/ && !/@SQ/){print $0} $2~/^SN:[1-9]|^SN:X|^SN:Y|^SN:MT/{print $0} $3~/^[1-9]|X|Y|MT/{$3="chr"$3; print $0} ' | sed 's/SN:/SN:chr/g' | sed 's/chrMT/chrM/g' | samtools view -@ $THREADS -b - > $i.tmp
	$STAR --runMode inputAlignmentsFromBAM --inputBAMfile $i --outWigType bedGraph $extra_flags_star_wig --outFileNamePrefix $SCRATCHDIR/star_signal/${i%.bam} # --outWigReferencesPrefix chr suitable for UCSC
#	rm $i.tmp
done

# Prepare .fai for conversion
#	We need to change chromosome names to visualize the output in UCSC Browser http://seqanswers.com/forums/archive/index.php/t-22504.html
cat $SCRATCHDIR/${GENOME}.fai | sed -e 's/^\([0-9XY]\)/chr\1/' -e 's/^MT/chrM/' | grep "^chr" > $SCRATCHDIR/${GENOME}.fai.tmp

for i in $SCRATCHDIR/star_signal/*.bg # https://ycl6.gitbooks.io/rna-seq-data-analysis/visualization.html
do
	echo $i
	cat $i | sed -e 's/^\([0-9XY]\)/chr\1/' -e 's/^MT/chrM/' | grep "^chr" > ${i%.*}_chr.bg #	We need to change chromosome names to visualize the output in UCSC Browser http://seqanswers.com/forums/archive/index.php/t-22504.html

	LC_COLLATE=C sort -k1,1 -k2,2n -o ${i}.tmp ${i%.*}_chr.bg # Sort to proper order
	mv ${i}.tmp ${i%.*}_chr.bg
	$BEDGRAPH2BIGWIG ${i%.*}_chr.bg $SCRATCHDIR/${GENOME}.fai.tmp ${i%.*}_chr.wg
#	rm ${i%.*}_chr.bg
done

rm $SCRATCHDIR/${GENOME}.fai.tmp

mkdir -p $SCRATCHDIR/star_signal/bg; mkdir $SCRATCHDIR/star_signal/wg; mkdir $SCRATCHDIR/star_signal/log
mv $SCRATCHDIR/star_signal/*.bg $SCRATCHDIR/star_signal/bg/; mv $SCRATCHDIR/star_signal/*.wg $SCRATCHDIR/star_signal/wg/; mv $SCRATCHDIR/star_signal/*.out $SCRATCHDIR/star_signal/log/

for i in $SCRATCHDIR/star_signal/bg/*.bg # Compress .bg files
do 
	pigz $i
done &
####################################################################################################
### Cleaning results and indexing

# Move transcriptome mapping
mkdir $SCRATCHDIR/transcriptome
mv $SCRATCHDIR/*Transcriptome.out.bam $SCRATCHDIR/transcriptome/

# Sort transcriptome BAMs
# Prepare for RSEM: sort transcriptome BAM to ensure the order of the reads, to make RSEM output (not pme) deterministic
cd $SCRATCHDIR/transcriptome

for i in *Transcriptome.out.bam
do
	cat <( $SAMTOOLS view -H $i ) <( $SAMTOOLS view -@ $THREADS $i | sort -S ${MEM_LIMIT}G -T ./ ) | \
		$SAMTOOLS view -@ $THREADS -b - > $i.tmp
	mv $i.tmp $i
done

# Move STAR logs
mkdir $SCRATCHDIR/star_log
mv $SCRATCHDIR/*.out $SCRATCHDIR/star_log/

$MULTIQC -o $SCRATCHDIR/star_log $SCRATCHDIR/star_log/ &

# Chimeric to BAM
cd $SCRATCHDIR/

mkdir $SCRATCHDIR/chimeric

for i in *.sam
do
	$SAMTOOLS view -@ $THREADS -b $i | $SAMTOOLS sort -@ $THREADS -T $SCRATCHDIR/tmp.sort \
	-o $SCRATCHDIR/chimeric/${i%.*}.bam -
	
	rm $i

	$SAMTOOLS index -@ $THREADS $SCRATCHDIR/chimeric/${i%.*}.bam &
done

# Move STAR gene counts - should be mainly used for RSEM as stradness determination help/confirmation
mkdir -p $SCRATCHDIR/star_gc
mv $SCRATCHDIR/*ReadsPerGene.out.tab $SCRATCHDIR/star_gc

$MULTIQC -o $SCRATCHDIR/star_gc $SCRATCHDIR/star_gc/

####################################################################################################
### Finalize and copy results

mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR_QC

mkdir -p $SCRATCHDIR/alignment/genome
mv $SCRATCHDIR/*.bam $SCRATCHDIR/alignment/genome/
mv $SCRATCHDIR/*.bai $SCRATCHDIR/alignment/genome/
mv $SCRATCHDIR/transcriptome $SCRATCHDIR/alignment/
mkdir -p $SCRATCHDIR/other/junction
mv $SCRATCHDIR/*junction $SCRATCHDIR/chimeric/
mv $SCRATCHDIR/*SJ.out.tab $SCRATCHDIR/other/junction/
mv $SCRATCHDIR/chimeric $SCRATCHDIR/other/

rm $SCRATCHDIR/*flexbar.fastq
rm $SCRATCHDIR/*.fa.gz
rm $SCRATCHDIR/$GENOME*
rm $SCRATCHDIR/$GTF*
rm -r $SCRATCHDIR/$GEN_DIR
rm -r $SCRATCHDIR/*_STARpass1
rm -r $SCRATCHDIR/*_STARgenome
rm -r $SCRATCHDIR/*_STARtmp

cp -r $SCRATCHDIR/star_signal $OUTPUT_DIR_QC/
cp -r $SCRATCHDIR/star_log $OUTPUT_DIR_QC/
cp -r $SCRATCHDIR/star_gc $OUTPUT_DIR_QC/
rm -r $SCRATCHDIR/star_log $SCRATCHDIR/star_gc $SCRATCHDIR/star_signal
cp -r $SCRATCHDIR/* $OUTPUT_DIR/

rm -r $SCRATCHDIR/*
