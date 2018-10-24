#!/bin/bash
#PBS -l select=1:ncpus=6:mem=36gb:scratch_local=150gb
#PBS -l walltime=24:00:00
#PBS -q default
#PBS -N 04_alignment_se_run
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
INPUT_DIR=/storage/brno3-cerit/home/opplatek/biocore/tichy_rna-seq_workshop/data/preprocessed
OUTPUT_DIR=/storage/brno3-cerit/home/opplatek/biocore/tichy_rna-seq_workshop
OUTPUT_DIR=${OUTPUT_DIR}/results
OUTPUT_DIR_QC=${OUTPUT_DIR}/qc

STRANDED="true" # "true" or "false" - STAR parameters: strandedness, affects bedGraph (wiggle) files and XS tag in BAM 

#APPENDIX1="_R1_trim.fastq.gz"
#APPENDIX2="_R2_trim.fastq.gz"
APPENDIX="_trim.fastq.gz"

NO_MISMATCHES=999 # Number of mismatches per read - since we do not trim; 14 is recommended in QuantSeq FWD protocol for SE and 16 for PE; to turn this off set it to 999
PER_MISMATCHES=0.05 # Percent of mismatches per read; I usually keep this at 0.05 or 0.1; for untrimmed QuantSeq FWD it might be 0.3 -> for 50 bp single end it gives max. 15 mismatches
MAX_INTRON=1000000 # Default used by ENCODE is 1000000; to turn this off set 1
MAX_MATE_DIST=1000000 # Default used by ENCODE is 1000000; For "normal" fragments 1000 should be enough but for special cases, like chimeric we should increase this

# Read length for sjdbOverhang; --sjdbOverhang 100 should work fine for most of the data, but more specific setting based on the real read length should be more sensitive https://groups.google.com/forum/#!msg/rna-star/h9oh10UlvhI/BfSPGivUHmsJ
RD_LENGTH=100 # Read length from the sequencing. Illumina sometimes reports N+1 http://seqanswers.com/forums/archive/index.php/t-31154.html; in case you change this value uncomment next line as well
#RD_LENGTH=$[$RD_LENGTH-$TRIM_LEFT-$TRIM_RIGHT-1] # Modified read length for the index creation and mapping - should be read length -1

# Genome and annotation
GENOME=/storage/brno2/home/opplatek/genomes/human/ensembl87/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
GENOME_INDEX=/storage/brno2/home/opplatek/genomes/human/ensembl87/STAR_index/Homo_sapiens.GRCh38.dna.primary_assembly # Default genome index with --sjdbOverhang 100 - STAR default
GTF=/storage/brno2/home/opplatek/genomes/human/ensembl87/Homo_sapiens.GRCh38.87.gtf.gz

RAM=34 # Max RAM memory for Samtools sort
THREADS=$PBS_NUM_PPN

# Binaries
STAR=/storage/brno2/home/opplatek/tools/STAR-2.5.2b/bin/Linux_x86_64/STAR
module add samtools-1.4
SAMTOOLS=$(which samtools)
BEDGRAPH2BIGWIG=/storage/brno2/home/opplatek/tools/ucsc/bedGraphToBigWig
module add python27-modules-gcc
#PYTHONPATH=$PYTHONPATH:/storage/brno2/home/opplatek/tools/MultiQC-1.0/lib/python2.7/site-packages/
#MULTIQC=/storage/brno2/home/opplatek/tools/MultiQC-1.0/bin/multiqc
MULTIQC=$(which multiqc)

# Check the tools versions
which $STAR
which $SAMTOOLS
which $BEDGRAPH2BIGWIG
which $MULTIQC

####################################################################################################
### Copy inputs
cp $GTF $SCRATCH/
cp $GENOME $SCRATCH/

cp -r $GENOME_INDEX $SCRATCH/

cp $INPUT_DIR/*$APPENDIX $SCRATCH/

cd $SCRATCH/

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
GENOME=$(basename $GENOME)
GTF=$(basename $GTF)

unpigz -p $THREADS $GTF
unpigz -p $THREADS $GENOME

GENOME=${GENOME%.gz*}
GTF=${GTF%.gz*}
GEN_DIR=$(basename $GENOME_INDEX)

# Set MAX ram
RAM=$[$RAM-1] # lower it by one to leave some space
MEM_LIMIT=$RAM # $[$RAM/$THREADS] # In case we need GB per thread

$SAMTOOLS faidx $SCRATCH/$GENOME & # For bedGraphToBigWig

####################################################################################################
### Alignment
for i in *$APPENDIX
do
	READ_FOR=$i
	echo "Now I am processing SE reads $READ_FOR - alignment"

	$STAR --runMode alignReads --runThreadN $THREADS --genomeDir $GEN_DIR \
	--readFilesIn $READ_FOR \
	--readFilesCommand zcat --sjdbOverhang $RD_LENGTH --sjdbGTFfile $GTF \
	--outFileNamePrefix ${i%.fastq.gz*}.${GTF%.*} \
	--outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
	--outFilterMismatchNmax $NO_MISMATCHES --outFilterMismatchNoverReadLmax $PER_MISMATCHES \
	--alignIntronMin 20 --alignIntronMax $MAX_INTRON --alignMatesGapMax 1000000 \
	--outFilterMatchNmin 0 --outFilterMatchNminOverLread 0.66 \
	--outSAMheaderHD @HD VN:1.4 SO:coordinate --chimSegmentMin 30 --chimOutType SeparateSAMold \
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
mkdir $SCRATCH/star_signal

for i in *Aligned.sortedByCoord.out.bam 
do
	echo $i
#	We need to change chromosome names to visualize the output in UCSC Browser http://seqanswers.com/forums/archive/index.php/t-22504.html
#	$SAMTOOLS view -@ $THREADS -h $i | awk 'BEGIN{FS=OFS="\t"} (/^@/ && !/@SQ/){print $0} $2~/^SN:[1-9]|^SN:X|^SN:Y|^SN:MT/{print $0} $3~/^[1-9]|X|Y|MT/{$3="chr"$3; print $0} ' | sed 's/SN:/SN:chr/g' | sed 's/chrMT/chrM/g' | samtools view -@ $THREADS -b - > $i.tmp
	$STAR --runMode inputAlignmentsFromBAM --inputBAMfile $i --outWigType bedGraph $extra_flags_star_wig --outFileNamePrefix $SCRATCH/star_signal/${i%.bam} # --outWigReferencesPrefix chr suitable for UCSC
#	rm $i.tmp
done

# Prepare .fai for conversion
#	We need to change chromosome names to visualize the output in UCSC Browser http://seqanswers.com/forums/archive/index.php/t-22504.html
cat $SCRATCH/${GENOME}.fai | sed -e 's/^\([0-9XY]\)/chr\1/' -e 's/^MT/chrM/' | grep "^chr" > $SCRATCH/${GENOME}.fai.tmp

for i in $SCRATCH/star_signal/*.bg # https://ycl6.gitbooks.io/rna-seq-data-analysis/visualization.html
do
	echo $i
	cat $i | sed -e 's/^\([0-9XY]\)/chr\1/' -e 's/^MT/chrM/' | grep "^chr" > ${i%.*}_chr.bg #	We need to change chromosome names to visualize the output in UCSC Browser http://seqanswers.com/forums/archive/index.php/t-22504.html

	LC_COLLATE=C sort -k1,1 -k2,2n -o ${i}.tmp ${i%.*}_chr.bg # Sort to proper order
	mv ${i}.tmp ${i%.*}_chr.bg
	$BEDGRAPH2BIGWIG ${i%.*}_chr.bg $SCRATCH/${GENOME}.fai.tmp ${i%.*}_chr.wg
#	rm ${i%.*}_chr.bg
done

rm $SCRATCH/${GENOME}.fai.tmp

mkdir -p $SCRATCH/star_signal/bg; mkdir $SCRATCH/star_signal/wg; mkdir $SCRATCH/star_signal/log
mv $SCRATCH/star_signal/*.bg $SCRATCH/star_signal/bg/; mv $SCRATCH/star_signal/*.wg $SCRATCH/star_signal/wg/; mv $SCRATCH/star_signal/*.out $SCRATCH/star_signal/log/

for i in $SCRATCH/star_signal/bg/*.bg # Compress .bg files
do 
	pigz $i
done &
####################################################################################################
### Cleaning results and indexing

# Move transcriptome mapping
mkdir $SCRATCH/transcriptome
mv $SCRATCH/*Transcriptome.out.bam $SCRATCH/transcriptome/

# Sort transcriptome BAMs
# Prepare for RSEM: sort transcriptome BAM to ensure the order of the reads, to make RSEM output (not pme) deterministic
cd $SCRATCH/transcriptome

for i in *Transcriptome.out.bam
do
	cat <( $SAMTOOLS view -H $i ) <( $SAMTOOLS view -@ $THREADS $i | sort -S ${MEM_LIMIT}G -T ./ ) | \
		$SAMTOOLS view -@ $THREADS -b - > $i.tmp
	mv $i.tmp $i
done

# Move STAR logs
mkdir $SCRATCH/star_log
mv $SCRATCH/*.out $SCRATCH/star_log/

$MULTIQC -o $SCRATCH/star_log $SCRATCH/star_log/ &

# Chimeric to BAM
cd $SCRATCH/

mkdir $SCRATCH/chimeric

for i in *.sam
do
	$SAMTOOLS view -@ $THREADS -b $i | $SAMTOOLS sort -@ $THREADS -T $SCRATCH/tmp.sort \
	-o $SCRATCH/chimeric/${i%.*}.bam -
	
	rm $i

	$SAMTOOLS index -@ $THREADS $SCRATCH/chimeric/${i%.*}.bam &
done

# Move STAR gene counts - should be mainly used for RSEM as stradness determination help/confirmation
mkdir -p $SCRATCH/star_gc
mv $SCRATCH/*ReadsPerGene.out.tab $SCRATCH/star_gc

$MULTIQC -o $SCRATCH/star_gc $SCRATCH/star_gc/
####################################################################################################
### Finalize and copy results

mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR_QC

mkdir -p $SCRATCH/alignment/genome
mv $SCRATCH/*.bam $SCRATCH/alignment/genome/
mv $SCRATCH/*.bai $SCRATCH/alignment/genome/
mv $SCRATCH/transcriptome $SCRATCH/alignment/
mkdir -p $SCRATCH/other/junction
mv $SCRATCH/*junction $SCRATCH/chimeric/
mv $SCRATCH/*SJ.out.tab $SCRATCH/other/junction/
mv $SCRATCH/chimeric $SCRATCH/other/

rm $SCRATCH/*$APPENDIX
rm $SCRATCH/*.fa.gz
rm $SCRATCH/$GENOME*
rm $SCRATCH/$GTF*
rm -r $SCRATCH/$GEN_DIR
rm -r $SCRATCH/*_STARpass1
rm -r $SCRATCH/*_STARgenome
rm -r $SCRATCH/*_STARtmp

cp -r $SCRATCH/star_signal $OUTPUT_DIR_QC/
cp -r $SCRATCH/star_log $OUTPUT_DIR_QC/
cp -r $SCRATCH/star_gc $OUTPUT_DIR_QC/
rm -r $SCRATCH/star_log $SCRATCH/star_gc $SCRATCH/star_signal
cp -r $SCRATCH/* $OUTPUT_DIR/

rm -r $SCRATCH/*
