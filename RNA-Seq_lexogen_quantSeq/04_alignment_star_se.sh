#!/bin/bash
#
# Alignment - Lexogen QuantSeq SE (RNA-Seq)
# Pipeline to align RNA-Seq experiment as recommended in GENCODE project https://github.com/ENCODE-DCC/long-rna-seq-pipeline/tree/master/dnanexus/align-star-pe; https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/dnanexus/Readme.md with few moditifications: --twopassMode Basic; --outSAMattributes All --outFilterMismatchNoverReadLmax 0.05
# To increase sensitivity you might try to add --seedSearchStartLmax 30 or even lower
#
# Should handle both Lexogen QuantSeq FWD (splicing should be on) and REV (splicing doesn't have to be on)
#
# Lexogen QuantSeq FWD is different than REV. For more details see https://www.lexogen.com/quantseq-data-analysis/ for FWD and https://www.lexogen.com/quantseq-data-analysis-rev/ for REV
# If we preprocess QuantSeq FWD it could be handled the same way as REV. REV doesn't require (most likely) splicing abnd two-pass mapping but it just takes more time so it could be included
#
# We are NOT using transcriptome mapping here, yet. Commands for transcriptome are commented out. This is the only difference between this alignment and general RNA-Seq alignment - only to save some time
#
# Requires STAR, pigz, samtools, multiQC
#
# qsub -l walltime=24:0:0 -q default -l select=1:ncpus=6:mem=36gb:scratch_local=150gb -I
#
# IN PROGRESS
#
############################################################################################
### Variables

INPUT_DIR=/storage/brno3-cerit/home/opplatek/biocore/karin_quantSeq/data/preprocessed
OUTPUT_DIR=/storage/brno3-cerit/home/opplatek/biocore/karin_quantSeq/results
OUTPUT_DIR_QC=${OUTPUT_DIR}/qc

#APPENDIX1="_R1_trim.fastq.gz"
#APPENDIX2="_R2_trim.fastq.gz"
APPENDIX="_trim.fastq.gz"

NO_MISMATCHES=999 # Number of mismatches per read - since we do not trim; 14 is recommended in QuantSeq FWD protocol for SE and 16 for PE; to turn this off set it to 999
PER_MISMATCHES=0.05 # Percent of mismatches per read; I usually keep this at 0.05 or 0.1; for untrimmed QuantSeq FWD it might be 12 + ((READ LENGTH-12) * 0.05) -> for 50 bp single end it should be at least 14 mismatches
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
SAMTOOLS=samtools
module add python27-modules-gcc
PYTHONPATH=$PYTHONPATH:/storage/brno2/home/opplatek/tools/MultiQC-1.0/lib/python2.7/site-packages/
MULTIQC=/storage/brno2/home/opplatek/tools/MultiQC-1.0/bin/multiqc

############################################################################################
### Copy inputs
cp $GTF $SCRATCH/

cp -r $GENOME_INDEX $SCRATCH/

cp $INPUT_DIR/*$APPENDIX $SCRATCH/

cd $SCRATCH/
###################################################################################################
### Genome and annotation preparation
GENOME=$(basename $GENOME)
GTF=$(basename $GTF)

unpigz -p $THREADS $GTF

GENOME=${GENOME%.gz*}
GTF=${GTF%.gz*}
GEN_DIR=${GENOME%.fa*}

# Set MAX ram
RAM=$[$RAM-1] # lower it by one to leave some space
MEM_LIMIT=$RAM # $[$RAM/$THREADS] # In case we need GB per thread

############################################################################################
### Alignment
for i in *$APPENDIX
do
	READ_FOR=$i
	echo "Now I am processing SE reads $READ_FOR - alignment"

	$STAR --runMode alignReads --runThreadN $THREADS --genomeDir $GEN_DIR \
	--readFilesIn $READ_FOR --readFilesCommand zcat --sjdbOverhang $RD_LENGTH --sjdbGTFfile $GTF \
	--outFileNamePrefix ${i%.fastq.gz*}.${GTF%.*} \
    --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax $NO_MISMATCHES --outFilterMismatchNoverReadLmax $PER_MISMATCHES \
    --alignIntronMin 20 --alignIntronMax $MAX_INTRON --alignMatesGapMax 1000000 \
	--outFilterMatchNmin 0 --outFilterMatchNminOverLread 0.66 \
    --outSAMheaderHD @HD VN:1.4 SO:coordinate --chimSegmentMin 30 --chimOutType SeparateSAMold \
    --outSAMunmapped Within --outFilterType BySJout --outSAMattributes All \
    --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts --sjdbScore 1 --twopassMode Basic --outMultimapperOrder Random # \ --outWigType bedGraph # --quantMode GeneCounts TranscriptomeSAM # Transcriptome mapping is not used here, yet
    
	echo "Done processing SE reads $READ_FOR - alignment"
done

############################################################################################
### Cleaning results and indexing

### Transcriptome mapping is not used here, yet ###
# Move transcriptome mapping
#mkdir $SCRATCH/transcriptome
#mv $SCRATCH/*Transcriptome.out.bam $SCRATCH/transcriptome/

# Sort transcriptome BAMs
#cd $SCRATCH/transcriptome
#
#for i in *Transcriptome.out.bam
#do
#	cat <( $SAMTOOLS view -H $i ) \
#	    <( $SAMTOOLS view -@ $THREADS $i | sort -S ${MEM_LIMIT}G -T ./ ) | \
#	    $SAMTOOLS view -@ $THREADS -b - > $i.tmp
#	mv $i.tmp $i
#done
### Transcriptome mapping is not used here, yet ###

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
#mv $SCRATCH/transcriptome $SCRATCH/alignment/ # Transcriptome mapping is not used here, yet
mkdir $SCRATCH/other
mv $SCRATCH/chimeric $SCRATCH/other/
mv $SCRATCH/*junction $SCRATCH/other/
mv $SCRATCH/*SJ.out.tab $SCRATCH/other/

rm $SCRATCH/*$APPENDIX
rm $SCRATCH/*.fa.gz
rm $SCRATCH/$GENOME*
rm $SCRATCH/$GTF*
rm -r $SCRATCH/$GEN_DIR
rm -r $SCRATCH/*_STARpass1
rm -r $SCRATCH/*_STARgenome
rm -r $SCRATCH/*_STARtmp

cp -r $SCRATCH/star_log $OUTPUT_DIR_QC/
cp -r $SCRATCH/star_gc $OUTPUT_DIR_QC/
rm -r $SCRATCH/star_log $SCRATCH/star_gc
cp -r $SCRATCH/* $OUTPUT_DIR/

rm -r $SCRATCH/*
