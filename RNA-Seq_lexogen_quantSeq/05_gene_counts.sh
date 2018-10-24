#!/bin/bash
#PBS -l select=1:ncpus=6:mem=25gb:scratch_local=200gb
#PBS -l walltime=12:00:00
#PBS -N 05_gene_counts
#
# Gene counting - Lexogen QuantSeq SE (RNA-Seq) 
#
# Can handle both SE and PE but you should now the strandness of the experiment (either by kit name or by Picard)
# Input is from STAR mapping to genome and transcriptome (genome for featureCounts and transcriptome for RSEM); RSEM is NOT used here yet
# featureCounts should be prefered option (for QuantSeq)
#
# Should handle both Lexogen QuantSeq FWD (counting by exons) and REV (counting by 3'UTRs)
#
# We are NOT using transcriptome mapping here, yet. Commands for transcriptome are commented out.
#
# Requires featureCounts, multiQC
#
# IN PROGRESS
#
############################################################################################
### Variables
INPUT_DIR=/storage/brno3-cerit/home/opplatek/biocore/karin_quantSeq/results/alignment
OUTPUT_DIR=/storage/brno3-cerit/home/opplatek/biocore/karin_quantSeq/results/gene_counts
OUTPUT_DIR_QC=/storage/brno3-cerit/home/opplatek/biocore/karin_quantSeq/results/qc

COUNT_OVER="exon" # [exon, three_prime_utr] what feature to use for the summary? For QuantSeq FWD is should be "exon", for QuantSeq REV it should be "three_prime_utr"
PAIRED="false" # [true, false] True or false for PE reads 
STRANDNESS="rev" # [fwd, rev, none] strandedness of read. For QuantSeq FWD is should be "rev" 

GENOME=/storage/brno2/home/opplatek/genomes/human/ensembl87/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
GTF=/storage/brno2/home/opplatek/genomes/human/ensembl87/Homo_sapiens.GRCh38.87.gtf.gz
#RSEM_GENOME_INDEX=/storage/brno2/home/opplatek/genomes/human/ensembl87/RSEM_index

THREADS=$PBS_NUM_PPN

# Binaries
FEATURE_COUNTS=/storage/brno2/home/opplatek/tools/subread-1.5.2-Linux-x86_64/bin/featureCounts
#RSEM_BIN=/storage/brno2/home/opplatek/tools/RSEM-1.3.0 # Transcriptome mapping is not used here, yet
module add samtools-1.4
SAMTOOLS=samtools
module add python27-modules-gcc
PYTHONPATH=$PYTHONPATH:/storage/brno2/home/opplatek/tools/MultiQC-1.0/lib/python2.7/site-packages/
MULTIQC=/storage/brno2/home/opplatek/tools/MultiQC-1.0/bin/multiqc

############################################################################################
### Copy inputs
cp $GTF $SCRATCH/

#cp -r $RSEM_GENOME_INDEX $SCRATCH/ # Transcriptome mapping is not used here, yet

cp -r $INPUT_DIR/* $SCRATCH/ 

cd $SCRATCH/

####################################################################################################
### Genome and annotation preparation
GENOME=$(basename $GENOME) 
GTF=$(basename $GTF)

unpigz -p $THREADS $GTF

GENOME=${GENOME%.gz*}
GTF=${GTF%.gz*}
GEN_DIR=${GENOME%.fa*}

# Set extra flags for featureCounts and RSEM
msg=""
extra_flags_feature=""
#extra_flags_rsem="" # Transcriptome mapping is not used here, yet
if [ "$paired_end" == "true" ]; then
	msg="PE experiment"
	extra_flags_feature="-p" # For featureCounts
#	extra_flags_rsem="--paired-end" # for RSEM # Transcriptome mapping is not used here, yet
else
	msg="SE experiment"
fi
# Set if it is stranded or not (fwd, rev, none)
if [ "$STRANDNESS" == "fwd" ] || [ "$STRANDNESS" == "+" ] || [ "$STRANDNESS" == "ScriptSeq" ]; then
	extra_flags_feature="$extra_flags_feature -s 1" # For featureCounts
#	extra_flags_rsem="$extra_flags_rsem --forward-prob 1" # for RSEM # Transcriptome mapping is not used here, yet
	msg="$msg forward strand"
elif [ "$STRANDNESS" == "rev" ] || [ "$STRANDNESS" == "-" ] || [ "$STRANDNESS" == "TruSeq" ]; then
	extra_flags_feature="$extra_flags_feature -s 2" # For featureCounts
#	extra_flags_rsem="$extra_flags_rsem --forward-prob 0" # for RSEM # Transcriptome mapping is not used here, yet
	msg="$msg reverse strand"
else
	extra_flags_feature="$extra_flags_feature -s 0" # For featureCounts
#	extra_flags_rsem="$extra_flags_rsem --forward-prob 0.5" # for RSEM # Transcriptome mapping is not used here, yet
	msg="$msg unstranded"
fi

echo "Running as $msg with counts over $COUNT_OVER"

############################################################################################
### Gene/isoform counting
# Count reads using featureCounts
mkdir -p $SCRATCH/gene_counts/featureCounts

cd $SCRATCH/genome/
$FEATURE_COUNTS -t $COUNT_OVER -g gene_id $extra_flags_feature -T $THREADS -a $SCRATCH/$GTF -o $SCRATCH/gene_counts/featureCounts/${STRANDNESS}.featureCounts *.bam

mkdir $SCRATCH/featureCounts_gc

$MULTIQC -o $SCRATCH/featureCounts_gc $SCRATCH/gene_counts/featureCounts/

### Transcriptome mapping is not used here, yet ###
# Count reads using RSEM
#mkdir -p $SCRATCH/gene_counts/rsem

#cd $SCRATCH/transcriptome/
#for i in *.bam
#do
#	echo "Started RSEM counting $i"
#	$RSEM_BIN/rsem-calculate-expression --bam --estimate-rspd --calc-ci --seed $RANDOM -p $THREADS \
#    --no-bam-output --ci-memory 30000 ${extra_flags_rsem} $i $SCRATCH/RSEM_index/$GEN_DIR $SCRATCH/gene_counts/rsem/${i%.*}.rsem
#	echo "Done RSEM counting $i"
#done
### Transcriptome mapping is not used here, yet ###

############################################################################################
### Copying results
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR_QC

cp -r $SCRATCH/gene_counts/* $OUTPUT_DIR/
cp -r $SCRATCH/featureCounts_gc $OUTPUT_DIR_QC/

rm -r $SCRATCH/*
