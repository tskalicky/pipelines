#!/bin/bash
#PBS -l walltime=96:0:0 
#PBS -q default@wagap-pro.cerit-sc.cz 
#PBS -l select=1:ncpus=10:mem=50gb:scratch_local=250gb
#PBS -j oe
#PBS -N 05_gene_counts
#
## initialize the required application
module add samtools-1.6
module add python27-modules-intel #required by multiQC
module add rsem-1.2.8
#
# Gene/isoform counting - RNA-Seq
#
# Can handle both SE and PE but you should now the strandness of the experiment (either by kit name or by Picard)
# Input is from STAR mapping to genome and transcriptome (genome for featureCounts and transcriptome for RSEM)
# RSEM should be prefered option
#
# Requires Samtools, MultiQC, featureCounts, RSEM
#
############################################################################################
### Variables
INPUT_DIR=/mnt/nfs/home/323639/999993-Bioda/projects/honza/smida_rnaseq/results/alignment # Pointing to general mapping directory containing genome/ and transcriptome/ mapping
OUTPUT_DIR=/mnt/nfs/home/323639/999993-Bioda/projects/honza/smida_rnaseq
OUTPUT_DIR=${OUTPUT_DIR}/results
OUTPUT_DIR_QC=${OUTPUT_DIR}/qc
OUTPUT_DIR=${OUTPUT_DIR}/gene_counts

COUNT_OVER="exon" # [exon, three_prime_utr] what feature to use for the summary? For QuantSeq it might be 3 UTR ("three_prime_utr" is for Ensembl annotation
PAIRED="false" # [true, false] "true" for PE reads, "false" for SE reads
STRANDNESS="rev" # [fwd, rev, none] strandedness of read 

RSEM_GENOME_INDEX="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/RSEM_index" # Reference RSEM index
MY_GENOME="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" # Reference fasta (.gz)
MY_GTF="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/Homo_sapiens.GRCh38.91.gtf.gz" # Reference GTF (.gz)

THREADS=$PBS_NUM_PPN

RSEM_RANDOM=12345

# Binaries
FEATURE_COUNTS="/storage/brno3-cerit/home/tskalicky/home/tskalicky/tools/subread-1.6.0/bin/featureCounts"
RSEM_BIN="/storage/brno3-cerit/home/tskalicky/tools/RSEM-1.3.0/bin"
SAMTOOLS=$(which samtools)
MULTIQC=$(which multiqc)

# Check the tools versions
which $FEATURE_COUNTS
which $RSEM
which $SAMTOOLS
which $MULTIQC

############################################################################################
### Copy inputs
cd $INPUT_DIR
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
cp -av $MY_GENOME $MY_GTF $SCRATCHDIR
cp -avr $RSEM_GENOME_INDEX $INPUT_DIR $SCRATCHDIR/ 
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR
echo "Following reads were copied to scratch:"
ls -c1

####################################################################################################
### Genome and annotation preparation
GENOME_NAME=$(basename $MY_GENOME)
GTF_NAME=$(basename $MY_GTF)
#GEN_DIR=$(basename $RSEM_GENOME_INDEX)

unpigz -p $THREADS $GENOME_NAME
unpigz -p $THREADS $GTF_NAME

GENOME=${GENOME_NAME%.*.*}
GTF=${GTF_NAME%.*.*}
GEN_DIR=$${GENOME_NAME%.*.*}

# Set extra flags for featureCounts and RSEM
msg=""
extra_flags_feature=""
extra_flags_rsem=""
if [ "$PAIRED" == "true" ]; then
	msg="PE experiment"
	extra_flags_feature="-p -C" # For featureCounts; -B
	extra_flags_rsem="--paired-end" # for RSEM
else
	msg="SE experiment"
fi
# Set if it is stranded or not (fwd, rev, none)
if [ "$STRANDNESS" == "fwd" ] || [ "$STRANDNESS" == "+" ] || [ "$STRANDNESS" == "ScriptSeq" ]; then
	extra_flags_feature="$extra_flags_feature -s 1" # For featureCounts
	extra_flags_rsem="$extra_flags_rsem --forward-prob 1" # for RSEM
	msg="$msg forward strand"
elif [ "$STRANDNESS" == "rev" ] || [ "$STRANDNESS" == "-" ] || [ "$STRANDNESS" == "TruSeq" ]; then
	extra_flags_feature="$extra_flags_feature -s 2" # For featureCounts
	extra_flags_rsem="$extra_flags_rsem --forward-prob 0" # for RSEM
	msg="$msg reverse strand"
else
	extra_flags_feature="$extra_flags_feature -s 0" # For featureCounts
	extra_flags_rsem="$extra_flags_rsem --forward-prob 0.5" # for RSEM
	msg="$msg unstranded"
fi
echo "Running as $msg"

############################################################################################
### Gene/isoform counting
# Count reads using featureCounts
mkdir -p $SCRATCHDIR/gene_counts/featureCounts

cd $SCRATCHDIR/genome/
$FEATURE_COUNTS -t $COUNT_OVER -g gene_id $extra_flags_feature -T $THREADS -a $SCRATCHDIR/$GTF -o $SCRATCHDIR/gene_counts/featureCounts/${STRANDNESS}.featureCounts *.bam # Be careful about --primary - only manual says "If specified, only primary alignments will be counted. Primary and secondary alignments are identified using bit 0x100 in the Flag field of SAM/BAM files. All primary alignments in a dataset will be counted no matter they are from multimapping reads or not (ie. ‘-M’ is ignored)." The second part is not mentioned in command line version!

mkdir $SCRATCHDIR/featureCounts_gc

$MULTIQC -o $SCRATCHDIR/featureCounts_gc $SCRATCHDIR/gene_counts/featureCounts/ &

# Count reads using RSEM
mkdir -p $SCRATCHDIR/gene_counts/rsem

cd $SCRATCHDIR/transcriptome/

#RSEM_RANDOM=$RANDOM
echo "Using $RSEM_RANDOM as --seed"

for i in *.bam
do
	echo "Started RSEM counting $i"
	$RSEM --bam --estimate-rspd --calc-ci --seed $RSEM_RANDOM -p $THREADS --no-bam-output --ci-memory 30000 ${extra_flags_rsem} $i $SCRATCHDIR/RSEM_index/$GEN_DIR $SCRATCHDIR/gene_counts/rsem/${i%.*}.rsem
	echo "Done RSEM counting $i"
done

mkdir $SCRATCHDIR/RSEM_gc

$MULTIQC -o $SCRATCHDIR/RSEM_gc $SCRATCHDIR/gene_counts/rsem/ &

for counts in $SCRATCHDIR/gene_counts/rsem/*.results
do
	pigz -p $THREADS $counts
done

wait

############################################################################################
### Copying results
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR_QC

cp -r $SCRATCHDIR/gene_counts/* $OUTPUT_DIR/
cp -r $SCRATCHDIR/featureCounts_gc $OUTPUT_DIR_QC/
cp -r $SCRATCHDIR/RSEM_gc $OUTPUT_DIR_QC/

rm -r $SCRATCHDIR/*