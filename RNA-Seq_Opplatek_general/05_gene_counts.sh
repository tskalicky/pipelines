#!/bin/bash
#PBS -l walltime=96:0:0 
#PBS -q default@wagap-pro.cerit-sc.cz 
#PBS -l select=1:ncpus=10:mem=50gb:scratch_local=250gb:os=debian9
#PBS -j oe
#PBS -N 05_gene_counts
#
## initialize the required application
module add samtools-1.6
module add python27-modules-intel #required by multiQC
module add piranha-1.2.1 #required by Piranha
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
### Variables
INPUT_DIR=/mnt/nfs/home/323639/999993-Bioda/projects/honza/smida_rnaseq/results/alignment # Pointing to general mapping directory containing genome/ and transcriptome/ mapping
OUTPUT_DIR=/mnt/nfs/home/323639/999993-Bioda/projects/honza/smida_rnaseq
OUTPUT_DIR=${OUTPUT_DIR}/results
OUTPUT_DIR_QC=${OUTPUT_DIR}/qc
OUTPUT_DIR=${OUTPUT_DIR}/gene_counts

COUNT_OVER="exon" # [exon, three_prime_utr] what feature to use for the summary? For QuantSeq it might be 3 UTR ("three_prime_utr" is for Ensembl annotation
PAIRED="false" # [true, false] "true" for PE reads, "false" for SE reads
STRANDNESS="rev" # [fwd, rev, none] strandedness of read 

GENOME=/storage/brno2/home/opplatek/genomes/human/ensembl87/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz # Reference fasta (.gz)
GTF=/storage/brno2/home/opplatek/genomes/human/ensembl87/Homo_sapiens.GRCh38.87.gtf.gz # Reference GTF (.gz)
RSEM_GENOME_INDEX=/storage/brno2/home/opplatek/genomes/human/ensembl87/RSEM_index

THREADS=$PBS_NUM_PPN

RSEM_RANDOM=12345

# Binaries
FEATURE_COUNTS=/storage/brno2/home/opplatek/tools/subread-1.5.2-Linux-x86_64/bin/featureCounts
RSEM=/storage/brno2/home/opplatek/tools/RSEM-1.3.0/rsem-calculate-expression
module add samtools-1.4
SAMTOOLS=$(which samtools)
module add python27-modules-gcc
#PYTHONPATH=$PYTHONPATH:/storage/brno2/home/opplatek/tools/MultiQC-1.0/lib/python2.7/site-packages/
#MULTIQC=/storage/brno2/home/opplatek/tools/MultiQC-1.0/bin/multiqc
MULTIQC=$(which multiqc)

# Check the tools versions
which $FEATURE_COUNTS
which $RSEM
which $SAMTOOLS
which $MULTIQC

############################################################################################
### Copy inputs
cp $GTF $SCRATCH/
#cp $GENOME $SCRATCH/
cp -r $RSEM_GENOME_INDEX $SCRATCH/

cp -r $INPUT_DIR/* $SCRATCH/ 

cd $SCRATCH/

####################################################################################################
### Genome and annotation preparation
GENOME=$(basename $GENOME)
GTF=$(basename $GTF)

#unpigz -p $THREADS $GENOME
unpigz -p $THREADS $GTF

GENOME=${GENOME%.gz*}
GTF=${GTF%.gz*}
GEN_DIR=${GENOME%.fa*}

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
mkdir -p $SCRATCH/gene_counts/featureCounts

cd $SCRATCH/genome/
$FEATURE_COUNTS -t $COUNT_OVER -g gene_id $extra_flags_feature -T $THREADS -a $SCRATCH/$GTF -o $SCRATCH/gene_counts/featureCounts/${STRANDNESS}.featureCounts *.bam # Be careful about --primary - only manual says "If specified, only primary alignments will be counted. Primary and secondary alignments are identified using bit 0x100 in the Flag field of SAM/BAM files. All primary alignments in a dataset will be counted no matter they are from multimapping reads or not (ie. ‘-M’ is ignored)." The second part is not mentioned in command line version!

mkdir $SCRATCH/featureCounts_gc

$MULTIQC -o $SCRATCH/featureCounts_gc $SCRATCH/gene_counts/featureCounts/ &

# Count reads using RSEM
mkdir -p $SCRATCH/gene_counts/rsem

cd $SCRATCH/transcriptome/

#RSEM_RANDOM=$RANDOM
echo "Using $RSEM_RANDOM as --seed"

for i in *.bam
do
	echo "Started RSEM counting $i"
	$RSEM --bam --estimate-rspd --calc-ci --seed $RSEM_RANDOM -p $THREADS --no-bam-output --ci-memory 30000 ${extra_flags_rsem} $i $SCRATCH/RSEM_index/$GEN_DIR $SCRATCH/gene_counts/rsem/${i%.*}.rsem
	echo "Done RSEM counting $i"
done

mkdir $SCRATCH/RSEM_gc

$MULTIQC -o $SCRATCH/RSEM_gc $SCRATCH/gene_counts/rsem/ &

for counts in $SCRATCH/gene_counts/rsem/*.results
do
	pigz -p $THREADS $counts
done

wait

############################################################################################
### Copying results
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR_QC

cp -r $SCRATCH/gene_counts/* $OUTPUT_DIR/
cp -r $SCRATCH/featureCounts_gc $OUTPUT_DIR_QC/
cp -r $SCRATCH/RSEM_gc $OUTPUT_DIR_QC/

rm -r $SCRATCH/*