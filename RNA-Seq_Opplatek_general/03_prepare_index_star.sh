#!/bin/bash
#PBS -l walltime=96:0:0 
#PBS -q uv@wagap-pro.cerit-sc.cz 
#PBS -l select=1:ncpus=1:mem=40gb:scratch_local=80gb
#PBS -j oe
#PBS -N 03v1_RSEM_index
#
# Prepare STAR alignment index
#
# Requires STAR, RSEM, pigz, bowtie or alernatively STAR aligner
#
## initialize the required application
module add rsem-1.2.8
module add star-2.5.2b
module add bowtie-1.0.0
module add bowtie2-2.3.0
############################################################################################
### Variables
DATADIR="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91"
GENOME="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GTF="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/Homo_sapiens.GRCh38.91.gtf.gz"
INDEX_OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/RSEM_index"

THREADS=$PBS_NUM_PPN

# Binaries
RSEM_BIN="/storage/brno3-cerit/home/tskalicky/tools/RSEM-1.3.0/bin"

RSEM-PREP=$(which rsem-prepare-reference)
STAR=$(which STAR)
BOWTIE2="$(which bowtie2)"
# Check the tools versions
echo $FASTQC
echo $FLEXBAR
echo $TRIMMOMATIC
#
# Prepare STAR alignment index
#
# Requires STAR, pigz
#
# "For small genomes, the parameter --genomeSAindexNbases must to be scaled down, with a typical
#       value of min(14, log2(GenomeLength)/2 - 1). For example, for 1 megaBase genome, this is equal
#       to 9 (8.965784), for 100 kiloBase genome, this is equal to 7."
############################################################################################
### Variables
GENOME=/storage/brno2/home/opplatek/genomes/athaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
GTF=/storage/brno2/home/opplatek/genomes/athaliana/Arabidopsis_thaliana.TAIR10.37.gtf.gz
INDEX_OUTPUT_DIR=/storage/brno2/home/opplatek/genomes/athaliana/ensembl37/TAIR10/STAR_index # Where you want to copy created index?

# Read length for sjdbOverhang; --sjdbOverhang 100 should work fine for most of the data, but more specific setting based on the real read length should be more sensitive https://groups.google.com/forum/#!msg/rna-star/h9oh10UlvhI/BfSPGivUHmsJ
RD_LENGTH=100 # Read length from the sequencing. Illumina sometimes reports N+1 http://seqanswers.com/forums/archive/index.php/t-31154.html; in case you change this value uncomment next line as well
#RD_LENGTH=$[$RD_LENGTH-$TRIM_LEFT-$TRIM_RIGHT-1] # Modified read length for the index creation and mapping - should be read length -1

THREADS=$PBS_NUM_PPN # Number of threads to use

# Binaries
STAR_BIN=/storage/brno2/home/opplatek/tools/STAR-2.5.2b/bin/Linux_x86_64

############################################################################################
### Copy inputs
cp $GTF $SCRATCH/
cp $GENOME $SCRATCH/

cd $SCRATCH/

####################################################################################################
### Genome and annotation preparation
GENOME=$(basename $GENOME)
GTF=$(basename $GTF)

unpigz -p $THREADS $GENOME
unpigz -p $THREADS $GTF

GENOME=${GENOME%.gz*}
GTF=${GTF%.gz*}
GEN_DIR=${GENOME%.fa*}

############################################################################################
### Indexing - RUN ONLY ONCE FOR ONE GENOME AND STAR VERSION AND READ LENGTH OR USE DEFAULT INDEX
mkdir $SCRATCH/$GEN_DIR

$STAR_BIN/STAR --runMode genomeGenerate --runThreadN $THREADS --genomeDir $SCRATCH/$GEN_DIR --genomeFastaFiles $GENOME \
--sjdbGTFfile $GTF --sjdbOverhang $RD_LENGTH # --sjdbOverhang 100 should work fine for most of the data, but more specific setting based on the real read length should be more sensitive https://groups.google.com/forum/#!msg/rna-star/h9oh10UlvhI/BfSPGivUHmsJ

############################################################################################
### Copying results
mkdir -p $INDEX_OUTPUT_DIR
cp -r $SCRATCH/$GEN_DIR $INDEX_OUTPUT_DIR/

rm -r $SCRATCH/*
