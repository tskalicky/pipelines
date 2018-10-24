#!/bin/bash
#PBS -l select=1:ncpus=12:mem=10gb:scratch_local=200gb
#PBS -l walltime=24:00:00
#PBS -q default
#PBS -N 06_fastq_screen
#
# Estimate rRNA, tRNA and other selected references/contamination content
# An alternative to Picard
# It could be run right after read preprocessing
#
# In default, fastq_screen ignores reads shorter than 20 bp - need to modify fastq_screen source at line 
#	~320 from 20bp->15bp to handle miRNA data
#
# Requires: fastq-screen; GD::Graph::bars (Perl); bowtie; bowtie2; multiqc; fastx-toolkit
#
####################################################################################################
### Variables
INPUT_DIR=/storage/brno2/home/opplatek/bioinf_projects/biocore/anzer_rnaseq/data/preprocessed
OUTPUT_DIR=/storage/brno2/home/opplatek/bioinf_projects/biocore/anzer_rnaseq
OUTPUT_DIR_QC=$OUTPUT_DIR/results/qc

APPENDIX="_trim.fastq.gz"

RRNA_INDEX=/storage/brno2/home/opplatek/genomes/drosophila/fastq_screen/index/bowtie2/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff.rRNA.fasta # Pointing to rRNA .bt2/.ebwt files
TRNA_INDEX=/storage/brno2/home/opplatek/genomes/drosophila/fastq_screen/index/bowtie2/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff.tRNA.fasta # Pointing to tRNA .bt2/.ebwt files
GENOME_INDEX=/storage/brno2/home/opplatek/genomes/drosophila/fastq_screen/index/bowtie2/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna # Pointing to genome .bt2/.ebwt files

FASTQ_SCREEN_CONF=/storage/brno2/home/opplatek/bioinf_projects/biocore/anzer_rnaseq/scripts/fastq_screen.conf # Config file that has to be modified to include references and bowtie versions/locations

FASTQ_SCREEN=/storage/brno2/home/opplatek/tools/fastq_screen_v0.11.3/fastq_screen
module rm bowtie2-2.3.0
module add bowtie2-2.3.0
module add bowtie-1.0.0
BOWTIE_BIN=/afs/ics.muni.cz/software/bowtie-1.0.0/bin
BOWTIE2_BIN=/software/bowtie/2.3.0/bin
FASTQSCREEN_R=/storage/brno2/home/opplatek/tools/scripts/fastqscreen.r # Preseq visualization script
module add fastx-0.0.14
FASTX_TOOLKIT=/software/fastx/0.0.14/bin
module add python27-modules-gcc
#PYTHONPATH=$PYTHONPATH:/storage/brno2/home/opplatek/tools/MultiQC-1.0/lib/python2.7/site-packages/
#MULTIQC=/storage/brno2/home/opplatek/tools/MultiQC-1.0/bin/multiqc
MULTIQC=$(which multiqc)

THREADS=$PBS_NUM_PPN

# Check binaries
which $FASTQ_SCREEN
which $BOWTIE_BIN
which $BOWTIE_BIN2
which $FASTQSCREEN_R
which $FASTX_TOOLKIT
which $MULTIQC

############################################################################################
### Copy inputs
cp $INPUT_DIR/*$APPENDIX $SCRATCH/ &

cp $GENOME_INDEX* $SCRATCH/ &
cp $RRNA_INDEX* $SCRATCH/ &
cp $TRNA_INDEX* $SCRATCH/ &
cp $FASTQ_SCREEN_CONF $SCRATCH/ &

#cp $GENOME $SCRATCH/
#cp $RRNA $SCRATCH/
#cp $TRNA $SCRATCH/

FASTQ_SCREEN_CONF=$(basename $FASTQ_SCREEN_CONF)

wait 

####################################################################################################
### Contamination check, rRNA/tRNA content check
cd $SCRATCH/

# Install cpanm and module for plotting - fastq_screen works even without this, so it's optional
# Possible source of problems,  might not work
curl -L https://cpanmin.us | perl - App::cpanminus
mkdir $HOME/perl5
~/perl5/bin/cpanm --local-lib=$HOME/perl5 local::lib && eval $(perl -I $HOME/perl5/lib/perl5/ -Mlocal::lib)
#export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/usr/lib/x86_64-linux-gnu/pkgconfig
~/perl5/bin/cpanm --local-lib=$HOME/perl5 --force GD::Graph::bars # For fastqc plotting

### Launch fastqc_screen
mkdir $SCRATCH/fastq_screen

for i in *$APPENDIX
do
	$FASTQ_SCREEN --subset 500000 --outdir $SCRATCH/fastq_screen --threads $THREADS --conf $FASTQ_SCREEN_CONF --nohits --aligner bowtie2 --force $i # Subset just 500k reads, to use all reads put there 0
done

cd $SCRATCH/fastq_screen/

# Will not work if GDlib is not installed
convert *.png all_screen.pdf & # Merge .png outputs together
$MULTIQC .

module add R-3.4.0-gcc
R --no-save < $FASTQSCREEN_R
module rm R-3.4.0-gcc

# Get the most frequent unampped reads - merge them and collapse them
# Each separately
for i in *.tagged_filter.fastq.gz
do
	unpigz -c -p $THREADS $i | $FASTX_TOOLKIT/fastx_collapser -Q33 | pigz -c -p $THREADS - > ${i%.fastq*}.collapsed.fasta.gz
	unpigz -c -p $THREADS ${i%.fastq*}.collapsed.fasta.gz | head -20 > ${i%.fastq*}.collapsed.top10.fasta
done

# All of them together
unpigz -c -p $THREADS *.tagged_filter.fastq.gz | $FASTX_TOOLKIT/fastx_collapser -Q33 | pigz -c -p $THREADS - > all.collapsed.fasta.gz
unpigz -c -p $THREADS all.collapsed.fasta.gz | head -20 > all.collapsed.top10.fasta

############################################################################################
### Copy outputs
mkdir -p $OUTPUT_DIR_QC
cp -r $SCRATCH/fastq_screen $OUTPUT_DIR_QC/

rm -r $SCRATCH/*
