#!/bin/bash
## Variables
THREADS="3"
INPUT_DIR="/storage/brno3-cerit/home/tskalicky/Dis3L2/Dasa_spikein/mapping/human_genome/star/results_job_1052371/alignment/genome"
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/Dis3L2/Dasa_spikein/mapping/human_genome/star/results_job_1052371/alignment/genome/mapped_extracted"
EXTRACT="/home/skalda/ownCloud/git/pipelines/RNA-seq_general/Dis3L2_uridylation_project/extract_uridylated.py"
#
APPENDIX1=".bam"
APPENDIX2=".fastq"
APPENDIX3=".fa"
# Binaries
SAMTOOLS=$(which samtools)
MULTIQC=$(which multiqc)
# Check the tools versions
which $SAMTOOLS

####################################################################################################
## Commands
cd $INPUT_DIR
mkdir -p $OUTPUT_DIR
for a in .fa.gz
do
	SAMPLE=$a
	SAMPLENAME=${a%.*.*}
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Start unpacking reads from file $SAMPLE"
	unpigz -vp $THREADS $SAMPLE
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Start extracting uridylated reads from file $SAMPLE"
	python extract_uridylated_out_fasta.py $SAMPLENAME".fa"
	pigz -vp $THREADS "*.fa"
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Start packing output files from file $SAMPLE"
done
wait
