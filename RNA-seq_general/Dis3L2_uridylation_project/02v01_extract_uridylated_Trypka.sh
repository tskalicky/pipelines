#!/bin/bash
## Variables
THREADS="10"
INPUT_DIR="/home/tomas/Bioinformatics/CEITEC/Dis3L2/Dasa_spikein/uridylation/fixed_pairs_fastq"
OUTPUT_DIR="/home/tomas/Bioinformatics/CEITEC/Dis3L2/Dasa_spikein/uridylation/sorted_uridyl/fastq"
EXTRACT_FASTA="/home/tomas/Bioinformatics/CEITEC/Dis3L2/Dasa_spikein/uridylation/extract_uridylated_out_fasta.py"
EXTRACT="/home/tomas/Bioinformatics/CEITEC/Dis3L2/Dasa_spikein/uridylation/fixed_pairs_fastq/extract_uridylated.py"
#
APPENDIX1=".bam"
APPENDIX2=".fastq"
APPENDIX3=".fa"
APPENDIX4=".fa.gz"
# Binaries
SAMTOOLS=$(which samtools)
MULTIQC=$(which multiqc)
# Check the tools versions
which $SAMTOOLS
which $MULTIQC
which python
#
####################################################################################################
## Commands
cd $INPUT_DIR
mkdir -p $OUTPUT_DIR
find "$INPUT_DIR" -maxdepth 1 -name "*.fastq.gz" -exec cp -vt "$OUTPUT_DIR" {} +
cp -av $EXTRACT $OUTPUT_DIR
cd $OUTPUT_DIR
set +f # enable globbing, because is needed for other parts of the script, like variable expansions !!!
for a in *.fastq.gz
do
	SAMPLE=$a
	SAMPLENAME=${a%.*.*}
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Start unpacking reads from file $SAMPLE"
	unpigz -vp $THREADS $SAMPLE
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Start extracting uridylated reads from file $SAMPLE"
	python extract_uridylated.py $SAMPLENAME".fastq" &
done
#wait
wait
# date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# echo "Start packing output files"
# pigz -vp $THREADS *.fastq
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Script finished."


