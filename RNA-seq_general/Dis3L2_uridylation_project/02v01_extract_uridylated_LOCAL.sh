#!/bin/bash
## Variables
THREADS="4"
INPUT_DIR="/home/tomas/ownCloud/CEITEC_lab/Dis3L2/Dasa_spikein/uridylation/extract_mapped_unmapped_from_orig_mapped/fixed_pairs_fastq/"
OUTPUT_DIR="/home/tomas/ownCloud/CEITEC_lab/Dis3L2/Dasa_spikein/uridylation/sorted_uridyl/fastq/"
EXTRACT_FASTA="/home/tomas/ownCloud/git/pipelines/RNA-seq_general/Dis3L2_uridylation_project/extract_uridylated_out_fasta.py"
EXTRACT="/home/tomas/ownCloud/git/pipelines/RNA-seq_general/Dis3L2_uridylation_project/extract_uridylated.py"
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
# cd $INPUT_DIR
# mkdir -p $OUTPUT_DIR
# find "$INPUT_DIR" -maxdepth 1 -name "*.fastq.gz" -exec cp -vt "$OUTPUT_DIR" {} +
# cp -av $EXTRACT $OUTPUT_DIR
cd $OUTPUT_DIR
set +f # enable globbing, because is needed for other parts of the script, like variable expansions !!!
for a in *.fastq
do
	SAMPLE=$a
	SAMPLENAME=${a%.*.*}
	# date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	# echo "Start unpacking reads from file $SAMPLE"
	# unpigz -vp $THREADS $SAMPLE
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Start extracting uridylated reads from file $SAMPLE"
	python extract_uridylated.py $SAMPLE
done
#wait
wait
# date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# echo "Start packing output files"
# pigz -vp $THREADS *.fastq
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Script finished."


