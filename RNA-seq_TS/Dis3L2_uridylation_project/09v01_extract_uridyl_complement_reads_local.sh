#!/bin/bash
############################################################################################
# Uridylated RNAseq reads processing
# Requires samtools and bbmap
############################################################################################
### Variables
THREADS=$PBS_NUM_PPN
INPUT_DIR="/home/tomas/CEITEC_lab/Dis3L2/Dasa_spikein/data/trimmed"
URI_EXTRACT="/home/tomas/ownCloud/git/pipelines/RNA-seq_general/Dis3L2_uridylation_project/extract_uri_complement_AAA.py"
OUTPUT_DIR="/home/tomas/ownCloud/CEITEC_lab/Dis3L2/Dasa_spikein/uridylation/complement_reads"
EXTRACT_FASTA="/storage/brno3-cerit/home/tskalicky/skripty/pipelines/RNA-seq_general/Dis3L2_uridylation_project/extract_uridylated_out_fasta.py"
EXTRACT="/storage/brno3-cerit/home/tskalicky/skripty/pipelines/RNA-seq_general/Dis3L2_uridylation_project/extract_uridylated.py"
export PATH="$PATH:/home/users/tskalicky/tools/bbmap/"
echo "PATH after modification is:"
echo '$PATH | tr ":" "\n" | nl'
#
APPENDIX1=".bam"
APPENDIX2=".fastq"
APPENDIX3=".fa"
# Binaries
SAMTOOLS=$(which samtools)
REPAIR=$(which repair.sh)
BBSPLITPAIRS=$(which bbsplitpairs.sh)
MULTIQC=$(which multiqc)
# Check the tools versions
which $SAMTOOLS
which $MULTIQC
which $REPAIR
which $BBSPLITPAIRS

####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
cd $INPUT_DIR
cp -av $URI_EXTRACT $INPUT_DIR
for a in *$APPENDIX2
do
	SAMPLE=$a
	SAMPLENAME=${a%.*}
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Start extracting complement uridylated reads from file $SAMPLE"
	python extract_uri_complement_AAA.py $SAMPLE
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Finnished extracting complement uridylated reads from file $SAMPLE"
done
#
wait
#
for b in "Uridyl_compl_"*.fastq
do
	SAMPLE2=$a
	SAMPLENAME2=${a%.*}	
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Start counting complement uridylated reads from file $SAMPLE2"
	echo "$SAMPLE2 number of complementary uridyl reads:"
	cat "Uridyl_compl_"$SAMPLE2 | echo $((`wc -l`/4)) >> Complement_uridyl_numbers.txt
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Done counting complement uridylated reads from file $SAMPLE2"
done
#
wait
#
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"



