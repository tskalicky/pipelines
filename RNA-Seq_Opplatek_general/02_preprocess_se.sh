#!/bin/bash
#PBS -l select=1:ncpus=6:mem=25gb:scratch_local=150gb
#PBS -l walltime=12:00:00
#PBS -q default
#PBS -N 02_preprocess_se
#
# Preprocessing - RNA-Seq SE
#
# Requires FastQC, Trimmomatic (java), multiQC, seqtk
#
#qsub -l walltime=24:0:0 -q default -l select=1:ncpus=6:mem=15gb:scratch_local=150gb
#
# Added seqtk trimming - testing
#
# TODO - check bbduk and seqtk trimming https://www.biostars.org/p/97848/
# TODO - add filtering based on N content - prinseq could do that
#
############################################################################################
### Variables
INPUT_DIR=/mnt/nfs/home/323639/999993-Bioda/projects/honza/blazek_rnaseq/data/raw
OUTPUT_DIR=/mnt/nfs/home/323639/999993-Bioda/projects/honza/blazek_rnaseq
OUTPUT_DIR_QC=${OUTPUT_DIR}/results/qc
OUTPUT_DIR=${OUTPUT_DIR}/data/preprocessed

TRIM_AD="false" # Trim also adapters? true or false
ADAPTERS="/storage/brno2/home/opplatek/tools/Trimmomatic-0.36/adapters/TruSeq3-SE.fa" # Adapters to trim, required if $TRIM_AD is true. Fasta format where /1 or /2 at the end \
#of sequence name signalizes which read is scanned for the adapter. If nothing is specified both reads are checked. You can also use _rc to scan reverse-complement. Please see Trimmomatic manual http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
RD_LENGTH=250 # Read length from the sequencing OR trimming down to this length. Illumina sometimes reports N+1 http://seqanswers.com/forums/archive/index.php/t-31154.html
TRIM_LE="false" # Force trim 5' end? true (QuantSeq FWD) 12 bases or false (QuantSeq REV), for QuantSeq SENSE it should be true and 10 bases (and possibly 6 bases from the right)
TRIM_RI="false" # Force trim of 3' end?

TRIM_LEFT=10 # Applied only if trim left is true
TRIM_RIGHT=0 # Applied only if trim right is true

APPENDIX=".fastq.gz"

PHRED_TRIM=5 # Trim the 3' end of read if four consequent bases have average base quality smaller than this value
LEN_FILTER=35 # Filter sequences shorter than this value after quality trimming

THREADS=$PBS_NUM_PPN # Number of threads to use

module add fastQC-0.11.5
FASTQC=$(which fastqc)
module add python27-modules-gcc
#PYTHONPATH=$PYTHONPATH:/storage/brno2/home/opplatek/tools/MultiQC-1.0/lib/python2.7/site-packages/
#MULTIQC=/storage/brno2/home/opplatek/tools/MultiQC-1.0/bin/multiqc
MULTIQC=$(which multiqc)
SEQTK="/storage/brno2/home/opplatek/tools/seqtk-1.2/seqtk"
module add jdk-8
TRIMMOMATIC="java -jar /storage/brno2/home/opplatek/tools/Trimmomatic-0.36/trimmomatic-0.36.jar"

# Check the tools versions
echo $FASTQC
echo $MULTIQC
echo $TRIMMOMATIC
echo $SEQTK

############################################################################################
### Copy inputs
cp $INPUT_DIR/*$APPENDIX $SCRATCH/

cd $SCRATCH/

############################################################################################
### Preparation and extra flags for Trimmomatic
msg=""
extra_flags=""
extra_flags_seqtk=""

if [ "$TRIM_AD" == "true" ]; then
	msg="Adapter and quality trimming"
	extra_flags="ILLUMINACLIP:${ADAPTERS}:2:30:10:3"
	cp $ADAPTERS $SCRATCH/
	ADAPTERS=$(basename $ADAPTERS)
else
	msg="Quality trimming"
fi

#if [ "$TRIM_LE" == "true" ]; then
#	msg="$msg, $TRIM_LEFT left"
#	extra_flags="$extra_flags HEADCROP:$TRIM_LEFT"
#else
#	msg="$msg, nothing on left"
#fi
if [ "$TRIM_LE" == "true" ]; then
	msg="$msg, $TRIM_LEFT left"
	extra_flags_seqtk="$extra_flags_seqtk -b $TRIM_LEFT"
else
	msg="$msg, nothing on left of"
fi
if [ "$TRIM_RI" == "true" ]; then
	msg="$msg, $TRIM_RIGHT right"
	extra_flags_seqtk="$extra_flags_seqtk -e $TRIM_RIGHT"
else
	msg="$msg, nothing on right"
fi

echo "Running $msg."

############################################################################################
### Preprocess and qc
mkdir -p $SCRATCH/fastqc/preprocessed

# Force trim on left/right if requested
if [ "$TRIM_LE" == "true" ] || [ "$TRIM_RI" == "true" ]; then
	for i in *$APPENDIX
	do
		$SEQTK trimfq $extra_flags_seqtk $i | gzip > ${i}.trim
		mv ${i}.trim $i
	done 
else
	echo "Not running any forced trimming"
fi	

mkdir $SCRATCH/trimmomatic

for i in *$APPENDIX
do
	READ_FOR=$i

	echo "Now I am processing SE reads $READ_FOR - Trimming"

	$TRIMMOMATIC SE -threads $THREADS -phred33 $READ_FOR $SCRATCH/${READ_FOR%$APPENDIX}_trim.fastq.gz \
	CROP:$RD_LENGTH LEADING:3 TRAILING:3 SLIDINGWINDOW:4:$PHRED_TRIM MINLEN:$LEN_FILTER $extra_flags &> $SCRATCH/trimmomatic/${READ_FOR%$APPENDIX}_trim.log # -trimlog $SCRATCH/trimmomatic/${READ_FOR%$APPENDIX}_trim.log

	echo "Done processing SE reads $READ_FOR - Trimming"

done

$FASTQC --outdir $SCRATCH/fastqc/preprocessed --format fastq --threads $THREADS $SCRATCH/*_trim.fastq.gz

cd $SCRATCH/fastqc/preprocessed

echo "Number of reads after preprocessing is in fastqc/preprocessed/total_sequences.txt"
for a in *zip
do
	unzip -q $a
	grep "Total Sequences" ${a%.zip*}/*.txt
done > total_sequences.txt
rm -r *fastqc

$MULTIQC -o $SCRATCH/fastqc/preprocessed $SCRATCH/fastqc/preprocessed/
$MULTIQC -o $SCRATCH/trimmomatic $SCRATCH/trimmomatic/

############################################################################################
### Copy results
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR_QC/fastqc

echo "Output will be stored in $OUTPUT_DIR and $OUTPUT_DIR_QC"

cp -r $SCRATCH/fastqc/preprocessed $OUTPUT_DIR_QC/fastqc/ &
cp -r $SCRATCH/trimmomatic $OUTPUT_DIR_QC/ &
cp $SCRATCH/*_trim.fastq.gz $OUTPUT_DIR/

rm -r $SCRATCH/*
