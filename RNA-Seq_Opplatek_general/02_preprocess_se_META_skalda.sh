#!/bin/bash
##PBS -l select=1:ncpus=12:mem=35gb:scratch_local=150gb
##PBS -l walltime=24:00:00
##PBS -q ceitecmu
##PBS -j oe
##PBS -N 02_preprocess_SE_METTL16
#
# Preprocessing - RNA-Seq SE
#
# Requires FastQC, Trimmomatic (java), flexbar, multiQC, seqtk
#
#qsub -l walltime=24:0:0 -q default -l select=1:ncpus=12:mem=35gb:scratch_local=150gb
#
# Added seqtk trimming - testing
#
# TODO - check bbduk and seqtk trimming https://www.biostars.org/p/97848/
# TODO - add filtering based on N content - prinseq could do that
#
############################################################################################
## initialize the required application
module add fastQC-0.11.5
module add flexbar-3.0.3
module add trimmomatic-0.36
### Variables
INPUT_DIR=/storage/brno9-ceitec/home/tskalicky/METTL16/raw
OUTPUT_DIR=/storage/brno9-ceitec/home/tskalicky/METTL16
OUTPUT_DIR_QC=${OUTPUT_DIR}/results/qc
OUTPUT_DIR=${OUTPUT_DIR}/data/preprocessed

ADAPTERS=/storage/brno9-ceitec/home/tskalicky/tools/adapters_v2.fasta
APPENDIX=".fastq.gz"
PHRED_TRIM=25 # Trim the 3' end of read if four consequent bases have average base quality smaller than this value
LEN_FILTER=17 # Filter sequences shorter than this value after quality trimming
THREADS=$PBS_NUM_PPN # Number of threads to use

FASTQC=$(which fastqc)
FLEXBAR=$(which flexbar)
TRIMMOMATIC="java -jar /storage/brno9-ceitec/home/tskalicky/tools/Trimmomatic-0.36"
# Check the tools versions
echo $FASTQC
echo $FLEXBAR
echo $TRIMMOMATIC


############################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
trap 'clean_scratch' TERM EXIT
cp -avr $INPUT_DIR/*$APPENDIX $SCRATCHDIR
cd $SCRATCHDIR
# commands
date +"%m/%d/%Y %H:%M:%S $HOSTNAME"
mkdir -p $SCRATCHDIR/fastqc/trimmed
mkdir -p $SCRATCHDIR/data/trimmed
QCDIR=$SCRATCHDIR/fastqc/trimmed/
TRIMDIR=$SCRATCHDIR/data/trimmed/
#trim reads witch flexbar
for i in *$APPENDIX
do
	READ_FOR=$i
	SAMPLENAME=${i%.*.*}
	echo "Now I am processing SE reads $READ_FOR - Flexbar Trimming"
	
	$FLEXBAR -t $TRIMDIR$SAMPLENAME"_flexbar" -r $READ_FOR -a $ADAPTERS -ac -ao 3 -u 2 -m 17 -n $THREADS -z GZ -l MOD 
	$FLEXBAR -t $TRIMDIR$SAMPLENAME"_flexbar2" -r $TRIMDIR$SAMPLENAME"_flexbar"$APPENDIX -a $ADAPTERS -ac -ae LEFT -ao 3 -u 2 -m 17 -n $THREADS -z GZ -l MOD 

	echo "Done processing SE reads $READ_FOR - Flexbar Trimming"

done

$FASTQC --outdir $SCRATCHDIR/fastqc/preprocessed --format fastq --threads $THREADS $SCRATCHDIR/*_trim.fastq.gz

cd $SCRATCHDIR/fastqc/preprocessed

echo "Number of reads after preprocessing is in fastqc/preprocessed/total_sequences.txt"
for a in *zip
do
	unzip -q $a
	grep "Total Sequences" ${a%.zip*}/*.txt
done > total_sequences.txt
rm -r *fastqc



############################################################################################
### Copy data from scratch back to home dir and clean scratch
cp -avr $SCRATCHDIR $DATADIR || export CLEAN_SCRATCH=false
cd $DATADIR
echo "Script finished on:"
date +"%m/%d/%Y %H:%M:%S $HOSTNAME"






