#!/bin/bash
#PBS -l walltime=96:0:0 
#PBS -q uv@wagap-pro.cerit-sc.cz
#PBS -l select=1:ncpus=20:mem=50gb:scratch_local=100gb
#PBS -j oe
#PBS -N 02v12_trimmomatic_fastqc_ABH8_RIPseq
#
# Preprocessing - RNA-Seq PE
#
# Requires FastQC, flexbar or Trimmomatic (java)
#
## initialize the required application
module add fastQC-0.11.5
module add flexbar-3.0.3
module add trimmomatic-0.36
### Variables
INPUT="/storage/brno3-cerit/home/tskalicky/ABH8/RIP-seq/preprocessed/semi_trimmed_job_1395695.wagap-pro.cerit-sc.cz/data/trimmed"
# INPUT="/storage/brno3-cerit/home/tskalicky/ABH8/RIP-seq/preprocessed/semi_trimmed/data/trimmed"
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/ABH8/RIP-seq/preprocessed"
export PATH="$PATH:/storage/brno3-cerit/home/tskalicky/anaconda2/bin/trimmomatic"
#
ADAPTERS="/storage/brno3-cerit/home/tskalicky/tools/adapters_v4.fa"
APPENDIX=".fastq.gz"
THREADS=$PBS_NUM_PPN # Number of threads to use
#
FASTQC=$(which fastqc)
FLEXBAR=$(which flexbar)
TRIMMOMATIC=$(which trimmomatic)
# TRIMMOMATIC="java -jar /software/trimmomatic/0.36/dist/jar/trimmomatic-0.36.jar"
echo $FASTQC
echo $FLEXBAR
echo $TRIMMOMATIC
#
####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
cd $INPUT
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
cp -av $INPUT/*.fastq.gz $ADAPTERS $SCRATCHDIR
#
cd $SCRATCHDIR
#
echo "Following files were copied to scratch:"
ls -c1 | tr " " "\n" | nl

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR

############################################################################################
# commands
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
mkdir -p $SCRATCHDIR/{data/trimmed,fastqc/trimmed}  #creates whole subdirectory tree using -p and {}
#
QCDIR=$SCRATCHDIR/fastqc/trimmed/
TRIMDIR=$SCRATCHDIR/data/trimmed/
echo "QCdir is:" $QCDIR
echo "TRIMdir is:" $TRIMDIR

#usage of $FLEXBAR and $FASTQC variables in function is not working in PBS "Illegal instruction" error
#loop over every file and process them in background simultaniously
#FLEXBAR is using only like 1/2 the cores assigned (no more like 3-4 cores)
# POZOR! Flexbar nebezi na Metacentru pri pouziti uzlu s os=debian9 !
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Started FastQC evaluation for all RAW libraries"
fastqc --outdir $SCRATCHDIR/fastqc --format fastq --threads $THREADS *.gz
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Finnished FastQC evaluation for all RAW libraries"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Now I am processing ABH8 RIPseq PE reads - Flexbar ANY and quality Trimming"
### Zvoleno trimovani na kvalitu az po trimovani adapteru! --qtrim-post-removal
# flexbar -t $TRIMDIR"ALKBH8_RIPseq_flexbar" -r ALKBH8_RIPseq_R1.fastq.gz -p ALKBH8_RIPseq_R2.fastq.gz -q TAIL -qf i1.8 -qt 30 --qtrim-post-removal -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n $THREADS -z GZ
# flexbar -t $TRIMDIR"ALKBH8_RIPseq_flexbar" -r ALKBH8_RIPseq_1.fastq.gz -p ALKBH8_RIPseq_2.fastq.gz -q TAIL -qf i1.8 -qt 30 -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n $THREADS -z GZ
### It seems that the quality based trimming is not performed. Another run with quality based trimming only. 
# flexbar -t $TRIMDIR"ALKBH8_RIPseq_flexbar" -r ALKBH8_RIPseq_1.fastq.gz -p ALKBH8_RIPseq_2.fastq.gz -q TAIL -qf i1.8 -qt 30 -m 17 -n $THREADS -z GZ
### Need to use trimmomatic which will trimm low quality reads which quality fals below selected number in certain window (kmer)
####################################################################################################
shopt nullglob # If set, Bash allows filename patterns which match no files to expand to a null string, rather than themselves.
IFS=$'\n' # split on newline only, needed for filling arrays with output from find
set -f    # disable globbing, needed for filling arrays with output from find
# declare an array variables
# prepared for mutiple libraries in future scripts
declare -a ABH8_RIPseq=($(find . -maxdepth 1 -name '*ALKBH8_RIPseq*.fastq.gz' -exec basename {} \; | sort -n))
# need to add other variable names if more libraries
declare -a ARRAY_NAMES=("ABH8_RIPseq")
# need to add other variables if more libraries
declare -a ALL_ARRAYS=("${ABH8_RIPseq[@]}")
set +f # enable globbing, because is needed for other parts of the script, like variable expansions !!!
# get length of an array
ABH8_RIPseq_num="${#ABH8_RIPseq[@]}"
lib_count="${#ARRAY_NAMES[@]}"
echo "There are $ABH8_RIPseq_num ABH8_RIPseq samples that will be trimmed."
echo "Sample names are: ${ABH8_RIPseq[@]}"
# Set MAX RAM and CPUs
RAM=$[$MY_RAM-1] # lower it by one to leave some space
MEM_LIMIT=$RAM # $[$RAM/$THREADS] # In case we need GB per thread
# In case we need to divide CPUs between several running jobs
# To round up (Function ceiling):
THREADS_LIMIT="$(echo "$THREADS $lib_count" | awk '{print int( ($1/$2) + 1 )}')"
# To round down (Function floor):
# THREADS_LIMIT="$(echo "$THREADS $quantity" | awk '{print int($1/$2)}')"
#
echo "There are $quantity LIBRARIES that will be trimmed."
echo "Sample names are:"
echo ${ALL_ARRAYS[@]} | tr " " "\n" | nl
####################################################################################################
### Trimming
# --limitSjdbInsertNsj 4000000 = need to increase limit for juction detection if STAR ends up after 1st pass with error
for (( w=0;w<${ABH8_RIPseq_num};w +=2 ));
do
	SAMPLE="${ALL_ARRAYS[$w]}"
	SAMPLENAME=${SAMPLE%"_1"*}
	SAMPLE2="${ALL_ARRAYS[$w+1]}"
	SAMPLENAME2=${SAMPLE2%"_2"*}
	if [[ -f $SAMPLE ]]; then
		date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
		echo "Now I am processing PE reads $SAMPLE and $SAMPLE2 - trimming"
		$TRIMMOMATIC PE -threads $THREADS -trimlog $SCRATCHDIR/$SAMPLENAME"_trimmomatic_log.txt" \
		$SAMPLE $SAMPLE2 \
		$TRIMDIR/$SAMPLENAME"_Paired1.fq.gz" $TRIMDIR/$SAMPLENAME"_Unpaired1.fq.gz" \
		$TRIMDIR/$SAMPLENAME2"_Paired2.fq.gz" $TRIMDIR/$SAMPLENAME2"_Unpaired2.fq.gz" \
		ILLUMINACLIP:adapters_v4.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:18
	else
		date +"%d/%m/%Y %H:%M:%S"
		echo "There is no $SAMPLE file for trimming!" && exit 1
	fi
done
#Wait for all jobs to finish before continuing the script
wait
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Done processing ABH8 RIPseq PE reads - Trimmomatic and quality Trimming"
echo ${ABH8_RIPseq[@]} | tr " " "\n" | nl
#
cd $TRIMDIR
#bylo nutne upravit volani souboru pro fastqc
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Started FastQC evaluation for all trimmed libraries"
fastqc --outdir $QCDIR --format fastq --threads $THREADS *.gz
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Finnished FastQC evaluation for all trimmed libraries"
#
cd $QCDIR
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Creating library summary info for all trimmed libraries"
echo "Number of reads after preprocessing is in fastqc/trimmed/total_sequences.txt"
for i in *zip
do
	unzip -q $i
	grep -hF "Total Sequences" ${i%.zip*}/*.txt
done > total_sequences.txt
rm -r ./*fastqc/
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Creating library summary info for all trimmed libraries"
echo "Number of reads in libraries are:"
cat total_sequences.txt

############################################################################################
### Copy data from scratch back to home dir and clean scratch
cd $SCRATCHDIR
rm $SCRATCHDIR/*$APPENDIX
rm $SCRATCHDIR/adapters_v4.fa
cp -avr $SCRATCHDIR $OUTPUT_DIR || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
