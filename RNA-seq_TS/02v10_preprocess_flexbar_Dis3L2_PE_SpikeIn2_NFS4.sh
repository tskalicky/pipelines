#!/bin/bash
#PBS -l walltime=96:0:0 
#PBS -q uv@wagap-pro.cerit-sc.cz
#PBS -l select=1:ncpus=8:mem=50gb:scratch_local=100gb
#PBS -j oe
#PBS -N 02v10_flexbar_fastqc_Dis3L2_SpikeIn
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
INPUT_D3L2_SpikeIn="/storage/brno3-cerit/home/tskalicky/Dis3L2/Dasa_spikein/new_upscale/data/RAW"
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/Dis3L2/Dasa_spikein/new_upscale/data"
#
ADAPTERS="/storage/brno3-cerit/home/tskalicky/tools/adapters_v4.fa"
APPENDIX=".fastq.gz"
THREADS=$PBS_NUM_PPN # Number of threads to use
#
FASTQC=$(which fastqc)
FLEXBAR=$(which flexbar)
TRIMMOMATIC="java -jar /software/trimmomatic/0.36/dist/jar/trimmomatic-0.36.jar"
# Check the tools versions
echo $FASTQC
echo $FLEXBAR
echo $TRIMMOMATIC
#
####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
cd $INPUT_D3L2_SpikeIn
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
cp -av $INPUT_D3L2_SpikeIn/*.fastq.gz $ADAPTERS $SCRATCHDIR
#
cd $SCRATCHDIR
#
echo "Following files were copied to scratch:"
ls -c1

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
echo "Now I am processing Dis3L2 WT and KO PE reads - Flexbar ANY and quality Trimming"
# # Zvoleno trimovani na kvalitu az po trimovani adapteru!
flexbar -t $TRIMDIR"DIS3l2_upscale_U12_nucl_flexbar" -r 65091_GTCCGC_CCCU1ANXX_7_20180615B_20180615_R1.fastq.gz -p 65091_GTCCGC_CCCU1ANXX_7_20180615B_20180615_R2.fastq.gz -q TAIL -qf i1.8 -qt 30 --qtrim-post-removal -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 4 -z GZ &
flexbar -t $TRIMDIR"DIS3l2_upscale_U12_nonfrac_flexbar" -r 65090_CCGTCC_CCCU1ANXX_7_20180615B_20180615_R1.fastq.gz -p 65090_CCGTCC_CCCU1ANXX_7_20180615B_20180615_R2.fastq.gz -q TAIL -qf i1.8 -qt 30 --qtrim-post-removal -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 4 -z GZ &
flexbar -t $TRIMDIR"DIS3l2_upscale_U12_cyto_flexbar" -r 65092_GTGAAA_CCCU1ANXX_7_20180615B_20180615_R1.fastq.gz -p 65092_GTGAAA_CCCU1ANXX_7_20180615B_20180615_R2.fastq.gz -q TAIL -qf i1.8 -qt 30 --qtrim-post-removal -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 4 -z GZ &
flexbar -t $TRIMDIR"DIS3l2_upscale_OAT_nucl_flexbar" -r 65085_AGTTCC_CCCU1ANXX_4_20180615B_20180615_R1.fastq.gz -p 65085_AGTTCC_CCCU1ANXX_4_20180615B_20180615_R2.fastq.gz -q TAIL -qf i1.8 -qt 30 --qtrim-post-removal -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 4 -z GZ &
flexbar -t $TRIMDIR"DIS3l2_upscale_OAT_nonfrac_flexbar" -r 65084_AGTCAA_CCCU1ANXX_4_20180615B_20180615_R1.fastq.gz -p 65084_AGTCAA_CCCU1ANXX_4_20180615B_20180615_R2.fastq.gz -q TAIL -qf i1.8 -qt 30 --qtrim-post-removal -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 4 -z GZ &
flexbar -t $TRIMDIR"DIS3l2_upscale_OAT_cyto_flexbar" -r 65086_ATGTCA_CCCU1ANXX_4_20180615B_20180615_R1.fastq.gz -p 65086_ATGTCA_CCCU1ANXX_4_20180615B_20180615_R2.fastq.gz -q TAIL -qf i1.8 -qt 30 --qtrim-post-removal -a $SCRATCHDIR"/adapters_v4.fa" -ac -ae ANY -ao 3 -u 2 -m 17 -n 4 -z GZ
#
#Wait for all jobs to finish before continuing the script
wait
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Done processing Dis3L2 WT and KO PE reads - Flexbar ANY and quality Trimming"

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
# rm $SCRATCHDIR/*$APPENDIX
cp -avr $SCRATCHDIR $OUTPUT_DIR || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
