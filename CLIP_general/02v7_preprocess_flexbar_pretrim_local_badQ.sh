#!/bin/bash
### Variables
INPUT_FlagPAR_UV="/Users/tomas/Data/owncloud/CEITEC_lab/METTL16/data/trimmed/FlagPAR_FlagUV_jobs/FlagPAR_FlagUV_3_job_563293_wagap"
INPUT_METTL16_PAR1="/Users/tomas/Data/owncloud/CEITEC_lab/METTL16/data/trimmed/METTL16_PAR1_job_4260487_arien/METTL16PAR-1_SRR6048658_flexbar_qual.fastq.gz"
INPUT_METTL16_UV1="/Users/tomas/Data/owncloud/CEITEC_lab/METTL16/data/trimmed/METTL16_UV1_job_4260489_arien/METTL16UV-1_SRR6048655_flexbar_qual.fastq.gz"
OUTPUT_DIR="/Users/tomas/Data/owncloud/CEITEC_lab/METTL16/data/trimmed/left_cut"

APPENDIX=".fastq.gz"
THREADS=4 # Number of threads to use

############################################################################################
# commands
cd $OUTPUT_DIR
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
#mkdir -p $OUTPUT_DIR/{data/trimmed,fastqc/trimmed}  #creates whole subdirectory tree using -p and {}
QCDIR=$OUTPUT_DIR/fastqc/trimmed/
TRIMDIR=$OUTPUT_DIR/data/trimmed/
echo $QCDIR
echo $TRIMDIR

#usage of $FLEXBAR and $FASTQC variables in function is not working in PBS "Illegal instruction" error
#loop over every file and process them in background simultaniously
#FLEXBAR is using only like 1/2 the cores assigned
echo "Now I am processing METTL16 SE reads - Flexbar pre-trim-left Trimming"
#flexbar -t $TRIMDIR"FlagPAR_SRR6048657_left_cut_flexbar" -r $INPUT_FlagPAR_UV/FlagPAR_SRR6048657_flexbar_qual3.fastq.gz --pre-trim-left 9 --min-read-length 17 -n 4 -z GZ
#flexbar -t $TRIMDIR"FlagUV_SRR6048654_left_cut_flexbar" -r $INPUT_FlagPAR_UV/FlagUV_SRR6048654_flexbar_qual3.fastq.gz --pre-trim-left 10 --min-read-length 17 -n 4 -z GZ
flexbar -t $TRIMDIR"METTL16PAR-1_SRR6048658_left_cut_flexbar" -r $INPUT_METTL16_PAR1 --pre-trim-left 9 --min-read-length 17 -n 4 -z GZ
#flexbar -t $TRIMDIR"METTL16UV-1_SRR6048655_left_cut_flexbar" -r $INPUT_METTL16_UV1 --pre-trim-left 10 --min-read-length 17 -n 4 -z GZ
# Wait for all jobs to finish before continuing the script
wait
echo "Done processing METTL16 SE reads - Flexbar pre-trim-left Trimming"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"

cd $TRIMDIR

#bylo nutne upravit volani souboru pro fastqc
/Users/tomas/bin/FastQC.app/Contents/MacOS/fastqc --outdir $QCDIR --format fastq --threads 4 *.gz
#
cd $QCDIR
echo "Number of reads after preprocessing is in fastqc/trimmed/total_sequences.txt"
for i in *zip
do
	unzip -q $i
	grep -hF "Total Sequences" ${i%.zip*}/*.txt
done > total_sequences.txt
rm -r ./*fastqc/
#
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
