#!/bin/bash
#PBS -l select=1:ncpus=12:mem=100gb:scratch_local=400gb
#PBS -l walltime=48:00:00
#PBS -j oe
#PBS -N FastQC_trimming_FTO_data
# initialize the required application
module add fastQC-0.11.5
module add trimmomatic-0.36
module add jdk-8
# copy data to fast local scratch
# secure copy using scp (disks not connected via NFS)
DATADIR="storage-brno9-ceitec.metacentrum.cz:~/FTO_project"
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
trap 'clean_scratch' TERM EXIT
scp -c blowfish -r $DATADIR/raw_reads $DATADIR/adapters.fasta $SCRATCHDIR
# use "scp -R ..." in case of copying directories
cd $SCRATCHDIR
#define variables
trimc_dir='/software/trimmomatic/0.36/dist/jar/'
trim='java -jar '$trimc_dir'trimmomatic-0.36.jar'
adapters='$SCRATCHDIR/adapters.fasta'
threads=12
mkdir trimmed_reads trimming_logs
file_dir='$SCRATCHDIR'
trim_dir='$SCRATCHDIR/trimmed_reads'
log_dir='$SCRATCHDIR/trimming_logs'

# commands
cd $file_dir
#quality assesment
fastqc --threads 12 SRR3290146-TREX_WT_1_pass_1.fastq.gz SRR3290146-TREX_WT_1_pass_2.fastq.gz \
SRR3290147-TREX_WT_2_pass_1.fastq.gz SRR3290147-TREX_WT_2_pass_2.fastq.gz \
SRR3290148-TREX_WT_3_pass_1.fastq.gz SRR3290148-TREX_WT_3_pass_2.fastq.gz \
SRR3290143-FTO_KO_1_pass_1_fastq.gz SRR3290143-FTO_KO_1_pass_2_fastq.gz \
SRR3290144-FTO_KO_2_pass_1_fastq.gz SRR3290144-FTO_KO_2_pass_2_fastq.gz \
SRR3290145-FTO_KO_3_pass_1.fastq.gz SRR3290145-FTO_KO_3_pass_2.fastq.gz \
#
#fastqc --threads 4 SRR3290136-FTO_CLIP1_pass_1_fastq.gz SRR3290136-FTO_CLIP1_pass_2_fastq.gz \
#SRR3290137-FTO_CLIP2_pass_1.fastq.gz SRR3290137-FTO_CLIP2_pass_2.fastq.gz \
#SRR3290138-FTO_CLIP3_pass_1.fastq.gz SRR3290138-FTO_CLIP3_pass_2.fastq.gz \
#SRR3290139-FTO_Input1_pass_1.fastq.gz SRR3290139-FTO_Input1_pass_2.fastq.gz \
#SRR3290140-FTO_Input2_pass_1.fastq.gz SRR3290140-FTO_Input2_pass_2.fastq.gz \
#SRR3290142-FTO_Input3_pass_1.fastq.gz SRR3290142-FTO_Input3_pass_2.fastq.gz
#
#trimming specific for iCLIP data with minlen:20 settings
for f in *_1.fastq.gz:
    name = $f
    p_out_fw=$trim_dir$name'_trimmed_paired_out_fw.fastq.gz'
    u_out_fw=$trim_dir$name'_trimmed_unpaired_out_fw.fastq.gz'
    for i in *_2.fastq.gz:
    	name2 = $i
	    p_out_rv=$trim_dir$name2'_trimmed_paired_out_rv.fastq.gz'
    	u_out_rv=$trim_dir$name2'_trimmed_unpaired_out_rv.fastq.gz'

    	log=$log_dir$name'_trimming.log'

    	file_fw=$file_dir$name'_1.fastq.gz'
    	file_rv=$file_dir$name'_2.fastq.gz'
    	illuminaclip='ILLUMINACLIP:'$adapters':2:30:10'

    	$trim PE -threads $threads -trimlog $log $file_fw $file_rv $p_out_fw $u_out_fw $p_out_rv $u_out_rv $illuminaclip LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20


#java -jar /software/trimmomatic/0.36/dist/jar/trimmomatic-0.36.jar PE -threads 4 -trimlog trmming.log 
