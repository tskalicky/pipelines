#!/bin/bash
#PBS -l walltime=96:0:0 
#PBS -q uv@wagap-pro.cerit-sc.cz
#PBS -l select=1:ncpus=20:mem=150gb:scratch_local=450gb
#PBS -j oe
#PBS -tN Dis3L2_spikein_uniq_mapped_extract
#
# For 1 library you need cca 35GB RAM
#
## initialize the required application
# module add samtools-1.4 # Samtools v1.6 are broken! 
# module add python27-modules-gcc #required by multiQC
############################################################################################
### Variables
THREADS="30"
INPUT_DIR="/home/tomas/Bioinformatics/CEITEC/Dis3L2/Dasa_spikein/mapping"
OUTPUT_DIR="/home/tomas/Bioinformatics/CEITEC/Dis3L2/Dasa_spikein/uridylation"
#
APPENDIX1=".bam"
APPENDIX2=".fastq"
APPENDIX3=".fa"
# Binaries
SAMTOOLS=$(which samtools)
REPAIR=$(which repair.sh)
MULTIQC=$(which multiqc)
# Check the tools versions
which $SAMTOOLS
which $REPAIR
which $MULTIQC

####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
# mkdir -p $OUTPUT_DIR
# cd $INPUT_DIR
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
# find "$INPUT_DIR" -maxdepth 1 -tName "*.bam" -exec cp -vt "$OUTPUT_DIR" {} +
# find "$INPUT_DIR" -maxdepth 2 -tName "*.bam" -exec cp -vt "$OUTPUT_DIR" {} + # this is not working at every occasion!!
cd $OUTPUT_DIR

if [ ! -d "$OUTPUT_DIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "OUTPUT_DIR path is:" $OUTPUT_DIR
echo "Following files were copied to scratch:"
ls -Rc1
####################################################################################################
## COMMANDS
# Warning! This samtools script is written for version 1.6 and above! 
# Others do NOT have option -tN to force include of the /1 and /2 to read names
# for a in *$APPENDIX1
# do
# 	SAMPLE=$a
# 	SAMPLENAME=${a%.*}
# 	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# 	echo "Start extracting uniquely mapped PE reads from file $SAMPLE"
# 	# to extract reads mapped only 1 time concordantly
# 	samtools view --threads $THREADS -b -hf 0x2 -o $SAMPLENAME"_uniq_mapped.bam" $SAMPLE
# 	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# 	echo "Start extracting unmapped PE reads from file $SAMPLE"
# 	samtools view --threads $THREADS -b -hf 0x4 -o $SAMPLENAME"_unmapped.bam" $SAMPLE
# 	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# 	echo "Start converting uniquely mapped BAM file $SAMPLE to fastq"
# 	samtools fastq --threads $THREADS -tN $SAMPLENAME"_uniq_mapped.bam" > $SAMPLENAME"_uniq_mapped.fastq" 
# 	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# 	echo "Start converting unmapped BAM file $SAMPLE to fastq"
# 	samtools fastq --threads $THREADS -tN $SAMPLENAME"_unmapped.bam" > $SAMPLENAME"_unmapped.fastq"
# 	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# 	echo "Start converting uniquely mapped BAM file $SAMPLE to fasta"
# 	samtools fasta --threads $THREADS -tN $SAMPLENAME"_uniq_mapped.bam" > $SAMPLENAME"_uniq_mapped.fa"
# 	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# 	echo "Start converting unmapped BAM file $SAMPLE to fasta"
# 	samtools fasta --threads $THREADS -tN $SAMPLENAME"_unmapped.bam" > $SAMPLENAME"_unmapped.fa"
# 	# #Single_End_Layout:
# 	# samtools view --threads $THREADS -b -F 4 in.bam > mapped.bam
# 	# samtools view --threads $THREADS -b -f 4 in.bam > unmapped.bam
# 	# #Paired_End_Layout
# 	# samtools view --threads $THREADS -b -f 2 in.bam > mapped.bam
# 	# samtools view --threads $THREADS -b -F 2 in.bam > unmapped.bam
# 	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# 	echo "Done all processing reads from file $SAMPLE"
# done
# # wait
# wait
####################################################################################################
#Repair broken pairs in fastq files
for d in *.fastq
do
	SAMPLE4=$d
	SAMPLENAME4=${d%.*}
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Now I am repairing pairs in PE reads $SAMPLE4"
	$REPAIR in=$SAMPLE4 out="$SAMPLENAME4"_fixed_pairs.fastq outs="$SAMPLENAME4"_all_singletons.fastq repair
	echo "Done repairing pairs in PE reads $SAMPLE4"
done
#wait
wait
####################################################################################################
### Compress output files
# date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# echo "Compressing results"
# for b in *.fastq
# do
# 	SAMPLE2=$b
# 	SAMPLENAME2=${b%.*}
# 	echo "Now I am compressing PE reads $SAMPLENAME2"
# 	pigz -v -p $THREADS $SAMPLE2
# 	echo "Done compressing PE reads $SAMPLENAME2"
# done
# #wait
# wait
# for c in *.fa
# do
# 	SAMPLE3=$c
# 	SAMPLENAME3=${c%.*}
# 	echo "Now compressing PE reads $SAMPLENAME3"
# 	pigz -v -p $THREADS $SAMPLE3
# 	echo "Done compressing PE reads $SAMPLENAME3"
# done
# #wait
# wait
# #
####################################################################################################
### Finalize and copy results
# date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# echo "Finalize and copy results"
# mkdir -p $OUTPUT_DIR
# #rm -v $SCRATCHDIR/*.bam
# mv -v $SCRATCHDIR/*.fastq.gz $OUTPUT_DIR
# mv -v $SCRATCHDIR/*.fa.gz $OUTPUT_DIR
# mv -v $SCRATCHDIR/* $OUTPUT_DIR
# rm -rv $SCRATCHDIR/*
# echo "Script finished on:"
# date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
# # samtools view --threads 3 -b -hf 0x2 DIS3l2_OAT_cyto.GRCh38.91Aligned.sortedByCoord.out_uniq_mapped.bam > DIS3l2_OAT_cyto.GRCh38.91Aligned.sortedByCoord.out_uniq_mapped.sam
# # samtools view --threads 3 -hF 4 DIS3l2_OAT_cyto.GRCh38.91Aligned.sortedByCoord.out_uniq_mapped.bam > DIS3l2_OAT_cyto.GRCh38.91Aligned.sortedByCoord.out_unique_uniq_mapped.sam