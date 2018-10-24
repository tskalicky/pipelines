#!/bin/bash
#PBS -l walltime=96:0:0 
#PBS -q uv@wagap-pro.cerit-sc.cz
#PBS -l select=1:ncpus=20:mem=150gb:scratch_local=450gb
#PBS -j oe
#PBS -N Dis3L2_spikein_mapped_extract
#
############################################################################################
# Uridylated RNAseq reads processing
# Requires samtools and bbmap
############################################################################################
### Variables
THREADS=30
INPUT_DIR="/home/tomas/Bioinformatics/CEITEC/Dis3L2/Dasa_spikein/mapping/upscale"
OUTPUT_DIR="/home/tomas/Bioinformatics/CEITEC/Dis3L2/Dasa_spikein/uridylation/UPSCALE"
EXTRACT_FASTA="/home/tomas/Bioinformatics/CEITEC/Dis3L2/Dasa_spikein/uridylation/extract_uridylated_out_fasta.py"
EXTRACT="/home/tomas/Bioinformatics/CEITEC/Dis3L2/Dasa_spikein/uridylation/extract_uridylated.py"
TRIMT="/home/tomas/Bioinformatics/CEITEC/Dis3L2/Dasa_spikein/uridylation/trimT_final.py"
export PATH="$PATH:/home/tomas/miniconda2/bin/"
echo "PATH after modification is:"
echo $PATH | tr ":" "\n" | nl
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
echo "Required tools are installed here:"
which $SAMTOOLS
which $MULTIQC
which $REPAIR
which $BBSPLITPAIRS

####################################################################################################
## copy input data using OUTPUT_DIR storage which is shared via NFSv4
## clean the SCRATCH when job finishes (and data
## are successfully copied out) or is killed
## use cp -avr when copying directories
cd $INPUT_DIR
mkdir -p $OUTPUT_DIR
find "$INPUT_DIR" -maxdepth 1 -name "*.bam" -exec cp -vt "$OUTPUT_DIR" {} +
cp -av $EXTRACT $EXTRACT_FASTA $TRIMT $OUTPUT_DIR
cd $OUTPUT_DIR
#
if [ ! -d "$OUTPUT_DIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "OUTPUT_DIR path is:" $OUTPUT_DIR
echo "Following files were copied to scratch:"
ls -Rc1
####################################################################################################
## COMMANDS
set +f # enable globbing, because is needed for other parts of the script, like variable expansions !!!
# Samtools for extraction of unmapped and uniquelly mapped reads from original bam files created by STAR mapper.
# Will create both fastq and fasta files
# Warning! This samtools script is written for version 1.6 and above! 
# Others do NOT have option -tN to force include of the /1 and /2 to read names
## USAGE:
# #Single_End_Layout:
# samtools view --threads $THREADS -b -F 4 in.bam > mapped.bam
# samtools view --threads $THREADS -b -f 4 in.bam > unmapped.bam
# #Paired_End_Layout
# samtools view --threads $THREADS -b -f 2 in.bam > mapped.bam
# samtools view --threads $THREADS -b -F 2 in.bam > unmapped.bam
for a in *$APPENDIX1
do
	SAMPLE=$a
	SAMPLENAME=${a%.*}
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Start extracting uniquely mapped PE reads from file $SAMPLE"
	# to extract reads mapped only 1 time concordantly
	samtools view --threads $THREADS -b -hf 0x2 -o $SAMPLENAME"_uniq_mapped.bam" $SAMPLE
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Start extracting unmapped PE reads from file $SAMPLE"
	samtools view --threads $THREADS -b -hf 0x4 -o $SAMPLENAME"_unmapped.bam" $SAMPLE
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Start converting uniquely mapped BAM file $SAMPLE to fastq"
	samtools fastq --threads $THREADS -tN $SAMPLENAME"_uniq_mapped.bam" > $SAMPLENAME"_uniq_mapped.fastq" 
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Start converting unmapped BAM file $SAMPLE to fastq"
	samtools fastq --threads $THREADS -tN $SAMPLENAME"_unmapped.bam" > $SAMPLENAME"_unmapped.fastq"
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Start converting uniquely mapped BAM file $SAMPLE to fasta"
	samtools fasta --threads $THREADS -tN $SAMPLENAME"_uniq_mapped.bam" > $SAMPLENAME"_uniq_mapped.fa"
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Start converting unmapped BAM file $SAMPLE to fasta"
	samtools fasta --threads $THREADS -tN $SAMPLENAME"_unmapped.bam" > $SAMPLENAME"_unmapped.fa"
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Done all processing reads from file $SAMPLE"
done
# wait
wait
#
### Trimm uridylation (TTT) from UNMAPPED fastq reads and extract uridylated reads from fastq and fasta files
#
for b in *.fastq
do
	SAMPLE2=$b
	SAMPLENAME2=${b%.*}
	if [[ $SAMPLE2 = DIS3l2_*_unmapped.fastq ]] ; then
		# date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
		echo "Found UNMAPPED reads for trimming in file $SAMPLE2"
		date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
		echo "Start trimming uridylated reads from unmapped reads in file $SAMPLE2"
		python trimT_final.py $SAMPLE2 &
		echo "Start extracting uridylated reads from FASTQ file $SAMPLE2"
		python extract_uridylated.py $SAMPLE2 &
		echo "Start extracting uridylated reads from FASTA and fasta file $SAMPLENAME2.fa"
		python extract_uridylated_out_fasta.py $SAMPLENAME2".fa" &
	elif [[ $SAMPLE2 = DIS3l2_*_uniq_mapped.fastq ]]; then
		echo "Found MAPPED reads for extracting in file $SAMPLE2"
		date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
		echo "Start extracting uridylated reads from FATSQ file $SAMPLE2"
		python extract_uridylated.py $SAMPLE2 &
		echo "Start extracting uridylated reads from FASTA and fasta file $SAMPLENAME2.fa"
		python extract_uridylated_out_fasta.py $SAMPLENAME2".fa" &
	fi
done
#wait
wait
#
for c in trimT_*.fastq
do
	SAMPLE3=$c
	SAMPLENAME3=${c%.*}
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Now I am repairing pairs and discarding too short reads in $SAMPLE3"
	$BBSPLITPAIRS  in=$SAMPLE3 out="$SAMPLENAME3"_good_pairs.fastq outsingle="$SAMPLENAME3"_singletons.fastq minlen=17
	$REPAIR in="$SAMPLENAME3"_good_pairs.fastq out="$SAMPLENAME3"_fixed_good_pairs_R1.fastq out2="$SAMPLENAME3"_fixed_good_pairs_R2.fastq outs="$SAMPLENAME3"_singletons2.fastq repair
	echo "Done repairing pairs and discarding too short reads in $SAMPLE3"
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Now I am merging singletons from original file $SAMPLE3"
	cat "$SAMPLENAME3"_singletons.fastq "$SAMPLENAME3"_singletons2.fastq >> "$SAMPLENAME3"_all_singletons.fastq
done
#wait
wait
#
#############################################################################################################
### Compress output files
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Compressing results"
for d in *.fa
do
	SAMPLE4=$d
	echo "Now compressing PE reads $SAMPLE4"
	pigz -v -p $THREADS $SAMPLE4
	echo "Done compressing PE reads $SAMPLE4"
done
#wait
wait
for e in *.fastq
do
	SAMPLE5=$e
	echo "Now I am compressing PE reads $SAMPLE5"
	pigz -v -p $THREADS $SAMPLE5
	echo "Done compressing PE reads $SAMPLE5"
done
#wait
wait
#####################################################################################################
### Finalize and copy results
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Finalize and copy results"
mkdir -p $OUTPUT_DIR/{extract_from_orig/{bam,fastq,fasta},trimT,sorted/{fastq,fasta}}  #creates whole subdirectory tree using -p and {}
mv -v $OUTPUT_DIR/*.bam $OUTPUT_DIR/extract_from_orig/bam
mv -v $OUTPUT_DIR/trimT_*.fastq.gz $OUTPUT_DIR/trimT
mv -v $OUTPUT_DIR/Uridyl_*.fastq.gz $OUTPUT_DIR/sorted/fastq
mv -v $OUTPUT_DIR/NOT_uridyl_*.fastq.gz $OUTPUT_DIR/sorted/fastq
mv -v $OUTPUT_DIR/Uridyl_*.fa.gz $OUTPUT_DIR/sorted/fasta
mv -v $OUTPUT_DIR/NOT_uridyl_*.fa.gz $OUTPUT_DIR/sorted/fasta
mv -v $OUTPUT_DIR/*_uniq_mapped.fastq.gz $OUTPUT_DIR/extract_from_orig/fastq
mv -v $OUTPUT_DIR/*_unmapped.fastq.gz $OUTPUT_DIR/extract_from_orig/fastq
mv -v $OUTPUT_DIR/*_uniq_mapped.fa.gz $OUTPUT_DIR/extract_from_orig/fasta
mv -v $OUTPUT_DIR/*_unmapped.fa.gz $OUTPUT_DIR/extract_from_orig/fasta
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
