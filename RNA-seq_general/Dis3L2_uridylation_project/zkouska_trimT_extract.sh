#!/bin/bash
for b in *.fastq
do
	SAMPLE2=$b
	SAMPLENAME2=${b%.*}
	if [[ $SAMPLE2 = DIS3l2_*_unmapped.fastq ]] ; then
		# date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
		echo "Found UNMAPPED reads for trimming in file $SAMPLE2"
		echo "Start trimming uridylated reads from unmapped reads in file $SAMPLE2"
		# python trimT_final.py $SAMPLE2 &
		# echo "Start extracting uridylated reads from FASTQ file $SAMPLE2"
		# python extract_uridylated.py $SAMPLE2 &
		# echo "Start extracting uridylated reads from FASTA and fasta file $SAMPLENAME2.fa"
		# python extract_uridylated_out_fasta.py $SAMPLENAME2".fa" &
	elif [[ $SAMPLE2 = DIS3l2_*_uniq_mapped.fastq ]]; then
		echo "Found MAPPED reads for extracting in file $SAMPLE2"
		# echo "Start extracting uridylated reads from FATSQ file $SAMPLE2"
		# python extract_uridylated.py $SAMPLE2 &
		# echo "Start extracting uridylated reads from FASTA and fasta file $SAMPLENAME2.fa"
		# python extract_uridylated_out_fasta.py $SAMPLENAME2".fa" &
	fi
done