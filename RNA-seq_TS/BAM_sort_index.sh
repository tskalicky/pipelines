#!/bin/bash
#PBS -l walltime=4:0:0 
#PBS -q uv@wagap-pro.cerit-sc.cz
#PBS -l select=1:ncpus=10:mem=150gb:scratch_local=250gb
#PBS -j oe
#PBS -N alignment_PE_Dis3L2_KO1
#
## initialize the required application
module add samtools-1.4 # Samtools v1.6 are broken! 
module add python27-modules-gcc #required by multiQC
#
# Sorting and index creation for BAM files
#
# Requires samtools, multiQC
#
############################################################################################
### Variables
THREADS=$PBS_NUM_PPN
WT2_input="/storage/brno3-cerit/home/tskalicky/Dis3L2/mapping/STAR_genome/WT2_results/alignment/genome"
WT3_input="/storage/brno3-cerit/home/tskalicky/Dis3L2/mapping/STAR_genome/WT3_results/alignment/genome"
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/Dis3L2/mapping/STAR_genome_tRNA/KO1"
#
cd $WT2_input
for a in *.bam
do
	SAMPLE=$a
	SAMPLENAME=${a%.*}
	$SAMTOOLS sort --threads $THREADS -T tmp.sort -o ${a%.*}.bam -
	$SAMTOOLS index --threads $THREADS ${a%.*}.bam
done
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Finnished sorting and indexing file BAM."
cd $WT3_input
for b in *.bam
do
	SAMPLE=$b
	SAMPLENAME=${b%.*}
	$SAMTOOLS sort --threads $THREADS -T tmp.sort -o ${b%.*}.bam -
	$SAMTOOLS index --threads $THREADS ${b%.*}.bam
done
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Finnished sorting and indexing file BAM."
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"