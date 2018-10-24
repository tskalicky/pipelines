#!/bin/bash
#PBS -l nodes=1:ppn=20
#PBS -l mem=30gb
#PBS -k oe
#PBS -N 01v02_packing_libs
#PBS -M tomas.skalicky@seznam.cz
#PBS -m abe
#
#PARU KRTECEK server is using TORQUE scheduling system !!!
THREADS=$PBS_NUM_PPN # Number of threads to use
INPUT="/home/users/tskalicky/CEITEC/genomes/human/ensembl91/STAR_index"
INPUT_index="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/STAR_index/GRCh38_with_RD_LENGTH_100_STAR2.6c"
#
cd $INPUT
# pigz -v -p $THREADS $INPUT_D3L2_SpikeIn/*.fastq
# tar -I pigz -p $THREADS -cvf GRCh38_with_RD_LENGTH_100_STAR2.6c.tar.gz $INPUT_index
tar -I unpigz -p $THREADS -xvf $INPUT/"GRCh38_with_RD_LENGTH_100_STAR2.6c.tar.gz"
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"