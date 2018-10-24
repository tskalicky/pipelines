#!/bin/bash
#PBS -l walltime=4:0:0 
#PBS -q default@wagap-pro.cerit-sc.cz
#PBS -l select=1:ncpus=5:mem=5gb:scratch_local=10gb
#PBS -j oe
#PBS -N 01v02_packing_libs
THREADS=$PBS_NUM_PPN # Number of threads to use
INPUT="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/STAR_index"
INPUT_index="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/STAR_index/GRCh38_with_RD_LENGTH_100_STAR2.6c"
#
cd $INPUT
# pigz -v -p $THREADS $INPUT_D3L2_SpikeIn/*.fastq
tar -I pigz -p $THREADS -cvf GRCh38_with_RD_LENGTH_100_STAR2.6c.tar.gz $INPUT_index
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"