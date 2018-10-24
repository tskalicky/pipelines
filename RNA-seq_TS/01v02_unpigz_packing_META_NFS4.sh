#!/bin/bash
#PBS -l walltime=4:0:0 
#PBS -q default@wagap-pro.cerit-sc.cz
#PBS -l select=1:ncpus=5:mem=4gb:scratch_local=10gb
#PBS -j oe
#PBS -N 01v02_packing_fastqc_Dis3L2_SpikeIn
#
# Packing RAW reads
#
# Requires pigz
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
APPENDIX=".fastq"
THREADS=$PBS_NUM_PPN # Number of threads to use
#
PIGZ=$(which pigz)
# Check the tools versions
echo $PIGZ
#
####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
cd /storage/brno3-cerit/home/tskalicky/Dis3L2/Dasa_spikein/RAW
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
#
# Copy to scratch will not be performed this time !
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR

############################################################################################
# commands
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
#
cd $INPUT_D3L2_SpikeIn
pigz -v -p $THREADS $INPUT_D3L2_SpikeIn/*.fastq
# tar -I pigz -p $THREADS -cvf GRCh38_with_RD_LENGTH_100_STAR2.6c.tar.gz $INPUT_index
# tar -I unpigz -p $THREADS -xvf $INPUT/"GRCh38_with_RD_LENGTH_100_STAR2.6c.tar.gz"
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"