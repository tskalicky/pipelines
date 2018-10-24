#!/bin/bash
#PBS -l walltime=24:0:0 
#PBS -q default@wagap-pro.cerit-sc.cz 
#PBS -l select=1:ncpus=2:mem=50gb:scratch_local=50gb:os=debian9
#PBS -j oe
#PBS -N 07_peak_calling_METTL16
export PBS_SERVER=wagap-pro.cerit-sc.cz
#
## initialize the required application
# module add gsl-1.16-intel #required by Piranha
# module add bamtools #required by Piranha
module add piranha-1.2.1
############################################################################################
### Variables
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/METTL16/peak_calling/PAR/pooled"
PARCRAC="/storage/brno3-cerit/home/tskalicky/METTL16/peak_calling/PAR/pooled/METTL16_PAR1-2.pool.tag.uniq.rgb.bed"
CONTROL="/storage/brno3-cerit/home/tskalicky/METTL16/peak_calling/PAR/pooled/FlagPAR.out.md.tag.bed"
THREADS=$PBS_NUM_PPN
APPENDIX=".bed"
#APPENDIX2=".fq"
#APPENDIX3=".sam"
####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
cp -av $PARCRAC $CONTROL $SCRATCHDIR
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR
echo "Following files/folders were copied to scratch:"
ls -c1

# Binaries
PIRANHA=$(which Piranha)

# Check if tools are installed
which $PIRANHA
####################################################################################################
### Filenames prep
PARCRAC=$(basename $PARCRAC)
CONTROL=$(basename $CONTROL)

PARCRAC_NAME=${PARCRAC%.*.*.*.*.*}
CONTROL_NAME=${CONTROL%.*.*.*.*}
### Commands
echo "Peak Calling from pooled data $PARCRAC"
$PIRANHA -s -z 20 -p 0.05 $PARCRAC $CONTROL -o $PARCRAC_NAME"_pool_peaks"
echo "Finnished Calling from pooled data $PARCRAC"

############################################################################################
### Copy data from scratch back to home dir and clean scratch
mkdir -p $OUTPUT_DIR
cp -avr $SCRATCHDIR $OUTPUT_DIR || export CLEAN_SCRATCH=false
echo "Script finished on:"
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"

