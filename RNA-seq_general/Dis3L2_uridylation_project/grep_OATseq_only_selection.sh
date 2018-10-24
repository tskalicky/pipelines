#!/bin/bash
#PBS -l walltime=24:0:0 
#PBS -q default@wagap-pro.cerit-sc.cz
#PBS -l select=1:ncpus=10:mem=200gb:scratch_local=200gb
#PBS -j oe
#PBS -N repair_broken_pairs_01
# For interactive jobs on Metacentrum use qsub -I -l select=...
qsub -I -l walltime=24:0:0 -q default@wagap-pro.cerit-sc.cz -l select=1:ncpus=4:mem=40gb:scratch_local=40gb
#
# ## initialize the required application
# module add samtools-1.4 # Samtools v1.6 are broken! 
# module add python27-modules-gcc #required by multiQC
# #
# # Sorting and index creation for BAM files
# #
# # Requires samtools, multiQC
# #
# ############################################################################################
# ### Variables
# INPUT="/storage/brno3-cerit/home/tskalicky/Dis3L2/data/trimmed/repaired_pairs"
# OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/Dis3L2/data/trimmed/interleaved"
# MY_RAM=200 # Max RAM memory for Samtools sort
# THREADS=$PBS_NUM_PPN
# # Adding bbmap binaries into the PATH
# # Will work ONLY on NFS4 connected servers
# export PATH="/storage/brno3-cerit/home/tskalicky/tools/bbmap:$PATH"
# export PATH="/storage/brno3-cerit/home/tskalicky/anaconda2/bin:$PATH"
# #
# # Binaries
# SAMTOOLS=$(which samtools)
# REPAIR=$(which repair.sh)
# MULTIQC=$(which multiqc)
# PYTHON="/storage/brno3-cerit/home/tskalicky/anaconda2/bin/python"
# # Check the tools versions
# which $SAMTOOLS
# which $MULTIQC
# which $REPAIR
# which $PYTHON
# #
# ####################################################################################################
# # copy input data using SCRATCHDIR storage which is shared via NFSv4
# # clean the SCRATCH when job finishes (and data
# # are successfully copied out) or is killed
# # use cp -avr when copying directories
# cd $INPUT
# trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
# find "$INPUT" -maxdepth 2 -name "*.fastq.gz" -exec cp -vt "$SCRATCHDIR" {} +
# cd $SCRATCHDIR
# 
# if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
# echo "SCRATCHDIR path is:" $SCRATCHDIR
# echo "Following files were copied to scratch:"
# ls -Rc1
####################################################################################################
## Commands
grep -hiwE "^GGCGTCACTGTTGCGCTTCATAGACGCCGCGT" 