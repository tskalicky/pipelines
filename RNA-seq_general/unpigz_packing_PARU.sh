#!/bin/bash
#PBS -l nodes=1:ppn=20
#PBS -l mem=80gb
#PBS -k oe
#PBS -N packing_Dis3L2_SpikeIn_unmapped_reads
#PBS -M tomas.skalicky@seznam.cz
#PBS -m abe
#
#PARU KRTECEK server is using TORQUE scheduling system !!!
#
## initialize the required application
INPUT_DIR="/home/users/tskalicky/CEITEC/Dis3L2/Dasa_SpikeIn/mapping/Scerevisiae_genome/unmapped_reads"
#
# Set number of CPU
THREADS=$PBS_NUM_PPN
# Set MAX ram
MY_RAM="80"
RAM=$[$MY_RAM-1] # lower it by one to leave some space
MEM_LIMIT=$RAM # $[$RAM/$THREADS] # In case we need GB per thread
#
export PATH="/usr/bin/:$PATH"
# Binaries
UNPIGZ=$(which unpigz)

which $UNPIGZ
#
## Commands
cd $INPUT_DIR
for a in *.fq;
do
	FILE="$a"
	FILENAME="${a%.*}"
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Packing file $FILE"
	gzip -v $FILE &
done
wait

