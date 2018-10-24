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
#module add piranha-1.2.1
module add python27-modules-gcc #required by multiQC
module add ctk-1.0.7
### Variables
# awk to print only fields from file containing mutations in unique tags .tag.uniq.mutation.txt
INPUTDIR="/home/skalda/Dropbox/CEITEC_lab/FTO/mapped_PCR_collapsed/"
FTO_CLIP="/storage/brno3-cerit/home/tskalicky/FTO/mapping/PCR_collapsed-genome/FTO_CLIP/BED"
FTO_INPUT="/storage/brno3-cerit/home/tskalicky/FTO/mapping/PCR_collapsed-genome/FTO_Input/BED"
FTO_KO1-3="/storage/brno3-cerit/home/tskalicky/FTO/mapping/PCR_collapsed-genome/FTO_KO1-3/BED"
FTO_WT1-3="/storage/brno3-cerit/home/tskalicky/FTO/mapping/PCR_collapsed-genome/FTO_WT1-3/BED"
THREADS=$PBS_NUM_PPN
APPENDIX="*.tag.uniq.mutation.txt"
####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
echo "SCRATCHDIR path is:" $SCRATCHDIR
echo "Following files/folders were copied to scratch:"
cp -av $FTO_CLIP/*$APPENDIX $FTO_INPUT/*$APPENDIX $FTO_WT1-3/*$APPENDIX $FTO_WT1-3/*$APPENDIX
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
####################################################################################################
# Binaries
#PIRANHA=$(which Piranha)
MULTIQC=$(which multiqc)


# Check if tools are installed
#which $PIRANHA
which $MULTIQC
####################################################################################################
mkdir CIMS
# Get specific types of mutations, such as deletions, substitutions, and insertions
# around the cross-linked mutation site.
for a in *.tag.uniq.mutation.txt
do
	FILE=$a
	FILENAME=${a%.*.*}
	printf "Counting number of mutations, deletions and insertions for file $FILE \n"
	awk '{if($9=="-") {print $0}}' $FILE | cut 1-6 >> $FILENAME".del.mutations.bed"
	awk '{if($9==">") {print $0}}' $FILE | cut 1-6 >> $FILENAME".sub.mutations.bed"
	awk '{if($9=="+") {print $0}}' $FILE | cut 1-6 >> $FILENAME".ins.mutations.bed"
	SUB=$(egrep -c $FILENAME".sub.mutations.bed")
	DEL=$(egrep -c $FILENAME".del.mutations.bed")
	INS=$(egrep -c $FILENAME".ins.mutations.bed")
	#awk '{print $8,$9,$10,$11}' $FILE >> $FILENAME".summary.temp.txt"
	#SUB=$(egrep -c '>' $FILENAME".summary.temp.txt")
	#DEL=$(egrep -c '-' $FILENAME".summary.temp.txt")
	#INS=$(egrep -c '+' $FILENAME".summary.temp.txt")
	printf "$SUB is summary number of mutattions in file $FILE \n" >> $FILENAME".mutation.summary.txt"
	printf "$DEL is summary number of deletions in file $FILE \n" >> $FILENAME".mutation.summary.txt"
	printf "$INS is summary number of insertions in file $FILE \n" >> $FILENAME".mutation.summary.txt"
done
#wait
wait
cat *".tag.uniq.del.mutations.bed" > $FILENAME".pool.tag.uniq.del.mutations.bed"

# It is always a good practice to look at the number of each type of mutation to see, 
# e.g., compare the relative abundance of deletions to insertions.
# For example, to get the number of deletions of different sizes:
awk '{print $3-$2}' Fox.pool.tag.uniq.del.bed | sort -n | uniq -c


printf "Moving resulting files to CIMS folder"
mv -v *".del.mutations.txt" $SCRATCHDIR/CIMS
mv -v *".sub.mutations.txt" $SCRATCHDIR/CIMS
mv -v *".ins.mutations.txt" $SCRATCHDIR/CIMS
mv -v *".mutation.summary.txt" $SCRATCHDIR/CIMS
