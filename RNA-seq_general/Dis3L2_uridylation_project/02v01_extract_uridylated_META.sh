#!/bin/bash
#!/bin/bash
#PBS -l walltime=96:0:0 
#PBS -q uv@wagap-pro.cerit-sc.cz
#PBS -l select=1:ncpus=20:mem=150gb:scratch_local=450gb
#PBS -j oe
#PBS -N Dis3L2_spikein_mapped_extract
#
# For 1 library you need cca 35GB RAM
#
## initialize the required application
module add samtools-1.4 # Samtools v1.6 are broken! 
module add python27-modules-gcc #required by multiQC
############################################################################################
### Variables
THREADS=$PBS_NUM_PPN
INPUT_DIR="/storage/brno3-cerit/home/tskalicky/Dis3L2/Dasa_spikein/mapping/human_genome/star/results_job_1052371/alignment/genome"
OUTPUT_DIR="/storage/brno3-cerit/home/tskalicky/Dis3L2/Dasa_spikein/mapping/human_genome/star/results_job_1052371/alignment/genome/mapped_extracted"
EXTRACT="/home/skalda/ownCloud/git/pipelines/RNA-seq_general/Dis3L2_uridylation_project/extract_uridylated.py"
#
APPENDIX1=".bam"
APPENDIX2=".fastq"
APPENDIX3=".fa"
# Binaries
SAMTOOLS=$(which samtools)
MULTIQC=$(which multiqc)
# Check the tools versions
which $SAMTOOLS
which $MULTIQC
which python
#
####################################################################################################
# copy input data using SCRATCHDIR storage which is shared via NFSv4
# clean the SCRATCH when job finishes (and data
# are successfully copied out) or is killed
# use cp -avr when copying directories
cd $INPUT_DIR
trap 'clean_scratch' TERM EXIT # sets up scratch cleaning in case an error occurs
find "$INPUT_DIR" -maxdepth 1 -name "*.fa.gz" -exec cp -vt "$SCRATCHDIR" {} +
# find "$INPUT_DIR" -maxdepth 2 -name "*.bam" -exec cp -vt "$SCRATCHDIR" {} + # this is not working at every occasion!!
cd $SCRATCHDIR

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch directory is not created!" 1>&2; exit 1; fi #checks if scratch directory is created
echo "SCRATCHDIR path is:" $SCRATCHDIR
echo "Following files were copied to scratch:"
ls -Rc1
####################################################################################################
####################################################################################################
## Commands
cd $INPUT_DIR
mkdir -p $OUTPUT_DIR
for a in .fa.gz
do
	SAMPLE=$a
	SAMPLENAME=${a%.*.*}
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Start unpacking reads from file $SAMPLE"
	unpigz -vp $THREADS $SAMPLE
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Start extracting uridylated reads from file $SAMPLE"
	python extract_uridylated_out_fasta.py $SAMPLENAME".fa"
	pigz -vp $THREADS "*.fa"
	date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
	echo "Start packing output files from file $SAMPLE"
done
wait
