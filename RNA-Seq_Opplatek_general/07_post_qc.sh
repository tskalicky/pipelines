#!/bin/bash
#PBS -l select=1:ncpus=6:mem=28gb:scratch_local=200gb
#PBS -l walltime=24:00:00
#PBS -q default
#PBS -N 07_qc_alignment_run
#
# QC of RNA-Seq experiment (after alignment)
#
# Requires Picard (java), samtools, featureCounts, dupRadar, Preseq, R, multiqc, prepared references (see 06_prepare_qc_alignment.sh)
#
# Note: Preseq R script might not work for PE experiments (probable naming errors - now based on "R1" naming)
# Note: Check dupRadar script - it must be pointing to correct path to Rscript
# Note: Used GTF for generation of intervals for post-qc is different than the one used for alignment - added rRNAs (usually) from NCBI 
############################################################################################
### Variables
INPUT_DIR=/storage/brno3-cerit/home/opplatek/biocore/karin_rnaseq/results/alignment/genome # Pointing to genome/ mapping directory
OUTPUT_DIR=/storage/brno3-cerit/home/opplatek/biocore/karin_rnaseq
OUTPUT_DIR=${OUTPUT_DIR}/results/qc

COUNT_OVER="exon" # [exon, three_prime_utr] what feature to use for the summary? For QuantSeq FWD is should be "exon", for QuantSeq REV it should be "three_prime_utr"
PAIRED="false" # True or false for PE reads
STRANDNESS="rev" # [fwd, rev, none] strandedness of read. For QuantSeq FWD is should be "rev" 

RPKM_SATUR="true" # [true, false] Do you want to run RPKM-saturation analysis from RSeQC? It takes quite a long time so it is off by default

# CHECK STRAND_RULE at RPKM_saturation.py
#STRAND_RULE="'1++,1--,2+-,2-+'" # THIS ASSIGNMENT TO VARIABLE DOESN'T WORK For RSeQC - strand orientation; can be taken from Infer experiment; '++,--' in case of QuantSeq FWD; '1++,1--,2+-,2-+' in case of ScriptSeq; single quotes HAS TO BE THERE

APPENDIX=".sortedByCoord.out.bam"

THREADS=$PBS_NUM_PPN

GTF=/storage/brno2/home/opplatek/genomes/human/ensembl87/postQC/Homo_sapiens.GRCh38.87.rRNAadd.gtf.gz # Modified GTF with added rRNA from NCBI 
REF_FLAT=/storage/brno2/home/opplatek/genomes/human/ensembl87/postQC/Homo_sapiens.GRCh38.87.refFlat.txt # gtf to refFlat http://seqanswers.com/forums/showthread.php?t=50450 BUT we need to modify it to add gene name otherwise it causes error - see lines above
RIBOSOMAL_INT=/storage/brno2/home/opplatek/genomes/human/ensembl87/postQC/Homo_sapiens.GRCh38.87.rRNA.intervalListBody.txt # interval file information 
REF_BED=/storage/brno2/home/opplatek/genomes/human/ensembl87/postQC/Homo_sapiens.GRCh38.87.genes.bed12 # BED12 format for RSeQC

module add jdk-8
#module add picard-2.8.1
#PICARD_RUN="java -jar -Xmx8g $PICARD281"
PICARD_RUN="java -jar -Xmx8g /storage/brno2/home/opplatek/tools/picard/picard-2.12.1.jar" # Testing
module add samtools-1.4 # Needed if BAM files are not indexes for RSeQC
SAMTOOLS=$(which samtools)
FEATURE_COUNTS=/storage/brno2/home/opplatek/tools/subread-1.5.2-Linux-x86_64/bin/featureCounts
DUPRADAR=/storage/brno2/home/opplatek/tools/scripts/dupRadar.r # DupRadar R script
DUPRADAR_INSTALL=/auto/brno2/home/opplatek/R/x86_64-pc-linux-gnu-library/3.4 # dupRadar installation folder 
PRESEQ=/storage/brno2/home/opplatek/tools/preseq_v2.0.2/preseq
PRESEQ_R=/storage/brno2/home/opplatek/tools/scripts/preseq.r # Preseq visualization script
module add python27-modules-gcc # For MultiQC
#PYTHONPATH=$PYTHONPATH:/storage/brno2/home/opplatek/tools/MultiQC-1.0/lib/python2.7/site-packages/
#MULTIQC=/storage/brno2/home/opplatek/tools/MultiQC-1.0/bin/multiqc
MULTIQC=$(which multiqc)
FEATURE_COUNT_R=/storage/brno2/home/opplatek/scripts/featureCounts_biotypes.r # Feature gene biotype visualization

# Check the tools versions
echo $PICARD_RUN
echo $SAMTOOLS
echo $FEATURE_COUNTS
echo $FEATURE_COUNT_R
echo $DUPRADAR
echo $DUPRADAR_INSTALL
echo $PRESEQ
echo $PRESEQ_R
echo $MULTIQC
echo $FEATURE_COUNT_R

####################################################################################################
### Copying inputs and preparing
cp $GTF $SCRATCH/
cp $REF_FLAT $SCRATCH/
cp $RIBOSOMAL_INT $SCRATCH/
cp $REF_BED $SCRATCH/

cp $INPUT_DIR/*$APPENDIX $SCRATCH/

GTF=$(basename $GTF)
REF_FLAT=$(basename $REF_FLAT)
RIBOSOMAL_INT=$(basename $RIBOSOMAL_INT)
REF_BED=$(basename $REF_BED)

cd $SCRATCH/

unpigz -p $THREADS $GTF
GTF=${GTF%.gz}

# Set extra flags for Picard, featureCounts and dupRadar
msg=""
extra_flags_picard="FIRST_READ_TRANSCRIPTION_STRAND" # Default for SE as in the manual
extra_flags_preseq=""
extra_flags_feature=""
extra_flags_dupRadar="single" # Default for SE
if [ "$PAIRED" == "true" ]; then
	extra_flags_preseq="-pe -l 1500000" # For Preseq
	extra_flags_feature="-p" # For featureCounts
	extra_flags_dupRadar="paired" # for RSEM
	msg="PE experiment"
else
	msg="SE experiment"
fi
# Set if it is stranded or not (fwd, rev, none)
if [ "$STRANDNESS" == "fwd" ] || [ "$STRANDNESS" == "+" ] || [ "$STRANDNESS" == "ScriptSeq" ]; then
	extra_flags_picard="FIRST_READ_TRANSCRIPTION_STRAND"	
	extra_flags_feature="$extra_flags_feature -s 1" # For featureCounts
	extra_flags_dupRadar="$extra_flags_dupRadar 1" # for RSEM
	msg="$msg forward stranded"
elif [ "$STRANDNESS" == "rev" ] || [ "$STRANDNESS" == "-" ] || [ "$STRANDNESS" == "TruSeq" ]; then
	extra_flags_picard="SECOND_READ_TRANSCRIPTION_STRAND"
	extra_flags_feature="$extra_flags_feature -s 2" # For featureCounts
	extra_flags_dupRadar="$extra_flags_dupRadar 2" # for RSEM
	msg="$msg reverse stranded"
else
	extra_flags_picard="NONE"
	extra_flags_feature="$extra_flags_feature -s 0" # For featureCounts
	extra_flags_dupRadar="$extra_flags_dupRadar 0" # for RSEM
	msg="$msg unstranded"
fi
echo "Running as $msg"

# Set extra flags for RSeQC
extra_flags_rseqc="" # For RPKM_saturation
if [ "$PAIRED" == "true" ]; then
	msg="PE experiment"
	if [ "$STRANDNESS" == "fwd" ]; then
		extra_flags_rseqc=--strand='1++,1--,2+-,2-+'
		msg="$msg forward strand"
	elif [ "$STRANDNESS" == "rev" ]; then
		extra_flags_rseqc=--strand='1+-,1-+,2++,2--'
		msg="$msg reverse strand"
	else
		msg="$msg unstranded"
	fi
else
	msg="SE experiment"
	if [ "$STRANDNESS" == "fwd" ]; then
		extra_flags_rseqc=--strand='++,--'
		msg="$msg forward strand"
	elif [ "$STRANDNESS" == "rev" ]; then
		extra_flags_rseqc=--strand='+-,-+'
		msg="$msg reverse strand"
	else
		msg="$msg unstranded"
	fi
fi
echo "Running RSeQC RPKM_saturation run as $msg"

for i in *$APPENDIX
do
	$SAMTOOLS index -@ $THREADS $i
done

####################################################################################################
### Alignment QC
### QC by Picard
mkdir $SCRATCH/picard

for i in *$APPENDIX
do
	$PICARD_RUN CollectRnaSeqMetrics \
		  I=$i \
		  O=$SCRATCH/picard/${i%.*}.output.RNA_Metrics.txt \
		  REF_FLAT=$REF_FLAT \
		  STRAND=$extra_flags_picard \
		  RIBOSOMAL_INTERVALS=$RIBOSOMAL_INT \
		  MINIMUM_LENGTH=200 \
		  CHART_OUTPUT=${i%.*}.npc.pdf \
		  VALIDATION_STRINGENCY=LENIENT # This is present only because older versions of STAR has some issue with "Not primary alignment flag should not be set for unmapped read" https://github.com/alexdobin/STAR/issues/177; this gives us only warning instead of error
		  #STOP_AFTER=50000000
done

mv $SCRATCH/*.pdf $SCRATCH/picard

mkdir -p $OUTPUT_DIR

#pdfunite $SCRATCH/picard/*.out.npc.pdf $SCRATCH/picard/all.out.npc.pdf
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=$SCRATCH/picard/all.out.npc.pdf $SCRATCH/picard/*.out.npc.pdf # alternative to pdfunite

cp -r $SCRATCH/picard $OUTPUT_DIR/ &

### QC by featureCounts 
# Summary of mapped gene biotypes
# Get gene biotype counts
mkdir -p $SCRATCH/featureCounts

$FEATURE_COUNTS -t exon -g gene_biotype -O -M --fraction $extra_flags_feature -T $THREADS -a $GTF -o $SCRATCH/featureCounts/$STRANDNESS.biotype.featureCounts *.bam # should we keep --fraction, if it's not there rRNA content is overestimated

# Get just biotype counts
for i in $SCRATCH/featureCounts/*.biotype.featureCounts
do
        cut -f 1,7- $i > ${i%.*}.biotype_counts.txt # Gets first columns and then all columns startgin from 7th column till the end
done

cd $SCRATCH/featureCounts

R --no-save < $FEATURE_COUNT_R # Visualize the results
module rm R-3.4.0-gcc

$MULTIQC -o $SCRATCH/featureCounts $SCRATCH/featureCounts/

cd $SCRATCH/

cp -r $SCRATCH/featureCounts $OUTPUT_DIR/ &

### QC by RSeQC part 1/2
module add RSeQC-2.6.1

mkdir $SCRATCH/RSeQC

for i in *$APPENDIX
do
	echo "Processing file $i"

	# Distribution of reads along the genome
	echo "Read distribution"
	read_distribution.py -i $i -r $SCRATCH/$REF_BED  > $SCRATCH/RSeQC/${i%.bam*}.read_distribution.txt

	# Junction saturation - is the sequencing depth enough for alt. splicing analysis?
	echo "Junction saturation"
	junction_saturation.py -i $i -r $SCRATCH/$REF_BED -o $SCRATCH/RSeQC/${i%.bam*}.junction_saturation

	# Junction annotation - compare detected splice junctions to reference gene model
	echo "Junction annotation"
	junction_annotation.py -i $i -r $SCRATCH/$REF_BED -o $SCRATCH/RSeQC/${i%.bam*}.junction_annotation

	# bam stats - number of reads mapped, not mapped, etc.
	echo "Bam stat"
	bam_stat.py -i $i 2> $SCRATCH/RSeQC/${i%.bam*}.bam_stat.txt

	# Check strandeness of the experiment
	echo "Infer experiment"
	infer_experiment.py -r $SCRATCH/$REF_BED -i $i > $SCRATCH/RSeQC/${i%.bam*}.infer_experiment.txt

	# Check read duplication
	echo "Read duplication"
	read_duplication.py -i $i -o $SCRATCH/RSeQC/${i%.bam*}.read_duplication

	# Inner distance - fragment size
	echo "Inner distance"
	inner_distance.py -r $SCRATCH/$REF_BED -i $i -o $SCRATCH/RSeQC/${i%.bam*}.inner.distance
done

module rm RSeQC-2.6.1

module add python27-modules-gcc # For MultiQC
PYTHONPATH=$PYTHONPATH:/storage/brno2/home/opplatek/tools/MultiQC-1.0/lib/python2.7/site-packages/

$MULTIQC -o $SCRATCH/RSeQC $SCRATCH/RSeQC/

cp -r $SCRATCH/RSeQC $OUTPUT_DIR/ &

### DupRadar duplicates check
#Mark duplicates using Picard and check biological/technical duplication by dupRadar
module add picard-2.8.1

mkdir -p $SCRATCH/picard_dupl
mkdir -p $SCRATCH/picard

for i in *$APPENDIX
do
	$PICARD_RUN MarkDuplicates \
	    INPUT=$i \
	    OUTPUT=$SCRATCH/picard_dupl/${i%.*}.markDups.bam \
	    METRICS_FILE=$SCRATCH/picard/${i%.*}.markDups_metrics.txt \
	    REMOVE_DUPLICATES=false \
	    ASSUME_SORTED=true \
	    PROGRAM_RECORD_ID=null \
	    VALIDATION_STRINGENCY=LENIENT
done

cd $SCRATCH/picard_dupl
mkdir $SCRATCH/dupRadar

module add R-3.4.0-gcc

for i in *.markDups.bam
do
	echo "Processing $i by dupRadar"
	$DUPRADAR $i $SCRATCH/$GTF $extra_flags_dupRadar $DUPRADAR_INSTALL
done

#pdfunite *_duprateExpBoxplot.pdf $SCRATCH/dupRadar/all.duprateExpBoxplot.pdf
#pdfunite *_expressionHist.pdf $SCRATCH/dupRadar/all.expressionHist.pdf
#pdfunite *_duprateExpDens.pdf $SCRATCH/dupRadar/all.duprateExpDens.pdf
#pdfunite *_multimapPerGene.pdf $SCRATCH/dupRadar/all.multimapPerGene.pdf
#pdfunite *_readDist.pdf $SCRATCH/dupRadar/all.readDist.pdf
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=$SCRATCH/dupRadar/all.duprateExpBoxplot.pdf *_duprateExpBoxplot.pdf # alternative to pdfunite
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=$SCRATCH/dupRadar/all.expressionHist.pdf *_expressionHist.pdf # alternative to pdfunite
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=$SCRATCH/dupRadar/all.duprateExpDens.pdf *_duprateExpDens.pdf # alternative to pdfunite
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=$SCRATCH/dupRadar/all.multimapPerGene.pdf *_multimapPerGene.pdf # alternative to pdfunite
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=$SCRATCH/dupRadar/all.readDist.pdf *_readDist.pdf # alternative to pdfunite

#rm $SCRATCH/dupRadar/*_duprateExpBoxplot.pdf $SCRATCH/dupRadar/*_expressionHist.pdf $SCRATCH/dupRadar/*_duprateExpDens.pdf $SCRATCH/dupRadar/*_multimapPerGene.pdf $SCRATCH/dupRadar/*_readDist.pdf

mv $SCRATCH/picard_dupl/*.duprateExpBoxplot.pdf $SCRATCH/dupRadar/
mv $SCRATCH/picard_dupl/*.expressionHist.pdf $SCRATCH/dupRadar/
mv $SCRATCH/picard_dupl/*.duprateExpDens.pdf $SCRATCH/dupRadar/
mv $SCRATCH/picard_dupl/*.multimapPerGene.pdf $SCRATCH/dupRadar/
mv $SCRATCH/picard_dupl/*.readDist.pdf $SCRATCH/dupRadar/
mv $SCRATCH/picard_dupl/*.txt $SCRATCH/dupRadar/

$MULTIQC -o $SCRATCH/picard $SCRATCH/picard/ &

for i in $SCRATCH/picard_dupl/*
do 
	pigz -p $THREADS $i
done

cp -r $SCRATCH/picard/* $OUTPUT_DIR/picard/
cp -r $SCRATCH/dupRadar $OUTPUT_DIR/
cp -r $SCRATCH/picard_dupl $OUTPUT_DIR/

mv $SCRATCH/*.bai $INPUT_DIR/

### QC by Preseq
# Run Preseq complexity library test
mkdir -p $SCRATCH/preseq

for i in *$APPENDIX
do
	echo "Preseq - processing $i"
	$PRESEQ c_curve -B $extra_flags_preseq -o $SCRATCH/preseq/${i%.*}.estimates.txt $i # c_curve    generate complexity curve for a library; if paired-end add -P
	$PRESEQ lc_extrap -B $extra_flags_preseq -o $SCRATCH/preseq/${i%.*}.yield_estimates.txt $i # lc_extrap  predict the yield for future experiments; if paired-end add -P
done

cd $SCRATCH/preseq

module add R-3.4.0-gcc
R --no-save < $PRESEQ_R # Visualize the results
module rm R-3.4.0-gcc

PYTHONPATH=$PYTHONPATH:/storage/brno2/home/opplatek/tools/MultiQC-1.0/lib/python2.7/site-packages/
$MULTIQC -o $SCRATCH/preseq $SCRATCH/preseq/

cd $SCRATCH/

mv $SCRATCH/preseq $OUTPUT_DIR/ &

### QC by RSeQC part 2/2
# This take very long time (RPKM_saturation)!
cd $SCRATCH/

if [ "$RPKM_SATUR" == "true" ]; then
	echo "Running RPKM_saturation"
	module add RSeQC-2.6.1

	for i in *$APPENDIX
	do
		# Check saturation of expressed genes in 4 expression quartiles http://rseqc.sourceforge.net/#rpkm-saturation-py
		# Need to add stranded information from Infer experiment
		echo "RPKM saturation for $i"
		#	RPKM_saturation.py -r $SCRATCH/$REF_BED -i $i --strand='++,–' -o $SCRATCH/RSeQC/${i%.bam*}.RPKM_saturation # -d equals to --strand= -> Need to get this from Infer experiment or from the experimental procedures so it might fail for the first time
		RPKM_saturation.py -r $SCRATCH/$REF_BED -i $i $extra_flags_rseqc -o $SCRATCH/RSeQC/${i%.bam*}.RPKM_saturation # -d equals to --strand= -> Need to get this from Infer experiment or from the experimental procedures so it might fail for the first time
	done

	module rm RSeQC-2.6.1

	ls $SCRATCH/RSeQC/*.RPKM_saturation.saturation.pdf > $SCRATCH/RSeQC/all.RPKM_saturation.saturation.txt # Get names for the plots bellow
	#pdfunite $SCRATCH/RSeQC/*.RPKM_saturation.saturation.pdf $SCRATCH/RSeQC/all.RPKM_saturation.saturation.pdf
	gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=$SCRATCH/RSeQC/all.RPKM_saturation.saturation.pdf $SCRATCH/RSeQC/*.RPKM_saturation.saturation.pdf # alternative to pdfunite

	mkdir $SCRATCH/RSeQC/RPKM_saturation
	mv $SCRATCH/RSeQC/*.RPKM_saturation.saturation.r $SCRATCH/RSeQC/RPKM_saturation/
	mv $SCRATCH/RSeQC/*.RPKM_saturation.eRPKM.xls $SCRATCH/RSeQC/RPKM_saturation/
	mv $SCRATCH/RSeQC/*.RPKM_saturation.rawCount.xls $SCRATCH/RSeQC/RPKM_saturation/

	for i in $SCRATCH/RSeQC/RPKM_saturation/*
	do
		pigz -p $THREADS $i
	done
else
	echo "Not running RPKM_saturation"
fi

# Coverage of genes - extremely slow, better to do it with picard
#module add RSeQC-2.6.1
#INPUT_BAMS=$(ls -m *bam | tr -d \\n)
#echo "Gene body coverage"
#geneBody_coverage.py -i $INPUT_BAMS -r $SCRATCH/$REF_BED -l 100 -f pdf -o $SCRATCH/RSeQC/mrna.all
#module rm RSeQC-2.6.1

cp -r $SCRATCH/RSeQC $OUTPUT_DIR/
####################################################################################################
### Clean
rm -r $SCRATCH/*