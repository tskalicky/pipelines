#!/bin/bash
#PBS -l select=1:ncpus=1:mem=15gb:scratch_local=50gb
#PBS -l walltime=04:00:00
#PBS -N 06_qc_prep
#
# Prepare input files for QC - picard and RSeQC
# You need to get conversion of chromosome names from NCBI -> Ensembl, especially 
#	for those with rRNA regions and change it at rRNA extraction part - you can 
#	find the information for example at NCBI genome page
#
# Requires python3, UCSC scripts (http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/), gtf2bed12 (http://allaboutbioinfo.blogspot.cz/2011/08/converting-cufflinks-gtf-predictions-to.html)
#
# EXPERIMENTAL!!!
#
# Note: Used GTF for generation of intervals for post-qc is different than the one used for alignment - added rRNAs (usually) from NCBI 
#
# TODO: automaticaly load chromosome names and make the NCBI->Ensembl conversion
####################################################################################################
### Variables
OUTPUT_DIR=/storage/brno2/home/opplatek/genomes/athaliana/ensembl37/TAIR10
OUTPUT_DIR=${OUTPUT_DIR}/postQC

GTF=/storage/brno2/home/opplatek/genomes/athaliana/ensembl37/TAIR10/Arabidopsis_thaliana.TAIR10.37.gtf.gz # Ensembl annotation (gtf)
GFF_NCBI=/storage/brno2/home/opplatek/genomes/athaliana/fastq_screen/ncbi/TAIR10/GCF_000001735.3_TAIR10_genomic.gff.gz # NCBI annotation (gff)
REF_SEQ=/storage/brno2/home/opplatek/genomes/athaliana/ensembl37/TAIR10/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz # Ensembl genome
EXTRA_RRNA=/storage/brno2/home/opplatek/genomes/human/ensembl87/GRCh38_rRNAadd.gtf # Additional manual additions to rRNA intervals - more info in GTF header
#HEADER=/storage/brno2/home/opplatek/genomes/human/ensembl87/Homo_sapiens.GRCh38.dna.primary_assembly.header # samtools view -H input.bam > header from the genome used for alignment and bam file

THREADS=$PBS_NUM_PPN

module add samtools-1.4
SAMTOOLS=$(which samtools)
module add cufflinks-2.2.1
GFFREAD=$(which gffread)
GTF_TO_BED12=/storage/brno2/home/opplatek/tools/scripts/gtf2bed12.py # http://allaboutbioinfo.blogspot.cz/2011/08/converting-cufflinks-gtf-predictions-to.html
UCSC_SCRIPTS=/storage/brno2/home/opplatek/tools/ucsc # http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
module add python-3.4.1-gcc

TMP=$SCRATCH

echo $SAMTOOLS
echo $GFFREAD
echo $GTF_TO_BED12
echo $UCSC_SCRIPTS
####################################################################################################
### Copying inputs
cp $REF_SEQ $SCRATCH/
cp $GTF $SCRATCH/
cp $GFF_NCBI $SCRATCH/
cp $EXTRA_RRNA $SCRATCH/
#cp $HEADER $SCRATCH/

REF_SEQ=$(basename $REF_SEQ)
GTF=$(basename $GTF)
GFF_NCBI=$(basename $GFF_NCBI)
EXTRA_RRNA=$(basename $EXTRA_RRNA)
#HEADER=$(basename $HEADER)

cd $SCRATCH/

unpigz -p $THREADS $REF_SEQ
REF_SEQ=${REF_SEQ%.gz}

unpigz -p $THREADS $GTF
GTF=${GTF%.gz}

unpigz -p $THREADS $GFF_NCBI
GFF_NCBI=${GFF_NCBI%.gz}

####################################################################################################
### Preparing
# Get the reference sequence header; this avoids getting if from SAM file; alternatives to http://slowkow.com/notes/ribosomal-rna/ or https://gist.github.com/adomingues/57e3c2f9f4aca8dfde7b
$SAMTOOLS faidx $REF_SEQ
cat $REF_SEQ.fai | awk 'OFS="\t" {print $1, $2}' > ${REF_SEQ%.fa*}.chrom.sizes # Get chrom sizes
#perl -lane 'print "\@SQ\tSN:$F[0]\tLN:$F[1]\tAS:TAIR10"' ${REF_SEQ%.fa*}.chrom.sizes | grep -v _ >> ${REF_SEQ%.fa*}.header # Sequence names and lengths. (Must be tab-delimited.)
# Don't know why there was "| grep -v "_ before!
perl -lane 'print "\@SQ\tSN:$F[0]\tLN:$F[1]\tAS:TAIR10"' ${REF_SEQ%.fa*}.chrom.sizes >> ${REF_SEQ%.fa*}.header # Sequence names and lengths. (Must be tab-delimited.)
HEADER=${REF_SEQ%.fa*}.header

### rRNA intervals
# Ensembl doesn't contain all rRNA annotations - we have to add them to the GTF - this is mainly applicable for 18S and 28S rRNA
# Ensembl annotations do not contain 18S and 28S rRNA annotation! You have to add it manually from Ensembl webserver/NCBI. Results in underestimation of rRNA content.

# Prepare rRNA intervals for Picard https://www.biostars.org/p/120145/
# NCBI (RefSeq) annotation contain all rRNAs - full rRNA annotation
#wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz
#unpigz -c -p $THREADS GRCh38_latest_genomic.gff.gz | grep "gbkey=rRNA" | grep -v "ribosomal RNA protein" > GRCh38_latest_genomic.rRNA.gff # All rRNAs
cat $GFF_NCBI | grep "gbkey=rRNA" | grep -v "ribosomal RNA protein" > ${GFF_NCBI%.*}.rRNA.gff # All rRNAs
#zcat GRCh38_latest_genomic.gff.gz | grep -P "\trRNA\t" > GRCh38_latest_genomic.rRNA.gff # All rRNAs

$GFFREAD ${GFF_NCBI%.*}.rRNA.gff -T -o ${GFF_NCBI%.*}.rRNA.gtf
cat ${GFF_NCBI%.*}.rRNA.gtf | sed 's/\texon\t/\tgene\t/g' > tmp_rRNA # Add gene lines
cat ${GFF_NCBI%.*}.rRNA.gtf | sed 's/\texon\t/\ttranscript\t/g' >> tmp_rRNA # Add transcript lines
cat $EXTRA_RRNA | grep -v "^#" >> tmp_rRNA # Manuall searched and added features
cat tmp_rRNA >> ${GFF_NCBI%.*}.rRNA.gtf
rm tmp_rRNA

# NCBI (RefSeq) uses different chromosomes than Ensembl - need to check when creating new rRNA intervals!
# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.3_TAIR10/GCF_000001735.3_TAIR10_assembly_report.txt for A.thaliana
sed -i 's/NC_003071\.7/2/g' ${GFF_NCBI%.*}.rRNA.gtf
sed -i 's/NC_003074\.8/3/g' ${GFF_NCBI%.*}.rRNA.gtf
sed -i 's/NC_001284\.2/Mt/g' ${GFF_NCBI%.*}.rRNA.gtf
sed -i 's/NC_000932\.1/Pt/g' ${GFF_NCBI%.*}.rRNA.gtf

sed -i -e 's/$/ gene_biotype \"rRNA\"; transcript_biotype \"rRNA\"; gene_source \"ncbi\"; transcript_source \"ncbi\";/' ${GFF_NCBI%.*}.rRNA.gtf # Add to end of each line

mv ${GFF_NCBI%.*}.rRNA.gtf ${GTF%.*}.rRNA.gtf
# Add also rRNA from Ensembl? There WILL be some overlap but it should be OK for Picard
grep "gene_biotype \"rRNA" $GTF >> ${GTF%.*}.rRNA.gtf # Ensembl only

# Intervals for rRNA transcripts.
cat $HEADER > ${GTF%.*}.rRNA.intervalListBody.txt
# cut -s -f 1,4,5,7,9 ${GTF%.*}.rRNA.gtf >> ${GTF%.*}.rRNA.intervalListBody.txt
cat ${GTF%.*}.rRNA.gtf | awk '$3 == "transcript"' | \
    cut -f1,4,5,7,9 | \
    perl -lane '
        /transcript_id "([^"]+)"/ or die "no transcript_id on $.";
        print join "\t", (@F[0,1,2,3], $1)
    ' | \
    sort -k1V -k2n -k3n \
>> ${GTF%.*}.rRNA.intervalListBody.txt

# Get rRNA which might be missing in Ensembl 
# You should check Ensembl if it contains 18S, 5.8S, 2S, 28S and 45S (45S might not be there even in NCBI)
grep rrn23S ${GTF%.*}.rRNA.gtf | grep "gene_source \"ncbi\"" > rRNA_add
grep rrn18 ${GTF%.*}.rRNA.gtf | grep "gene_source \"ncbi\"" >> rRNA_add
grep rrn16S ${GTF%.*}.rRNA.gtf | grep "gene_source \"ncbi\"" >> rRNA_add
cat $EXTRA_RRNA | grep -v "^#" >> rRNA_add
cat rRNA_add | sort | uniq > tmp; mv tmp rRNA_add
cat rRNA_add >> $GTF
(grep ^"#" $GTF; grep -v ^"#" $GTF | sort -T $TMP -k1,1 -k4,4n) > ${GTF%.*}.sorted.gtf
mv ${GTF%.*}.sorted.gtf $GTF
cat $GTF | uniq > tmp; mv tmp $GTF

# Prepare refFlat from gtf for Picard
$UCSC_SCRIPTS/gtfToGenePred -genePredExt $GTF ${GTF%.*}.refFlat.txt.tmp
paste <(cut -f 12 ${GTF%.*}.refFlat.txt.tmp) <(cut -f 1-10 ${GTF%.*}.refFlat.txt.tmp) > ${GTF%.*}.refFlat.txt # We have to add the gene name to the refFlat https://www.biostars.org/p/120145/; #cat ${GTF%.*}.refFlat.txt | awk 'BEGIN{FS="\t";} {OFS="\t";} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > ${GTF%.*}.refFlat.txt # Same as above but written in a different way
rm ${GTF%.*}.refFlat.txt.tmp

# Prepare bed12 for RSeQC
sed '/^#/d' $GTF > $GTF.tmp # Temporary remove # comments in GTF if present
python3 $GTF_TO_BED12 $GTF.tmp > ${GTF%.*}.genes.bed12 # Be careful about python3, doesn't work with python2

# sed -i 's/^chr//g' ${GTF%.*}.genes.bed12 # If downloaded from RSeQC by alignment done to Ensembl; TODO check mitochondrial genome

mv $SCRATCH/$GTF $SCRATCH/${GTF%.*}.rRNAadd.gtf
pigz -p $THREADS $SCRATCH/${GTF%.*}.rRNAadd.gtf
rm ${GTF}.tmp
####################################################################################################
### Copying outputs
rm $SCRATCH/$HEADER
rm $SCRATCH/$GFF_NCBI
rm $SCRATCH/$REF_SEQ*

mkdir -p $OUTPUT_DIR

cp -r $SCRATCH/* $OUTPUT_DIR/

rm -r $SCRATCH/*