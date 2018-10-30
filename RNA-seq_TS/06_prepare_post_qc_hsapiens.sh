#!/bin/bash
#!/bin/bash
#PBS -l walltime=24:0:0 
#PBS -q uv@wagap-pro.cerit-sc.cz
#PBS -l select=1:ncpus=2:mem=20gb:scratch_local=100gb
#PBS -j oe
#PBS -N 06_qc_prep_Dis3L2_spikeIn
#
# Prepare input files for QC - picard and RSeQC
# You need to get conversion of chromosome names from NCBI -> Ensembl, especially 
#	for those with rRNA regions and change it at rRNA extraction part - you can 
#	find the information for example at NCBI genome page
#
# Requires python3, UCSC scripts (http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/), gtf2bed12 (http://allaboutbioinfo.blogspot.cz/2011/08/converting-cufflinks-gtf-predictions-to.html)
#
# Note: Used GTF for generation of intervals for post-qc is different than the one used for alignment - added rRNAs (usually) from NCBI 
#
# TODO: automaticaly load chromosome names and make the NCBI->Ensembl conversion
# TODO: You might consider using picardmetrics https://slowkow.github.io/picardmetrics/
####################################################################################################
### Variables
OUTPUT_DIR=/mnt/storage-brno3-cerit/nfs4/home/tskalicky/Dis3L2/Dasa_spikein  # Main directory with references and indexes
OUTPUT_DIR=${OUTPUT_DIR}/post_mapping_QC

GTF="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl94/Homo_sapiens.GRCh38.94.gtf.gz" # Ensembl annotation (gtf)
GFF_NCBI="/storage/brno3-cerit/home/tskalicky/genomes/human/NCBI/GRCh38.p12_genomic.gff.gz" # NCBI annotation (gff)
REF_SEQ="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl94/Homo_sapiens.GRCh38.94.dna.primary_assembly.fa.gz" # Ensembl genome
EXTRA_TRNA="/storage/brno3-cerit/home/tskalicky/genomes/human/ensembl91/gencode.v27.tRNAs_modif_for_Pickard_RSeQC.gtf"
EXTRA_RRNA=/storage/brno3-cerit/home/tskalicky/genomes/human/NCBI/GRCh38_rRNAadd.gtf # Additional manual additions to rRNA intervals - more info in GTF header
#HEADER=/storage/brno2/home/opplatek/genomes/human/ensembl87/Homo_sapiens.GRCh38.dna.primary_assembly.header # samtools view -H input.bam > header from the genome used for alignment and bam file

THREADS=$PBS_NUM_PPN

module add samtools-1.4
SAMTOOLS=$(which samtools)
module add cufflinks-2.2.1
GFFREAD=$(which gffread)
GTF_TO_BED12=/storage/brno3-cerit/home/tskalicky/tools/gtf2bed12.py # http://allaboutbioinfo.blogspot.cz/2011/08/converting-cufflinks-gtf-predictions-to.html
UCSC_SCRIPTS=/storage/brno3-cerit/home/tskalicky/tools/ucsc # http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
module add python-3.4.1-gcc

TMP=$SCRATCHDIR

echo $SAMTOOLS
echo $GFFREAD
echo $GTF_TO_BED12
echo $UCSC_SCRIPTS
####################################################################################################
### Copying inputs
cp $REF_SEQ $SCRATCHDIR/
cp $GTF $SCRATCHDIR/
cp $GFF_NCBI $SCRATCHDIR/
cp $EXTRA_TRNA $SCRATCHDIR/
#cp $HEADER $SCRATCHDIR/

REF_SEQ=$(basename $REF_SEQ)
GTF=$(basename $GTF)
GFF_NCBI=$(basename $GFF_NCBI)
EXTRA_TRNA=$(basename $EXTRA_TRNA)
#HEADER=$(basename $HEADER)

cd $SCRATCHDIR/

unpigz -p $THREADS $REF_SEQ
REF_SEQ=${REF_SEQ%.gz}

unpigz -p $THREADS $GTF
GTF=${GTF%.gz}

unpigz -p $THREADS $GFF_NCBI
GFF_NCBI=${GFF_NCBI%.gz}

unpigz -p $THREADS $EXTRA_TRNA
EXTRA_TRNA=${EXTRA_TRNA%.gz}

####################################################################################################
### Preparing
# Get the reference sequence header; this avoids getting if from SAM file; alternatives to http://slowkow.com/notes/ribosomal-rna/ or https://gist.github.com/adomingues/57e3c2f9f4aca8dfde7b
$SAMTOOLS faidx $REF_SEQ
cat $REF_SEQ.fai | awk 'OFS="\t" {print $1, $2}' > ${REF_SEQ%.fa*}.chrom.sizes # Get chrom sizes
#perl -lane 'print "\@SQ\tSN:$F[0]\tLN:$F[1]\tAS:GRCh38"' ${REF_SEQ%.fa*}.chrom.sizes | grep -v _ >> ${REF_SEQ%.fa*}.header # Sequence names and lengths. (Must be tab-delimited.)
# Don't know why there was "| grep -v "_ before!
perl -lane 'print "\@SQ\tSN:$F[0]\tLN:$F[1]\tAS:GRCh38"' ${REF_SEQ%.fa*}.chrom.sizes >> ${REF_SEQ%.fa*}.header # Sequence names and lengths. (Must be tab-delimited.)
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
if [ -n "${EXTRA_TRNA}" ]; then cat $EXTRA_TRNA | grep -v "^#" >> tmp_rRNA; fi # Manuall searched and added features
cat tmp_rRNA >> ${GFF_NCBI%.*}.rRNA.gtf
rm tmp_rRNA

# NCBI (RefSeq) uses different chromosomes than Ensembl - need to check when creating new rRNA intervals!
# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.37_GRCh38.p11/GCF_000001405.37_GRCh38.p11_assembly_report.txt for human
sed -i 's/NC_000001\.11/1/g' ${GFF_NCBI%.*}.rRNA.gtf
sed -i 's/NC_012920\.1/MT/g' ${GFF_NCBI%.*}.rRNA.gtf
sed -i 's/NT_167214\.1/GL000220.1/g' ${GFF_NCBI%.*}.rRNA.gtf

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
grep RNA45S ${GTF%.*}.rRNA.gtf | grep "gene_source \"ncbi\"" > rRNA_add
grep RNA28S ${GTF%.*}.rRNA.gtf | grep "gene_source \"ncbi\"" >> rRNA_add
grep RNA18S ${GTF%.*}.rRNA.gtf | grep "gene_source \"ncbi\"" >> rRNA_add
if [ -n "${EXTRA_TRNA}" ]; then cat $EXTRA_TRNA | grep -v "^#" >> rRNA_add; fi
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

mv $SCRATCHDIR/$GTF $SCRATCHDIR/${GTF%.*}.rRNAadd.gtf
pigz -p $THREADS $SCRATCHDIR/${GTF%.*}.rRNAadd.gtf
rm ${GTF}.tmp
####################################################################################################
### Copying outputs
rm $SCRATCHDIR/$HEADER
rm $SCRATCHDIR/$GFF_NCBI
rm $SCRATCHDIR/$REF_SEQ*

mkdir -p $OUTPUT_DIR

cp -r $SCRATCHDIR/* $OUTPUT_DIR/

rm -r $SCRATCHDIR/*
