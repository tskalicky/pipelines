#!/bin/bash
#
# variables
ENSEMBL_DB="$HOME/ownCloud/CEITEC_lab/genomes/human/ensembl91"
GeneCodeDB="$HOME/ownCloud/CEITEC_lab/genomes/annotation/my_annot_DB/genecode/"
tRNA_ANNOT="$HOME/ownCloud/CEITEC_lab/genomes/annotationDB/gencode/gencode.v27.tRNAs.gtf.gz"
OUTPUT="$HOME/ownCloud/CEITEC_lab/genomes/annotation/my_annot_DB/genecode"
BED_OUTPUT="$HOME/ownCloud/CEITEC_lab/genomes/annotation/my_annot_DB/bed"
#
# Check if tools are installed
BEDTOOLS=$(which bedtools)
SORTBED=$(which sortBed)
MERGEBED=$(which mergeBed)
COMPLEMENTBED=$(which complementBed)
#
which $BEDTOOLS
which $SORTBED
which $MERGEBED
which $COMPLEMENTBED
#
## commands
#get latest annotation file from GeneCode
cd $GeneCodeDB
v=27
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_$v/gencode.v$v.annotation.gtf.gz
#what does this file look like?
zcat gencode.v27.primary_assembly.annotation.gtf.gz | grep -v "^#" | cut -f3 | sort | uniq -c | sort -k1rn > $HOME/ownCloud/CEITEC_lab/genomes/annotation/my_annot_DB/genecode/genecode_RNA_features.txt
#
head -n20 gencode.v27.primary_assembly.annotation.gtf | awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5}'
sed -er 's/gene_name \"//' -er 's/\"\t/1000\t/' GRCh38.91.five_prime_utr.sorted.bed >GRCh38.91.five_prime_utr.2.sorted.bed

#if file in .gz format use zcat
##zcat $ENSEMBL_DB/Homo_sapiens.GRCh38.91.sorted.gtf.gz | grep -v "^#" | cut -f3 | sort | uniq -c | sort -k1rn > $OUTPUT/ENSEMBL_RNA_features.txt
# cat /home/skalda/ownCloud/CEITEC_lab/genomes/human/ensembl91/Homo_sapiens.GRCh38.91.sorted.gtf \
# | grep -v "^#" | cut -f3 | sort | uniq -c | sort -k1rn \
# > /home/skalda/ownCloud/CEITEC_lab/genomes/annotation/my_annot_DB/ENSEMBL_RNA_features.txt
#
# /home/skalda/ownCloud/CEITEC_lab/genomes/human/ensembl91/Homo_sapiens.GRCh38.91.sorted.gtf  | grep -v "^#" | cut -f3 | sort | uniq -c | sort -k1rn >$OUTPUT/ ENSEMBL_RNA_features.txt

# Creating annotation file subset 
#cd $FULL_ANNOTATION
#grep -hwF "gene" $ENSEMBL_DB/Homo_sapiens.GRCh38.91.sorted.gtf >$OUTPUT/Homo_sapiens.GRCh38.91.sorted.genes_only.gtf
##
#grep -hwF "rRNA" $ENSEMBL_DB/Homo_sapiens.GRCh38.91.sorted.genes_only.gtf >$OUTPUT/Homo_sapiens.GRCh38.91.sorted.rRNA_only.gtf
## tRNA have to be downloaded from special DB like GeneCodeDB
#grep -hwF "snoRNA" $ENSEMBL_DB/Homo_sapiens.GRCh38.91.sorted.genes_only.gtf >$OUTPUT/Homo_sapiens.GRCh38.91.sorted.snoRNA_only.gtf
#grep -hwF "snRNA" $ENSEMBL_DB/Homo_sapiens.GRCh38.91.sorted.genes_only.gtf >$OUTPUT/Homo_sapiens.GRCh38.91.sorted.snRNA_only.gtf
#grep -hwF "miRNA" $ENSEMBL_DB/Homo_sapiens.GRCh38.91.sorted.genes_only.gtf >$OUTPUT/Homo_sapiens.GRCh38.91.sorted.miRNA_only.gtf
#grep -hwF "snRNA" $ENSEMBL_DB/Homo_sapiens.GRCh38.91.sorted.genes_only.gtf >$OUTPUT/Homo_sapiens.GRCh38.91.sorted.snRNA_only.gtf
## For mRNA we need to get all exons from full annotion file
#grep -hwF "exon" $ENSEMBL_DB/Homo_sapiens.GRCh38.91.sorted.gtf >$OUTPUT/Homo_sapiens.GRCh38.91.sorted.mRNA.exons_only.gtf
#grep -hwF "lincRNA" $ENSEMBL_DB/Homo_sapiens.GRCh38.91.sorted.genes_only.gtf >$OUTPUT/Homo_sapiens.GRCh38.91.sorted.lincRNA_only.gtf
## Antisense RNA are part of exon file
#grep -hwF "antisense_RNA" $OUTPUT/Homo_sapiens.GRCh38.91.sorted.mRNA.exons_only.gtf >$OUTPUT/GRCh38.91.antisense_mRNA.gtf
## mRNA in old annotation was CDS + 3UTR + 5UTR in new everything is in exon biotype
#grep -hwF "CDS" $ENSEMBL_DB/Homo_sapiens.GRCh38.91.sorted.gtf >$OUTPUT/Homo_sapiens.GRCh38.91.sorted.CDS_only.gtf
#grep -hwF "three_prime_utr" $ENSEMBL_DB/Homo_sapiens.GRCh38.91.sorted.gtf >$OUTPUT/Homo_sapiens.GRCh38.91.sorted.three_prime_utr_only.gtf
#grep -hwF "five_prime_utr" $ENSEMBL_DB/Homo_sapiens.GRCh38.91.sorted.gtf >$OUTPUT/Homo_sapiens.GRCh38.91.sorted.five_prime_utr_only.gtf
#grep -hwF "start_codon" $ENSEMBL_DB/Homo_sapiens.GRCh38.91.sorted.gtf >$OUTPUT/Homo_sapiens.GRCh38.91.sorted.start_codon_only.gtf
#grep -hwF "stop_codon" $ENSEMBL_DB/Homo_sapiens.GRCh38.91.sorted.gtf >$OUTPUT/Homo_sapiens.GRCh38.91.sorted.stop_only.gtf
##
# GTF to BED transformation
cd $OUTPUT
for a in GRCh38.91.antisense_mRNA.gtf
do
	FILE=$a
	FILENAME=${a%.*}
	printf "Converting GTF file $FILE to BED format \n" 
	awk -F'[\t|;]' -v OFS="\t" '{for(i=1;i<=NF;i++){if($i~/ gene_name /){print $1, $4, $5, $i, $7}}}' $FILE > $BED_OUTPUT/$FILENAME.bed
#	awk -F'[\t|;]' -v OFS="\t" '{print $1, $4, $5, $9, $7}' $FILE > $BED_OUTPUT/$FILENAME.bed
	printf "Done converting GTF file $FILE to BED format \n"
	printf "Sorting BED file $FILENAME \n"
	$SORTBED -i $BED_OUTPUT/$FILENAME.bed > $BED_OUTPUT/$FILENAME.sorted.bed
	printf "Done sorting BED file $FILENAME \n"
done
#wait
wait
## SED rename of UCSC chromosome names to ENSEMBL
## option -E undocumented for usage of extended regular expression
# sed -E s/^chr// ucsc_hg38.repetitive_elements.sorted.bed >ucsc_hg38.repetitive_elements.ren.sorted.bed
##
## Merge all annotated exone genome elements in order to get intronic regions only
# cd $BED_OUTPUT
# $SORTBED -i GRCh38.91.mRNA_exons.sorted.bed | \
# $MERGEBED -i stdin > GRCh38.91.mRNA_exons.sorted.merged.bed
#
## To define intronic regions, we just need to subtract the exonic region from the genic region.
# bedtools subtract -a GRCh38.91.genes.sorted.bed -b GRCh38.91.mRNA.exons.sorted.bed | \
# gzip > GRCh38.91.mRNA_introns.bed.gz
#
## Define intergenic regions, we use complementBed to find regions not covered by genes.
# date +"%d/%m/%Y %H:%M:%S" 
# printf "Converting whole genome annotation GTF file to BED format \n" 
# gzip -dc $ENSEMBL_DB/Homo_sapiens.GRCh38.91.sorted.gtf.gz \
# | awk -F'[\t|;]' -v OFS="\t" '{for(i=1;i<=NF;i++){if($i~/ gene_name /){print $1, $4, $5, $i, $7}}}'  >$BED_OUTPUT/Homo_sapiens.GRCh38.91.bed
# $SORTBED -i Homo_sapiens.GRCh38.91.bed > Homo_sapiens.GRCh38.91.sorted.bed
# date +"%d/%m/%Y %H:%M:%S" 
# printf "Done Converting whole genome annotation GTF file to BED format \n"

#
# chromosome sizes can be downloaded here (be aware that it is UCSC version where chromosomes are marked as chrXX not XX like in ENSEMBL)
# Also you SORTBED is sorting as follows 1,11,12,13..2,21,22..3,4... so you need to manualy resort the genome chromosome table!
# https://genome.ucsc.edu/goldenpath/help/hg38.chrom.sizes
## Original script:
# mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
#         "select chrom, size from hg19.chromInfo"  > hg19.genome
#  
# zcat gencode.v$v.annotation.gtf.gz |
# awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' |
# bedtools2/bin/sortBed |
# bedtools2/bin/complementBed -i stdin -g hg19.genome |
# gzip > gencode_v19_intergenic.bed.gz
#
#
# $SORTBED -i $BED_OUTPUT/Homo_sapiens.GRCh38.91.sorted.bed |
# $COMPLEMENTBED -i stdin -g $ENSEMBL_DB/hg38_chromosome_sizes_modif.genome |
# gzip > GRCh38.91.intergenic.bed.gz

# $SORTBED -i $BED_OUTPUT/Homo_sapiens.GRCh38.91.bed > $BED_OUTPUT/Homo_sapiens.GRCh38.91.sorted.bed #already done and modified be careful when overwritting
# $COMPLEMENTBED -i $BED_OUTPUT/Homo_sapiens.GRCh38.91.sorted.bed -g $ENSEMBL_DB/hg38_ENSEMBL_primary_ass_chromosome_sizes_modif.genome > GRCh38.91.intergenic.bed
##################################################################################################################################################################











