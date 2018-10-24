#!/bin/bash
# Usage: converts NCBI Genbank flat file to GFF3 format. bp_genbank2gff3.pl [options] filename(s)
# bp_genbank2gff3.pl is part of a Bioperl package that needs to be installed on your system!

############################################################################################
## Variables ##
export PATH="/home/tomas/anaconda2/bin:$PATH"
echo "PATH is:"
echo $PATH | tr ":" "\n" | nl
# MY_GBK="/home/tomas/ownCloud/CEITEC_lab/genomes/human/ensembl91/STAR_index/OAT_124310076-124444951.flat.gbk"
MY_GBK="/home/tomas/ownCloud/CEITEC_lab/genomes/human/ensembl91/STAR_index/RNU12_42615088-42616003.flat.gbk"
OUTPUT_DIR="/home/tomas/ownCloud/CEITEC_lab/genomes/human/ensembl91/STAR_index"
# Check the tools
GBK2GFF=$(which bp_genbank2gff3.pl)
which $GBK2GFF
### Commands ###
date +"%d/%m/%Y %H:%M:%S $HOSTNAME"
echo "Start converting Genbank flat file to GFF3"
$GBK2GFF --split --outdir $OUTPUT_DIR $MY_GBK
