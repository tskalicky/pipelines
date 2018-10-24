#!/bin/bash
##
# script to copy multiple files from different location on remote server using sftp
# for batch script use this command:
# sftp tomas@147.231.253.13 -b batchFile.txt
##
OUTPUT_DIR="/home/users/tskalicky/Hinxton/trimming"
OUTPUT_DIR2="/home/users/tskalicky/Hinxton/PacBio"
WALL_MISEQ1_PE1="/media/4TB3/Tomas_data/Hinxton/Trimming/18098_1_6_1_trimmed_Paired.fq.gz"
WALL_MISEQ1_PE2="/media/4TB3/Tomas_data/Hinxton/Trimming/18098_1_6_2_Paired.fq.gz"
WALL_MISEQ1_UP1="/media/4TB3/Tomas_data/Hinxton/Trimming/18098_1_6_1_trimmed_Unpaired.fq.gz"
WALL_MISEQ1_UP2="/media/4TB3/Tomas_data/Hinxton/Trimming/18098_1_6_2_Unpaired.fq.gz"
WALL_MISEQ1_UP3="/media/4TB3/Tomas_data/Hinxton/Trimming/18098_1_6_both_Unpaired.fq.gz"
#
WALL_MISEQ2_PE1="/media/4TB3/Tomas_data/Hinxton/Trimming/18021_1_6_1_trimmed_Paired.fq.gz"
WALL_MISEQ2_PE2="/media/4TB3/Tomas_data/Hinxton/Trimming/18021_1_6_2_Paired.fq.gz"
WALL_MISEQ2_UP1="/media/4TB3/Tomas_data/Hinxton/Trimming/18021_1_6_1_trimmed_Unpaired.fq.gz"
WALL_MISEQ2_UP2="/media/4TB3/Tomas_data/Hinxton/Trimming/18021_1_6_2_Unpaired.fq.gz"
WALL_MISEQ2_UP3="/media/4TB3/Tomas_data/Hinxton/Trimming/18021_1_6_both_Unpaired.fq"
#
WALL_HISEQ1_1="/media/4TB3/Tomas_data/Hinxton/Trimming/Mbr04_Wallacemonas_both_6_1_trimmed_Paired.fq.gz"
WALL_HISEQ1_2="/media/4TB3/Tomas_data/Hinxton/Trimming/Mbr04_Wallacemonas_both_6_2_trimmed_Paired.fq.gz"
WALL_HISEQ1_UP1="/media/4TB3/Tomas_data/Hinxton/Trimming/Mbr04_Wallacemonas_both_6_1_trimmed_Unpaired.fq.gz"
WALL_HISEQ1_UP2="/media/4TB3/Tomas_data/Hinxton/Trimming/Mbr04_Wallacemonas_both_6_2_trimmed_Unpaired.fq.gz"
WALL_HISEQ1_UP3="/media/4TB3/Tomas_data/Hinxton/Trimming/Mbr04_Wallacemonas_both_6_both_trimmed_Unpaired.fq.gz"
#
WALL_PACBIO="/media/4TB3/Tomas_data/Hinxton/Assembly/PacBio/HGAP/Mbr04_Wallacemonas/Mbr04_filtered_longreads.fasta"
TRUSTED_CONTIGS="/media/4TB3/Tomas_data/Hinxton/Assembly/PacBio/HGAP/Mbr04_Wallacemonas/HGAP_25/polished_assembly.fasta"
#
get $WALL_MISEQ1_PE1 $OUTPUT_DIR
get $WALL_MISEQ1_PE2 $OUTPUT_DIR
get $WALL_MISEQ1_UP1 $OUTPUT_DIR
get $WALL_MISEQ1_UP2 $OUTPUT_DIR
get $WALL_MISEQ1_UP3 $OUTPUT_DIR
get $WALL_MISEQ2_PE1 $OUTPUT_DIR
get $WALL_MISEQ2_PE2 $OUTPUT_DIR
get $WALL_MISEQ2_UP1 $OUTPUT_DIR
get $WALL_MISEQ2_UP2 $OUTPUT_DIR
get $WALL_MISEQ2_UP3 $OUTPUT_DIR
get $WALL_HISEQ1_1 $OUTPUT_DIR
get $WALL_HISEQ1_2 $OUTPUT_DIR
get $WALL_HISEQ1_UP1 $OUTPUT_DIR
get $WALL_HISEQ1_UP2 $OUTPUT_DIR
get $WALL_HISEQ1_UP3 $OUTPUT_DIR
get $WALL_PACBIO $OUTPUT_DIR2
get $TRUSTED_CONTIGS $OUTPUT_DIR2
EOF
#for command line only:
# sftp tomas@147.231.253.13 << EOF
#   get $WALL_MISEQ1_PE1 $OUTPUT_DIR
#   get $WALL_MISEQ1_PE2 $OUTPUT_DIR
#   get $WALL_MISEQ1_UP1 $OUTPUT_DIR
#   get $WALL_MISEQ1_UP2 $OUTPUT_DIR
#   get $WALL_MISEQ1_UP3 $OUTPUT_DIR
#   get $WALL_MISEQ2_PE1 $OUTPUT_DIR
#   get $WALL_MISEQ2_PE2 $OUTPUT_DIR
#   get $WALL_MISEQ2_UP1 $OUTPUT_DIR
#   get $WALL_MISEQ2_UP2 $OUTPUT_DIR
#   get $WALL_MISEQ2_UP3 $OUTPUT_DIR
#   get $WALL_HISEQ1_1 $OUTPUT_DIR
#   get $WALL_HISEQ1_2 $OUTPUT_DIR
#   get $WALL_HISEQ1_UP1 $OUTPUT_DIR
#   get $WALL_HISEQ1_UP2 $OUTPUT_DIR
#   get $WALL_HISEQ1_UP3 $OUTPUT_DIR
#   get $WALL_PACBIO $OUTPUT_DIR2
#   get $TRUSTED_CONTIGS $OUTPUT_DIR2
# EOF