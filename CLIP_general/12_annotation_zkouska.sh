#!/bin/bash
#PBS -l walltime=24:0:0 
#PBS -q default@wagap-pro.cerit-sc.cz 
#PBS -l select=1:ncpus=2:mem=250gb:scratch_local=200gb:os=debian9
#PBS -j oe
#PBS -N sortBED_ABH8
#export PBS_SERVER=wagap-pro.cerit-sc.cz # needed only when executing from arien frontends
#
## initialize the required application
module add bedtools-2.26.0
############################################################################################
### Variables
INPUTDIR="/storage/brno3-cerit/home/tskalicky/ABH8/anotace/ABH8_RNAseq1-3"
# manual merging
cd $INPUTDIR
sortBed -i temp_6_mRNA_exons.ABH8_RNAseq1-3_GRCh38.pcr_dedupl.rg.merged.bed > temp_6_mRNA_exons.ABH8_RNAseq1-3_GRCh38.pcr_dedupl.rg.sorted.bed 
mergeBed -c 10,12 -o distinct,collapse -i temp_6_mRNA_exons.ABH8_RNAseq1-3_GRCh38.pcr_dedupl.rg.sorted.bed > temp_6_mRNA_exons.merged.ABH8_RNAseq1-3_GRCh38.pcr_dedupl.rg.sorted.bed 
#
#sortBed -i temp_6_mRNA_exons.ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.bed > temp_6_mRNA_exons.ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.sorted.bed 
#mergeBed -c 10,12 -o distinct,collapse -i temp_6_mRNA_exons.ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.sorted.bed > temp_6_mRNA_exons.merged.ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.sorted.bed
#
#sortBed -i temp_2_tRNA.METTL16PAR1-2_GRCh38.pcr_dedupl.rg.bed > temp_2_tRNA.METTL16PAR1-2_GRCh38.pcr_dedupl.rg.sorted.bed 
#mergeBed -c 10,12 -o distinct,collapse -i temp_2_tRNA.METTL16PAR1-2_GRCh38.pcr_dedupl.rg.sorted.bed > temp_2_tRNA.merged.METTL16PAR1-2_GRCh38.pcr_dedupl.rg.sorted.bed 
#
#sortBed -i temp_2_tRNA.METTL16UV1-2_GRCh38.pcr_dedupl.rg.bed > temp_2_tRNA.METTL16UV1-2_GRCh38.pcr_dedupl.rg.sorted.bed 
#mergeBed -c 10,12 -o distinct,collapse -i temp_2_tRNA.METTL16UV1-2_GRCh38.pcr_dedupl.rg.sorted.bed > temp_2_tRNA.merged.METTL16UV1-2_GRCh38.pcr_dedupl.rg.sorted.bed
#
sortBed -i temp_2_tRNA.ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.bed > temp_2_tRNA.ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.sorted.bed
mergeBed -c 10,12 -o distinct,collapse -i temp_2_tRNA.ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.sorted.bed > temp_2_tRNA.mrgd.ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.sorted.bed

# Manual count one-liner
for b in temp*.bed; do wc -l $b; done  | sort -s -k1,1 -rn >> Summary_ABH8-Input2.txt