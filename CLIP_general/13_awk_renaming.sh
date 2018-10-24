#!/bin/bash
#mkdir /home/tomas/CEITEC_lab/ABH8/anotace/job_774301_repaired_beds/ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.bam.sorted/renamed_merged
REN="/home/tomas/CEITEC_lab/ABH8/anotace/job_774301_repaired_beds/ABH8_CLIP1-3_GRCh38.pcr_dedupl.rg.merged.bam.sorted/renamed_merged"
for a in *.bed
do
    FILE=$a
    FILENAME=${a%.*.*.*}
    date +"%d/%m/%Y %H:%M:%S"
    echo "Renaming reads in file $FILE"
    # Strand information in bed format has to be in 6th collumn!!!
    awk -F'[\t]' -v OFS='\t' '{print $1, $2, $3, $16, $6, $5}' $FILE > "$REN"/"$FILENAME".ren.bed
    date +"%d/%m/%Y %H:%M:%S"
    echo "Finnished renaming reads in file $FILE"
    #
    date +"%d/%m/%Y %H:%M:%S"
    echo "Sorting and merging reads in file $FILE"
    sortBed -i "$REN"/"$FILENAME".ren.bed | \
    mergeBed -s -c 4,5,6 -o distinct,distinct,distinct -i stdin > "$REN"/"$FILENAME".ren.merged.sorted.bed
    date +"%d/%m/%Y %H:%M:%S"
    echo "Finnished sorting and merging reads in file $FILE"
done
wait
## One liner to repair strand info position
# for a in *.ren.sorted.bed; do awk -F'[\t]' -v OFS='\t' '{print $1, $2, $3, $4, $6, $5}' $a | mergeBed -s -c 4,5,6 -o distinct,distinct,distinct -i stdin > $a".merged"; done
# for a in *.ren.sorted.bed; do awk -F'[\t]' -v OFS='\t' '{print $1, $2, $3, $4, $6, $5}' $a >$a".ren2" ; done
# for a in *.ren.sorted.bed; do mergeBed -s -c 4,5,6 -o distinct,median,distinct -i $a > $a".merged"; done
#
# local sort merge
# #!/bin/bash
# File="temp_2_tRNA.ABH8_CLIP1-3.ren.bed"
# Name="tRNA.ABH8_CLIP1-3"
# File2="temp_2_tRNA.ABH8_CLIP_vs_RNAseq_bamCompare_bin10_SES_log2_RPKM_smooth.bed"
# sortBed -i $File > $Name.sorted.bed
# mergeBed -s -c 4,5,6 -o distinct,distinct,distinct -i $Name.sorted.bed > $Name.sorted.merged.bed
