#!/bin/bash
#delete unwanted parts of seq. names, complete bodo name, delete extra Tgrayi seq
cd "/home/tomas/CEITEC_lab/Dis3L2/Dasa_spikein/data/trimmed"
awk '{gsub(/\/1\/1/, "\/1"); gsub(/\/1\/2/, "\/2");}' DIS3l2_OAT_cyto_interleaved.fastq
sed -e 's/\/1\/1/\/1/g' -e 's/\/1\/2/\/2/g' DIS3l2_OAT_cyto_interleaved.fastq >DIS3l2_OAT_cyto_interleaved2.fastq