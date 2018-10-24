#!/bin/bash
# Sequence library download from NCBI GEO archove using NCBI sra-toolkit
# -e and -p = num of CPUs
fasterq-dump -e 4 -p -x SRR500121 && pigz -v -p 4 SRR500121.fastq