#!/bin/bash
# One liner to rename multiple files from bash
for a in *.fastq.gz; do mv -i "${a}" "${a/_FW/_R1}"; done