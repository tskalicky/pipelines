#!/usr/bin/python
import sys
import os
from itertools import islice
#
"""
This script will take file containing just sequence names from originally 
unmapped uridylated reads that were TTT trimmed and mapped.
Usage python trimT.py mapped_trimmed_seq_fastq untrimmed_fastq output_filename
"""
#
os.chdir("/home/tomas/ownCloud/CEITEC_lab/Dis3L2/Dasa_spikein/uridylation/sorted_uridyl/fastq/")
#
# seq_mapped = sys.argv[1]
# untrimmed = sys.argv[2]
# output = sys.argv[3]
#
seq_mapped = "pokus_trimmed.fastq"
untrimmed = "pokus.fastq"
output = "pokus_original.fastq"
#
f1 = open(seq_mapped, 'r')
FILENAME = os.path.splitext(os.path.basename(seq_mapped))[0] # To get filename without extension
f2 = open("Seq_Names_%s.txt" % FILENAME, 'w')
f3 = open(untrimmed, "r")
f4 = open(output, "a")
print "Extracting sequence names from {}".format(seq_mapped) # using format standart way in python 2.7
name = ''
seq = ''
qual = ''
ID = ''
for line in f1:
    if line[0] == '@':
        name = line.replace('\n', '')
        f2.write(name + '\n')
f2.close()
#
AI_DICT = {}
f2 = open("Seq_Names_%s.txt" % FILENAME, 'r')
for line in f2:
    AI_DICT[line.replace('\n', '')] = 1
#
skip = 0
for line in f3:
    if line[0] == '@':
        ID = line.replace('\n', '')
        # print ID
        if ID in AI_DICT:
            print "Sequence ",ID," in list"
            f4.write(ID + '\n')
            for _ in range(3):
                f4.write(f3.next())
                skip = 0
        else:
            print "Sequence " ,ID," is not in the list"
            skip = 1
    else:
        if not skip:
            f4.write(ID + '\n')
            for _ in range(3):
                f4.write(f3.next())
                
                    
f1.close()
f2.close()
f3.close()
f4.close()

        
        

