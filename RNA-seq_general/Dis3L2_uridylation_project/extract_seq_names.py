#!/usr/bin/python
import sys
import os
#
"""Usage python trimT.py seq_names_fastq"""
#
seq_names = sys.argv[1]
FILENAME = os.path.splitext(os.path.basename(seq_names))[0] # To get filename without extension
print "Extracting sequence names from {}".format(seq_names) # using format standart way in python 2.7
#
name = ''
seq = ''
qual = ''
n = 1
with open('Seq_Names_' + seq_names, 'w') as res:
    for line in open(seq_names):
        if n == 1:
            name = line.replace('\n', '')
            res.write(seq + '\n')
            n += 1
        elif n == 2:
            seq = line.replace('\n', '')
            n += 1
        elif n == 3:
            n += 1
        elif n == 4:
            qual = line.replace('\n', '')