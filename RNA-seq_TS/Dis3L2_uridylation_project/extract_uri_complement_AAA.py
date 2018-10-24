#!/usr/bin/python
# from Bio import SeqIO
import sys
import os
import io
import re

"""Usage python extract_uri_complement_AAA.py fastq_file_FW fastq_file_RV"""
## Description ##
# This script takes fastq file and outputs just sequences in fasta sorted to uridylated and nonuridylated files
#####
#file_FW = "pokus.fastq"
file_FW = sys.argv[1]
file_RV = sys.argv[2]
FILENAME = os.path.splitext(os.path.basename(file_FW))[0] # To get filename without extension
print "Extracting uridylated sequences from {}".format(file_FW) # using format standart way in python 2.7
# print "Extracting uridylated sequences from", file_FW
# os.chdir("/home/tomas/ownCloud/CEITEC_lab/Dis3L2/Dasa_spikein/uridylation/NEW/sorted_original/fastq")

def extractT(seq, qual):
    length = -4
    if line.startswith('GGCGTCACTGTTGCGCTTCATAGAC'): # OAT primer
            return(length)
        elif line.startswith('ACCTTATTCACGCCTAAAAA'): # RNU12 primer
            return(length)
        elif seq[-4:].count('A') == 3:
                return(length)
            elif seq[-4:].count('A') == 4:
                while seq[length - 1] == 'A':
                    return(length)
                else:
                    try:
                        if seq[length - 2] == 'A':
                            length = length - 1
                            while seq[length - 1] == 'A':
                                length = length - 1
                            else:
                                return(length)
                        else:
                            return(length)
                    except IndexError:
                        return "jahody" 
            
name = ''
seq = ''
qual = ''
n = 1
#
URIDYL = open('Uridyl_compl_' + file_FW, "a")
NON = open('NOT_uridyl_compl_' + file_FW, "a")
#
for line in open(file_FW, "r"):
    if n == 1:
        name = line.replace('\n', '')
        n += 1
    elif n == 2:
        seq = line.replace('\n', '')
        n += 1
    elif n == 3:
        n += 1
    elif n == 4:
        qual = line.replace('\n', '')
        length = extractT(seq, qual)
        if length == "jahody":
            print "There are fully uridylated sequences:"
            print(name)
            print(seq)
            URIDYL.write(name + '\n')
            URIDYL.write(seq + '\n')
            URIDYL.write('+\n')
            URIDYL.write(qual + '\n')
            n = 1
        else:
            if length:
                URIDYL.write(name + '\n')
                URIDYL.write(seq + '\n')
                URIDYL.write('+\n')
                URIDYL.write(qual + '\n')
            else:
                NON.write(name + '\n')
                NON.write(seq + '\n')
                NON.write('+\n')
                NON.write(qual + '\n')
            n = 1
    elif n == 5:
        n += 1
    elif n == 6:
        n += 1
    elif n == 7:
        n += 1
    elif n == 8:
        n = 1

URIDYL.close()
NON.close()




