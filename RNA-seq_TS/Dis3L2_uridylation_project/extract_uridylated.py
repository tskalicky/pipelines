#!/usr/bin/python
# from Bio import SeqIO
import sys
import os

"""Usage python extract_uridylated.py fastq_file"""
## Description ##
# This script takes fastq file and outputs just sequences in fasta sorted to uridylated and nonuridylated files
#####
# file = "pokus.fastq"
file = sys.argv[1]
FILENAME = os.path.splitext(os.path.basename(file))[0] # To get filename without extension
print "Extracting uridylated sequences from {}".format(file) # using format standart way in python 2.7
# print "Extracting uridylated sequences from", file
# os.chdir("/home/tomas/ownCloud/CEITEC_lab/Dis3L2/Dasa_spikein/uridylation/extract_uridylated/")

def extractT(seq, qual):
    length = -4
    if seq[-4:].count('T') == 3:
        return(length)
    elif seq[-4:].count('T') == 4:
        while seq[length - 1] == 'T':
            return(length)
        else:
            try:
                if seq[length - 2] == 'T':
                    length = length - 1
                    while seq[length - 1] == 'T':
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
URIDYL = open('Uridyl_' + file, "a")
NON = open('NOT_uridyl_' + file, "a")
#
for line in open(file, "r"):
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
URIDYL.close()
NON.close()




