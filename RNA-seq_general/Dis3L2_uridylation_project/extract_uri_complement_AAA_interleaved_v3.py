#!/usr/bin/python
# from Bio import SeqIO
import sys
import os
from itertools import islice

"""Usage python extract_uri_complement_AAA.py fastq_file"""
## Description ##
# This script takes fastq file and outputs just sequences in fasta sorted to uridylated and nonuridylated files
#####
# file = sys.argv[1]
file = "pokus.fastq"
FILENAME = os.path.splitext(os.path.basename(file))[0] # To get filename without extension
print "Extracting uridylated sequences from {}".format(file) # using format standart way in python 2.7
# print "Extracting uridylated sequences from", file
os.chdir("/home/tomas/ownCloud/CEITEC_lab/Dis3L2/Dasa_spikein/uridylation/NEW/sorted_original/fastq")

def extractA(seq):
    length = -4
    if seq[-4:].count('A') == 3:
        try:
            if seq[length - 1] == 'A':
                while seq[length - 1] == 'A':
                    length = length - 1
                    return(length)
            else:
                return(length)
        except IndexError:
            return "jahody"                        
    elif seq[-4:].count('A') == 4:
        try:
            while seq[length - 1] == 'A':
                length = length -1
                #return(length)
            else:
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
n = 1
N = 8
#
COMPLEMENT = open('Compl_' + file, "w")
NON = open('NOT_compl_' + file, "w")
#
with open(file, "r") as interleaved:
    while True:
        next_n_lines = list(islice(interleaved, N))
        if not next_n_lines:
            break # will terminate loop before conditions are applied because there are no more sequences in file
        if str(next_n_lines[1]).startswith('GGCGTCACTGTTGCGCTTCATAGAC'): # OAT primer
            NON.writelines( "%s" % item for item in next_n_lines)
        elif str(next_n_lines[1]).startswith('ACCTTATTCACGCCTAAAAA'): # RNU12 primer
            NON.writelines( "%s" % item for item in next_n_lines)
        else:
            seq = str(next_n_lines[5])
            name = str(next_n_lines[4])
            length = extractA(seq)
            if length == "jahody":
                print "There are fully AAA sequences:"
                print(name)
                print(seq)
                COMPLEMENT.writelines( "%s" % item for item in next_n_lines)
            else:
                if length:
                    COMPLEMENT.writelines( "%s" % item for item in next_n_lines)
                else:
                    NON.writelines( "%s" % item for item in next_n_lines)

COMPLEMENT.close()
NON.close()




