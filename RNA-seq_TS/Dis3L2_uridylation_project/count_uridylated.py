#!/usr/bin/python
import sys
import os
from collections import defaultdict
from collections import OrderedDict
import csv
import pandas as pd

"""Usage python count_uridylated.py fasta_file
This script will take file containing sequences in fasta format from originally 
urydylated sequences and outputs uridylated lengths, TTT seqs and number of their occurences in dataset
Works only with Python v. 2.7.15+
"""
#
os.chdir("/home/tomas/ownCloud/CEITEC_lab/Dis3L2/Dasa_spikein/uridylation/NEW/sorted_original/fastq/")
file = "pokus.fasta"
# file = sys.argv[1]
FILENAME = os.path.splitext(os.path.basename(file))[0] # To get filename without extension
print "Counting uridylated sequences from {}".format(file) # using format standart way in python 2.7
# print "Counting uridylated sequences from", file
FULLY_Uridyl = open("Fully_Uridyl_%s.fasta" % FILENAME, 'w')
# Define function
def countT(seq):
    length = -4
    if seq[-4:].count('T') == 3:
        try:
            if seq[length - 1] == 'T':
                while seq[length - 1] == 'T':
                    length = length - 1
                    return(length)
            else:
                return(length)
        except IndexError:
            return "jahody"                        
    elif seq[-4:].count('T') == 4:
        try:
            while seq[length - 1] == 'T':
                length = length -1
                #return(length)
            else:
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
# Commands            
name = ''
seq = ''
TTT_seq = ''
n = 1
vyskyty = defaultdict(int)
vyskyty.default_factory


for line in open(file):
    if n == 1:
        name = line.replace('\n', '')
        n += 1
    elif n == 2:
        seq = line.replace('\n', '')
        length = countT(seq)
        if length == "jahody":
            # print "There are fully uridylated sequences:"
            FULLY_Uridyl.write(name + '\n')
            FULLY_Uridyl.write(seq + '\n')
            print(name)
            print(seq)
            # print(len(seq))
            TTT_seq = seq
            vyskyty[TTT_seq] += 1
            n = 1
        else:
            if length:
                TTT_seq = seq[length:]
                vyskyty[TTT_seq] += 1
            n = 1
d_sorted_by_key = OrderedDict(sorted(vyskyty.items(), key=lambda x: (len(x[0])), reverse=False))  # k = key, v = value, sort according to key length only
# d_sorted_by_key = OrderedDict(sorted(vyskyty.items(), key=lambda x: (x[1],len(x[0])), reverse=False)) # k = key, v = value, sort according to key 1. alphabetically and then according key length
### debug
print vyskyty.items()
# print sorted(vyskyty.items())
# keys = vyskyty.keys()
# values = vyskyty.values()
# print keys
# print values
### end debug
#
#writing CSV
with open('CountT_%s.csv' % FILENAME, 'wb') as f: # binary is better, avoids blank lines in some python 2 versions
    writer = csv.writer(f, delimiter="\t")
    keys=["uridylation","count", "length"] # define multiple keys if each key has multiple
    writer.writerow(keys) # write extra row with key names, add [""]+keys for extra empty first field for key name insertion
    for (k, v) in d_sorted_by_key.items(): # k = key, v = value
        writer.writerow([k] + [v] + [len(k)]) # k and v need to be in [] because we want them to input as they are
            
##reading CSV
#with open('CountT_' + file, 'r') as f:
#    reader = csv.reader(f)    
#    mydict = vyskyty.defaultdict(list)
#    for row in reader:
#        mydict[row[0]] = row[1:]
##            
# Read in data and examine first 10 rows
TTT_distrib = pd.read_csv('CountT_%s.csv' % FILENAME, delimiter="\t")
TTT_distrib.head(10)




