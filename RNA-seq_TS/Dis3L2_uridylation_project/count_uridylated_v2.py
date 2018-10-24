#!/usr/bin/python
# Import libraries
import sys
import os
from collections import defaultdict
from collections import OrderedDict
import csv
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
#
"""Usage python count_uridylated.py fasta_file
This script will take file containing sequences in fasta format from originally 
urydylated sequences and outputs uridylated lengths, TTT seqs and number of their occurences in dataset
Works only with Python v. 2.7.15+
"""
#
# file = "Uridyl_DIS3l2_U12_nonfrac_aligned_unmapped.fa"
file = sys.argv[1]
# OUT_DIR = os.path.dirname(sys.argv[1])
# file = os.path.basename(sys.argv[1])
# os.chdir("/home/tomas/ownCloud/CEITEC_lab/Dis3L2/Dasa_spikein/uridylation/NEW/sorted_original/fastq/")
# os.chdir(OUT_DIR)
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
# print vyskyty.items()
# print sorted(vyskyty.items())
# keys = vyskyty.keys()
# values = vyskyty.values()
# print keys
# print values
### end debug
uridylation = defaultdict(int)
uridylation.default_factory
#
## Creating dictionary and dataframe for histogram
with open('Raw_%s.csv' % FILENAME, 'wb') as f2: # binary is better, avoids blank lines in some python 2 versions
    writer = csv.writer(f2, delimiter="\t")
    keys=["uridylation", "uridyl_length"] # define multiple keys if each key has multiple
    writer.writerow(keys) # write extra row with key names, add [""]+keys for extra empty first field for key name insertion
    #dict = {}
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
                key, value = TTT_seq, len(TTT_seq)
                writer.writerow([key] + [value])
                n = 1
            else:
                if length:
                    TTT_seq = seq[length:]
                    vyskyty[TTT_seq] += 1
                    key, value = TTT_seq, len(TTT_seq)
                    writer.writerow([key] + [value])
                n = 1
    #print(dict)
#
#writing CSV
with open('CountT_%s.csv' % FILENAME, 'wb') as f: # binary is better, avoids blank lines in some python 2 versions
    writer = csv.writer(f, delimiter="\t")
    # keys=["uridylation","count"]
    keys=["uridylation","count", "uridyl_length"] # define multiple keys if each key has multiple
    writer.writerow(keys) # write extra row with key names, add [""]+keys for extra empty first field for key name insertion
    for (k, v) in d_sorted_by_key.items(): # k = key, v = value
        writer.writerow([k] + [v] + [len(k)]) # k and v need to be in [] because we want them to input as they are
        # writer.writerow([k] + [v])

#
##reading CSV
#with open('CountT_' + file, 'r') as f:
#    reader = csv.reader(f)    
#    mydict = vyskyty.defaultdict(list)
#    for row in reader:
#        mydict[row[0]] = row[1:]
##            
# Read in csv data and examine first 10 rows
#TTT_distrib = pd.read_csv('Raw_%s.csv' % FILENAME, delimiter="\t")
#TTT_distrib.head(10)
#print("Plotting Histogram from results stored in file {}").format('Raw_%s.csv' % FILENAME)
#
## matplotlib histogram
# plt.hist(TTT_distrib['uridyl_length'], color = 'blue', edgecolor = 'black',
#          bins = int(180/5))
# Seaborn histogram
#sns.distplot(TTT_distrib['uridyl_length'], hist=True, kde=False, 
#             bins=int(100), color = 'blue',
#             hist_kws={'edgecolor':'black'})
## Density Plot and Histogram
#sns.distplot(TTT_distrib['uridyl_length'], hist=True, kde=False, 
#             bins=int(100), color = 'darkblue', 
#             hist_kws={'edgecolor':'black'},
#             kde_kws={'linewidth': 4})
## Shaded Density plot
#sns.distplot(subset['arr_delay'], hist = False, kde = True,
#                 kde_kws = {'shade': True, 'linewidth': 3}, 
#                  label = airline)
## Density Plot with Rug Plot
## WARNING! Rugplot will get stuck when too many items to plot and will use huge amount of memory
#sns.distplot(TTT_distrib['uridyl_length'], kde = True, rug = True,
#             color = 'darkblue', 
#             kde_kws={'linewidth': 3, "label": "DIS3l2_U12_nonfrac"},
#             rug_kws={'color': 'black'})
# Add labels
#plt.title('Histogram of Uridylation lengths', size = 18)
#plt.xlabel('uridylation lengths', size = 16)
#plt.ylabel('uridylation counts', size = 16)




