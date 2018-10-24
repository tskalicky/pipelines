#!/usr/bin/python
from Bio import SeqIO
import sys

"""Usage python trimT.py fastq_file"""



file = sys.argv[1]


def trimT(seq, qual):
    length = -4
    if seq[-4:].count('T') == 3:
        if seq[-5] != 'T':
            return -3
        else:
            while seq[length - 1] == 'T':
                length = length - 1
            else:
                return(length)
    elif seq[-4:].count('T') == 4:
        while seq[length - 1] == 'T':
            length = length - 1
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
with open('trimT_' + file, 'w') as res:
    for line in open(file):
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
            length = trimT(seq, qual)
            if length == "jahody":
                print(name)
                print(seq)
                n = 1
            else:
                res.write(name +'\n')
                if length:
                    res.write(seq[:length] + '\n')
                    res.write('+\n')
                    res.write(qual[:length] + '\n')
                else:
                    res.write(seq + '\n')
                    res.write('+\n')
                    res.write(qual + '\n')
                n = 1




