#!/usr/bin/python
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import itertools
import os
os.chdir("/home/tomas/CEITEC_lab/Dis3L2/Dasa_spikein/data/trimmed")

#Setup variables (could parse command line args instead)
FILE_FW = "DIS3l2_OAT_cyto_flexbar_1.fastq"
FILE_RV = "DIS3l2_OAT_cyto_flexbar_1.fastq"
FILE_OUT = "DIS3l2_OAT_cyto_interleaved.fastq"
#
handle = open(FILE_OUT, "w")
count = 0
#
f_iter = FastqGeneralIterator(open(FILE_FW,"rU"))
r_iter = FastqGeneralIterator(open(FILE_RV,"rU"))
for (f_id, f_seq, f_q), (r_id, r_seq, r_q) \
in itertools.izip(f_iter,r_iter):
    assert f_id == r_id
    count += 2
    #Write out both reads with "/1" and "/2" suffix on ID
    handle.write("@%s/1\n%s\n+\n%s\n@%s/2\n%s\n+\n%s\n" 
                 % (f_id, f_seq, f_q, r_id, r_seq, r_q))
handle.close()
print "%i records written to %s" % (count, FILE_OUT)

# correction of interleaved file with bash script 
# I don|t know how to do it on fly with previous python script
# sed -e 's/\/1\/1/\/1/g' -e 's/\/1\/2/\/2/g' DIS3l2_OAT_cyto_interleaved.fastq >DIS3l2_OAT_cyto_interleaved2.fastq
# sed -e 's/\/1\/1/\/1/g' -e 's/\/1\/2/\/2/g' DIS3l2_OAT_nucl_interleaved.fastq >DIS3l2_OAT_nucl_interleaved2.fastq
# sed -e 's/\/1\/1/\/1/g' -e 's/\/1\/2/\/2/g' DIS3l2_OAT_nonfrac_interleaved.fastq >DIS3l2_OAT_nonfrac_interleaved2.fastq
#sed -e 's/\/1\/1/\/1/g' -e 's/\/1\/2/\/2/g' DIS3l2_U12_nonfrac_interleaved.fastq >DIS3l2_U12_nonfrac_interleaved2.fastq
# sed -e 's/\/1\/1/\/1/g' -e 's/\/1\/2/\/2/g' DIS3l2_U12_nucl_interleaved.fastq >DIS3l2_U12_nucl_interleaved2.fastq
# sed -e 's/\/1\/1/\/1/g' -e 's/\/1\/2/\/2/g' DIS3l2_U12_cyto_interleaved.fastq >DIS3l2_U12_cyto_interleaved2.fastq