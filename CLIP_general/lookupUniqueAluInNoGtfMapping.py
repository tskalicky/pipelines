#!/usr/bin/python

import sys
import subprocess

uniqueAluReads = []
multiMap = {k:0 for k in range(1,21)}
print multiMap

if not sys.stdin.isatty():
  for line in sys.stdin:
    line = line.rstrip().split("\t")
    uniqueAluReads.append(line[0])
    
print str(len(uniqueAluReads)) + " unique reads read from stdin"
 
output = subprocess.Popen(["samtools", "view","/homes2/marek/CLIP/FTO/tophat_noGTF/FTO_CLIP_3.bam"], stdout=subprocess.PIPE).communicate()[0]
print "bam file read OK"

for line in output.split("\n"):
  line = line.rstrip().split("\t")
#  print line[0] in uniqueAluReads
  if line[0] in uniqueAluReads:
    for i in line:
      if "NH:i:" in i:
        multiMap[int(i.split(":")[2])] += 1
        
for key in multiMap:
  print key,multiMap[key]/int(key)