#!/usr/bin/python

import sys
import itertools

class interval:
  def __init__(self,line):
    self.chrom = line[0]
    self.type = line[1]
    self.feature = line[2]
    self.start = int(line[3])
    self.end = int(line[4])
    self.score = line[5]
    self.strand = line[6]
    self.frame = line[7]
    self.attributes = line[8].split("; ")
    self.attributes = {attribute.split(" ")[0]: attribute.split(" ")[1].strip("\"").strip(";") for attribute in self.attributes}


  def __str__(self):
   return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.chrom,
                                                      self.type,
                                                      self.feature,
                                                      self.start,
                                                      self.end,
                                                      self.score,
                                                      self.strand,
                                                      self.frame,
                                                      self.__printAttrSorted__())


  def __printAttrSorted__(self):
    return "; ".join([" ".join([key,"\"{}\"".format(self.attributes[key])]) for key in ["gene_id",
                                                                                        "transcript_id",
                                                                                        "exon_number",
                                                                                        "intron_number",
                                                                                        "intron_number_end",
                                                                                        "gene_name",
                                                                                        "gene_source",
                                                                                        "gene_biotype",
                                                                                        "transcript_name",
                                                                                        "transcript_source",
                                                                                        "MB_intron_id",
                                                                                        "exon_id"] if key in self.attributes.keys()])


###################################################################################################################################################

def printEnsemblDefault(genes,transcripts,exons):
  for gene in genes:
    print genes[gene]
    for transcript in transcripts[gene]:
      print transcript
      for exon in exons[transcript.attributes["transcript_id"]]:
        print exon

def getLengthsOfConsecutive(s):
  return [len(tuple(g[1])) for g in itertools.groupby([i-j for i,j in enumerate(s)])]     # Calculates lengths of span of consecutive integers

def setToInterval(s):
  s = sorted(list(s))
  lngths = getLengthsOfConsecutive(s)
  cml = 0 
  starts = []
  ends = []
  for l in lngths:
    starts.append(s[cml])
    cml += l 
    ends.append(s[cml-1])
  return [(starts[x],ends[x]) for x in range(len(starts))]


def printIntrons(genes,transcripts,exons):
  i = 1   # Intron counter for MBIntron codes, don't reset
  for gene in genes:
    for transcript in transcripts[gene]:
      transcriptPositions = set([x for x in range(transcript.start,transcript.end + 1)])
                                                                                                      # Exon positions
      exonsPositions = []
      for exon in exons[transcript.attributes["transcript_id"]]:
        exonsPositions += set([x for x in range(exon.start,exon.end + 1)])
      
      intronIntervals = setToInterval(transcriptPositions.difference(exonsPositions))                 # Transforms set of integers into intervals
      j = 1                                                                                           # Intron counter for transcript -> Reset each itteration
      for intronInterval in intronIntervals:
        intronInterval = interval([transcript.chrom,
                                   transcript.type,
                                   transcript.feature,
                                   intronInterval[0],
                                   intronInterval[1],
                                   transcript.score,
                                   transcript.strand,
                                   transcript.frame,
                                   transcript.__printAttrSorted__()
                                   ])
        intronInterval.feature = "intron"
        intronInterval.attributes["MB_intron_id"] = "MBIntron%011d" % i
        intronInterval.attributes["intron_number"] = j
        intronInterval.attributes["intron_number_end"] = -(len(intronIntervals) + 1 - j)
        if intronInterval.strand == "-":    # Flip intron numbering for - strand transcripts
          intronInterval.attributes["intron_number"], intronInterval.attributes["intron_number_end"] = -intronInterval.attributes["intron_number_end"], -intronInterval.attributes["intron_number"]
        j += 1
        i += 1
        print intronInterval



def main():
  genes = {}
  transcripts = {}
  exons = {}
  for line in sys.stdin:# open(sys.argv[1]):
    if line.startswith("#"):
      print line.rstrip()
      continue
    line = line.rstrip().split("\t")
    line = interval(line)
    if line.feature == "gene":
      genes[line.attributes["gene_id"]] = line
      transcripts[line.attributes["gene_id"]] = []
      continue
    if line.feature == "transcript":
      transcripts[line.attributes["gene_id"]].append(line)
      exons[line.attributes["transcript_id"]] = []
      continue
    if line.feature == "exon":
      exons[line.attributes["transcript_id"]].append(line)
      continue
  
  #printEnsemblDefault(genes,transcripts,exons)
  printIntrons(genes,transcripts,exons)



if __name__ == "__main__":
  main()
