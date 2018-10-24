#!/usr/bin/python

import pybedtools as b
import numpy
from scipy import stats
import sys

CLIP1 = b.BedTool("../FTO_1_tophat_noGTF_combined_sorted_multi_random_selected.gff")
CLIP2 = b.BedTool("../FTO_2_tophat_noGTF_combined_sorted_multi_random_selected.gff")
CLIP3 = b.BedTool("../FTO_3_tophat_noGTF_combined_sorted_multi_random_selected.gff")
RNAseq1 = b.BedTool("../FTO_1_RNAseq_tophat_noGTF_combined_sorted_multi_random_selected.gff")
RNAseq2 = b.BedTool("../FTO_2_RNAseq_tophat_noGTF_combined_sorted_multi_random_selected.gff")
RNAseq3 = b.BedTool("../FTO_3_RNAseq_tophat_noGTF_combined_sorted_multi_random_selected.gff")

files = [CLIP1,CLIP2,CLIP3,RNAseq1,RNAseq2,RNAseq3]
names = ["CLIP1","CLIP2","CLIP3","RNAseq1","RNAseq2","RNAseq3"]

CLIP_experiments = [0,1,2]
RNAseq_experiments = [3,4,5]

min_RNAseq = 0.1
DECOMPRESS = False
THRESHOLD = 2
min_reads = 20

b.helpers.set_bedtools_path("/usr/local/bedtools/bedtools2.16.2/bin")

def check_format_gff(line):
  chrom = line[0]
  start = int(line[1])
  stop = int(line[2])
  strand = line[4]
  if "chr" in chrom and \
     start < stop and \
     start != 0 and \
     (strand == "+" or strand =="-"):
    return True
  return False


def create_RPKM_table():
  #merged = b.BedTool("/homes2/marek/annotation/genome_50bp_window.gff")
  merged = b.BedTool("/homes2/marek/annotation/genome_50bp_window.gff")
  result = []
  clusters = []
  res_RPKM = []
  res_raw_reads = []
  for file in files:
    reads_count_milion = float(len(file)) / 1000000
    RPKM = {}
    raw_reads = {}
    for line in merged.intersect(file,loj=True,s=True):  
      cluster_id = str(line).split()[8]
      read_id = str(line).split()[17]
      length_feature_kb = (int(str(line).split()[4]) - int(str(line).split()[3])) / 1000.0

      if read_id == ".":
        coverage_decompressed = 0
      else:
        if DECOMPRESS:
          coverage_decompressed = int(read_id[3:-1].split("_")[1])
        else:
          coverage_decompressed = 1
      if coverage_decompressed == 0:
        continue
      if cluster_id in RPKM:
        RPKM[cluster_id] += coverage_decompressed/(reads_count_milion * length_feature_kb)
        raw_reads[cluster_id] += coverage_decompressed
      else:
        RPKM[cluster_id] = coverage_decompressed/(reads_count_milion * length_feature_kb)
        raw_reads[cluster_id] = coverage_decompressed
        clusters.append(cluster_id)
    
    res_RPKM.append(RPKM)
    res_raw_reads.append(raw_reads)
    del RPKM,raw_reads

    
    
    # Report to file:
  with open("RPKM_clusters.table","w") as f:
    # Print header
    f.write ("{}\t".format("#cluster_name"))
    for i in names:
      f.write("{}\t".format(i))
    f.write("\n")   
    # report
    for cluster_id in clusters:
      RPKM_cluster = []
      raw_read_cluster = []
      for RPKM_table,raw_reads_table in zip(res_RPKM,res_raw_reads):
        try:
          RPKM_cluster.append(RPKM_table[cluster_id])
          raw_read_cluster.append(raw_reads_table[cluster_id])
        except KeyError:
          RPKM_cluster.append(0)
          raw_read_cluster.append(0)
      if sum([x for x in raw_read_cluster[min(CLIP_experiments):max(CLIP_experiments)+1]]) < min_reads:
        continue
      f.write("{}\t".format(cluster_id))
      f.write("{}\t".format("\t".join([str(x) for x in RPKM_cluster])))
      f.write("{}\n".format(str(sum([x for x in raw_read_cluster[min(CLIP_experiments):(max(CLIP_experiments)+1)]]) ))) 
    
def calculate_enrichment():
  numpy.seterr(all="ignore")
  with open("RPKM_enrichment.table","w") as f:
    for line in open("RPKM_clusters.table","r"):
      if line.startswith("#"):
        f.write ("{}\t".format("#Cluster_name"))
        for i in CLIP_experiments:
          f.write("{}\t".format(line.split()[i+1]))
        f.write("t_test\traw_reads\n")
      else:
        RNAseqRPKMs = [float(x) for x in line.split()[min(RNAseq_experiments)+1:(max(RNAseq_experiments)+2)]]
        ClipRPKMs = [float(x) for x in line.split()[min(CLIP_experiments)+1:(max(CLIP_experiments)+2)]]
        RNAseq_average = numpy.mean(RNAseqRPKMs)
        if RNAseq_average > min_RNAseq:
          f.write ("{}\t".format(line.split()[0]))
          for i in CLIP_experiments:
            f.write("{}\t".format(numpy.log2(numpy.divide(float(line.split()[i+1]),float(RNAseq_average)))))
          f.write("{}\t".format(stats.ttest_ind(ClipRPKMs,RNAseqRPKMs)[1]))
          f.write("{}\n".format(line.split()[-1]))



def filter(enrichment):
  if enrichment == "inf":
    return True
  elif enrichment == "-inf":
    return False
  elif enrichment == "nan":
    return False
  elif float(enrichment) > THRESHOLD:
    return True

def filter_ttest(ttest_tuple):
  if ttest_tuple[1] < 0.05:
    return True
  return False

def filter_enrichment():
  enrichments = {}
  for line in open("RPKM_enrichment.table","r"):
    if not line.startswith("#"):
      line = line.rstrip()
      cluster_name = line.split()[0]
      CLIP1 = line.split("\t")[1]
      CLIP2 = line.split("\t")[2]
      CLIP3 = line.split("\t")[3]
      ttest = line.split("\t")[4]
      read_count = line.split("\t")[5]
      
      cnt = 0
      if filter(CLIP1):
        cnt += 1
      if filter(CLIP2):
        cnt += 1
      if filter(CLIP3):
        cnt += 1
      if cnt > 1:
        enrichments[cluster_name] = CLIP1+";"+CLIP2+";"+CLIP3+";"+ttest+";"+read_count
  with open("clusters.gff","w") as f:
    for line in open("/homes2/marek/annotation/genome_50bp_window.gff","r"):
      cluster_name = line.split()[8]
      if cluster_name in enrichments:
        line = line.strip().split()
        line[7] = ";".join(enrichments[cluster_name].split(";")[0:-1])
        line[5] = enrichments[cluster_name].split(";")[-1]
        line = "\t".join(line)
        f.write("{}\n".format(line))




def main():
  create_RPKM_table()
  calculate_enrichment()
  filter_enrichment()

if __name__ == "__main__":
  main()



