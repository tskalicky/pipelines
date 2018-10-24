#!/usr/bin/python


def return_last_exon(transcript_info):
  last_exon = None
  for record in transcript_info:
    if record.split()[2] == "exon" and record.split()[1] == "protein_coding":
      last_exon = record
  return last_exon
 

 
with open("/homes2/marek/annotation/Ensembl_GRCh37.75.gtf") as f:
  line = f.readline().split("\t")
  while line[0]:
    if not line[0].startswith("#"):
      if line[1] == "protein_coding" and line[2] == "gene":
        ENSG = line[8].split()[1]
        gene_info = []
        last_exons_of_this_gene = []
        line = f.readline().split("\t")
        while line[8].split()[1] == ENSG: 
          ENST = line[8].split()[3]
          transcript_info = []
          while line[8].split()[3] == ENST :
            transcript_info.append("\t".join(line))
            line = f.readline().split("\t")
                      
          #gene_info.append(transcript_info)
          last_exons_of_this_gene.append(return_last_exon(transcript_info))
         
        for i in last_exons_of_this_gene:
          if i:         
            print i, 
          
      else:
        line = f.readline().split("\t")			# PASS
    else:
      line = f.readline().split("\t")			# PASS
    if not line[0]:					# while loop condition
      break
