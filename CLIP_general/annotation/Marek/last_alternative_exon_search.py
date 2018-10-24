#!/usr/bin/python


def return_last_exon(transcript_info):
  last_exon = None
  for record in transcript_info:
    if record.split()[2] == "exon" and record.split()[1] == "protein_coding":
      last_exon = record
  return last_exon
  
def last_exons_overlap(last_exons_of_this_gene):
  last_exons = [x for x in last_exons_of_this_gene if x is not None]
  positions = []
  for i in filter(lambda x: x != None,last_exons):
    start = int(i.split()[3])
    end = int(i.split()[4])
    positions.append([x for x in range(start,end)])
  for i in positions:	
    for j in positions:
      if not (set(i) & set(j)):
        #return last_exons[0].split()[9].replace(";","").replace("\"","")
        return last_exons

with open("/homes2/marek/annotation/Ensembl_GRCh37.75.gtf") as f:
  line = f.readline().split("\t")
  while line[0]:  
    if line[0].startswith("#"):
      line = f.readline().split("\t")
    else:
      if line[1] == "protein_coding" and line[2] == "gene":
        ENSG = line[8].split()[1]
        gene_info = []
        last_exons_of_this_gene = []
        line = f.readline().split("\t")
        while line[8].split()[1] == ENSG: 
          ENST = line[8].split()[3]
          transcript_info = []
          while line[8].split()[3] == ENST:
            transcript_info.append("\t".join(line))
            line = f.readline().split("\t")
          
          gene_info.append(transcript_info)
          last_exons_of_this_gene.append(return_last_exon(transcript_info))
        
        x = last_exons_overlap(last_exons_of_this_gene)
        if x is not None:
          for i in x:
            print i,
           
      else:
        line = f.readline().split("\t")
    if not line[0]:
      break     
      
    
