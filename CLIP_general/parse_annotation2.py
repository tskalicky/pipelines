#!/usr/bin/python

import os,sys

annotation = "/homes2/marek/annotation/Ensembl_GRCh37.75.gtf"

def sort_positions(list,strand):
	if strand == "+":
		return sorted(list)
	elif strand == "-":
		return reversed(sorted(list))		# Reverse if on - strand

def report_segment(positions,transcript,feature_type,output_gff_parsed):
	chrom = transcript[0]
        strand = transcript[6]
        transcript_start = float(transcript[3])
        transcript_end = float(transcript[4])
        gene_name = transcript[9].strip(";").strip("\"")
        type = transcript[1]
	i = 0
	with open(output_gff_parsed,"a") as f:
		for position in sort_positions(positions,strand):
			if feature_type == "5UTR":
				correction = 0			# 5'UTR is in range <0,1>
			elif feature_type == "CDS":
				correction = 1			# CDS is shifted by 1
			elif feature_type == "3UTR":
				correction = 2			# 3' UTR is shifted by 2
			
			percentage =  correction + (float(i) / len(positions))
			f.write( "{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom,position,strand,percentage,gene_name,type))
			i += 1

def parse_gff_to_percentages(annotation,output_gff_parsed):
	exons = []
	CDS = []
	UTR5 = []
	UTR3 = []
	total_3UTR = 0
	total_5UTR = 0
	total_CDS = 0
	for line in open(annotation,"r"):
		data = line.strip().split()
	
		if line[0] == "#":						# Pass comment lines 
			pass
	
		elif data[2] == "exon":                                         # Read exon information
	                for position in range(int(data[3]),int(data[4])):
	                        exons.append(position)
	
	        elif data[2] == "CDS":                                          # Read CDS information
	                for position in range(int(data[3]),int(data[4])):
	                        CDS.append(position)

		elif line.strip().split()[2] == "transcript":			# New transcript encountered - FUN BEGINS !!!!
			if exons and CDS:					# Make sure info about CDS and exons is present
				strand = transcript[6]
				for position in exons:
			                if position in CDS:					# CDS and exons positions
	       	 		                pass
	               			if strand == "+" and position < min(CDS):		# 5'UTR positions on +
			                        UTR5.append(position)
			                elif strand == "+" and position > max(CDS):		# 3' UTR positions on +
			                        UTR3.append(position)
			                elif strand == "-" and position > max(CDS):		# Vice vesa due to orientation
			                        UTR5.append(position)
			                elif strand == "-" and position < min(CDS):
			                        UTR3.append(position)
		
					total_3UTR += len(UTR3)					# Count total coverage of 3'UTR, 5'UTR and CDS for normalization
					total_5UTR += len(UTR5)
					total_CDS = len(CDS)
	
				report_segment(UTR5,transcript,"5UTR",output_gff_parsed)
				report_segment(CDS,transcript,"CDS",output_gff_parsed)  
				report_segment(UTR3,transcript,"3UTR",output_gff_parsed)  
	

		

			transcript = data
			exons = []
			CDS = []
			UTR5 = []
			UTR3 = []  # Reset variables after new transcript is encountered



def check_complete_parsing(output_gff_parsed):
	for line in open(output_gff_parsed):
		if "chrMT" in line:
			return True
	return False

def summarize_parsed_annotation(output_gff_parsed):
	pass

def main():
	output_gff_parsed = annotation.split("/")[-1].split(".")[0]+"_parsed.sgrs"
 	print sys.argv
	if "--new_parse" in sys.argv:
		if os.path.isfile(output_gff_parsed):
			os.remove(output_gff_parsed)
		parse_gff_to_percentages(annotation,output_gff_parsed)
	else:
		if not os.path.isfile(output_gff_parsed) or not check_complete_parsing(output_gff_parsed):
			sys.stderr.write("*** Error: Ensembl annotation file not parsed or parsed incompletely, please re-run with --new_parse argument ***\n")
			sys.exit(2)		

	summarize_parsed_annotation(output_gff_parsed)
		



	

if __name__ == "__main__":
	main()











