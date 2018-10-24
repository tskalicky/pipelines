#!/usr/bin/python

import sys

dict = {}
for line in sys.stdin:
	try:
		if line.strip().split("\t")[9] in dict:
			dict[line.strip().split("\t")[9]] = dict[line.strip().split("\t")[9]] + line.strip().split("\t")[8].split(";")[1]+";"+line.strip().split("\t")[0]+";"+line.strip().split("\t")[3]+";"+line.strip().split("\t")[6]+"\t"
		else:
			dict[line.strip().split("\t")[9]] = line.strip().split("\t")[8].split(";")[1]+";"+line.strip().split("\t")[0]+";"+line.strip().split("\t")[3]+";"+line.strip().split("\t")[6]+"\t"
	except IndexError:
		pass

print "#Ensembl_gene_id","Alus_on_plus", "Alus_on_minus", "Alus_details"

for key in dict:
	cnt_plus = 0
	cnt_minus = 0
	for item in dict[key].split("\t"):
		try:
			if item.split(";")[3] == "+":
				cnt_plus += 1
			if item.split(";")[3] == "-":
				cnt_minus += 1
		except IndexError:
			pass
	print key.replace(";","").replace("\"",""),cnt_plus,cnt_minus,dict[key]
