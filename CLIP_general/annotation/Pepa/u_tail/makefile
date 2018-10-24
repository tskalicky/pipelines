## Workflow in this folder
## The script calculates average U-tail length in uadd
## It uses trimmed U's - so the length may not be correct
## When 0 U's was trimmed, it skips the read
#################################################
DATA_MUT := $(wildcard ~/Projects/dis3l2/samples/mut*/*.gz)
DATA_RIB := $(wildcard ~/Projects/dis3l2/ribosomal/samples/mut*/*.gz)

all: tail_length_histogram.pdf


results_histogram_u_tail.txt: count_T.pl $(DATA_MUT) $(DATA_RIB)
	cd rRNA_data && ./count_T.pl
	perl count_T.pl
	tail -n 1 rRNA_data/results_average_u_tail.txt >> results_average_u_tail.txt
	paste results_histogram_u_tail.txt rRNA_data/results_histogram_u_tail.txt > tmp
	mv tmp results_histogram_u_tail.txt

## NOTE: This scripts counts only T's withnout a mismatch!

tail_length_histogram.pdf: results_histogram_u_tail.txt tail_length.gp
	gnuplot tail_length.gp
	epspdf tail_length_histogram.eps
	pdfcrop tail_length_histogram.pdf
	mv tail_length_histogram-crop.pdf tail_length_histogram.pdf

