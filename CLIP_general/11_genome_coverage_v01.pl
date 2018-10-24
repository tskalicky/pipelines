#!/bin/env perl
 ### GENOME COVERAGE ###
#
use strict;
use warnings;
#
# variables 
my $v ="GRCh38.91";
my $input = "/home/tomas/ownCloud/CEITEC_lab/genomes/annotation/my_annot_DB/bed";
 
my $exon_file       = "$input/${v}.mRNA_exons.sorted.merged.bed";
my $intergenic_file = "$input/${v}.intergenic.bed";
my $intron_file     = "$input/${v}.mRNA_introns.bed";
 
my $exon_coverage       = coverage($exon_file);
my $intergenic_coverage = coverage($intergenic_file);
my $intron_coverage     = coverage($intron_file);
 
my $total = $exon_coverage + $intergenic_coverage + $intron_coverage;
 
printf "Exon: %.2f\n", $exon_coverage*100/$total;
printf "Intron: %.2f\n", $intron_coverage*100/$total;
printf "Intergenic: %.2f\n", $intergenic_coverage*100/$total;
 
sub coverage {
   my ($infile) = @_;
   my $coverage = 0;
   # open(IN,'-|',"zcat $infile") || die "Could not open $infile: $!\n";
   open(IN,'-|',"cat $infile") || die "Could not open $infile: $!\n";
   while(<IN>){
      chomp;
      my ($chr, $start, $end) = split(/\t/);
      my $c = $end - $start;
      $coverage += $c;
   }
   close(IN);
   return($coverage);
}
 
exit(0);