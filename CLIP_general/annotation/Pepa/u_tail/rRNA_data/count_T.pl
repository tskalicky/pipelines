#!/usr/bin/perl
# Counts U-tail length in given RNA type
##############################################

use strict;
use warnings;

my $home      = "$ENV{'HOME'}/Projects/dis3l2/samples";
my @samples   = qw(mut1 mut2 mut3);
my %totalT    = ();
my %totalCopy = ();
my %histogram = ();

sub loadData {
  my $path = shift;
  open( F, "gunzip -c $path |") or die $!;
  while (my $line = <F>) {
    my @s = split("\t", $line);
    # s[7] = transcript multiplicity
    # s[14]= RNA type
    if( $s[3] =~ /\w+\-(\d+)\-(\w+)/ ) {
      my $cp = $1;
      my $ts = $2;
      $totalT{"rRNA"} += length($ts) * $cp / $s[7];
      $totalCopy{"rRNA"} += $cp / $s[7];
      $histogram{"rRNA"}{length($ts)} += $cp / $s[7];
    }
  }
  close(F);
}

foreach my $sample (@samples) {
  my $path = "$home/../ribosomal/samples/$sample/foreground.uadd.sel.bed.gz";
  loadData($path);
}

## Print results: average
open(OUT, ">results_average_u_tail.txt");
print OUT "# Average (uadd) U-tail length in:\n";
foreach my $key (keys %totalT) {
  my $average = $totalT{$key} / $totalCopy{$key};
  $average = sprintf("%.2f", $average);
  print OUT "$key: $average\n";
}
close(OUT);

## Print results: histogram
open(OUT, ">results_histogram_u_tail.txt");
my @allKeys = keys %histogram;
my $tempAllKeys = join("\t", @allKeys);
print OUT "\n$tempAllKeys\n";
for my $i (0..25) {
  foreach my $key (@allKeys) {
    if( defined $histogram{$key}{$i} ) {
      my $len = sprintf("%.2f", $histogram{$key}{$i} / $totalCopy{$key} * 100);
      print OUT "$len";
    } else {
      print OUT "0";
    }
  }
  print OUT "\n";
}
close(OUT);
