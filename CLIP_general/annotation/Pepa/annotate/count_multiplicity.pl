#!/usr/bin/perl
## USAGE: perl count_multiplicity.pl mut1 uadd
## Adds multiplicity of annotation to the mapping
#######################################################

use strict;
use warnings;
use Data::Dumper;
use File::Copy;

if( scalar(@ARGV) != 2 ) {
  warn "Not enough arguments. Use: mutX filename";
  exit 1;
}

my $home     = "$ENV{'HOME'}/Projects/dis3l2/samples";
my @samples  = qw(mut1 mut2 mut3);
my %multi    = (); # annotation multiplicity
my $mutX   = $ARGV[0];
my $filename = $ARGV[1];


## Load multiplicity
my $path = "$home/$mutX/$filename";
print "Counting multiplicity: $path\n";
open( F, "gunzip -c $path |" ) or die();
while (my $line = <F>) {
  chomp($line);
  my @s = split('\t', $line);
  my ($ID) = $s[3] =~ /([^\-]+)-\d+-/;
  $multi{$ID} += 1;
}
close(F);


## Print results
open( F, "gunzip -c $path |" ) or die();
open( OUT, ">$path.temp" );
while (my $line = <F>) {
  chomp($line);
  my @s = split('\t', $line);
  my ($ID) = $s[3] =~ /([^\-]+)-\d+-/;
  $s[7] *= $multi{$ID};
  my $out = join("\t", @s);
  print OUT "$out\n";
}
close(F);
close(OUT);


print("$path\n");
($path) = $path =~ /(.+)\.gz/;
move("$path.gz.temp", "$path");
system("gzip -f $path");
