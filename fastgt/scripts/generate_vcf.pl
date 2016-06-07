#!/usr/bin/perl

use strict;
use warnings;

#
# generate_vcf CALLS_FILE
#

my $calls = $ARGV[0];

my $sex = 0;

my @d = localtime();
my $year = $d[5];
my $month = $d[4];
my $day = $d[3];

printf (STDOUT "##fileformat=VCFv4.1\n");
printf (STDOUT "##fileDate=%4d%02d%02d\n", 1900 + $year, 1 + $month, $day);
printf (STDOUT "##source=%s\n", $calls);
printf (STDOUT "##reference=HumanNCBI37_UCSC\n");
printf (STDOUT "##phasing=none\n");
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=TI,Number=.,Type=String,Description="Transcript ID">
##INFO=<ID=GI,Number=.,Type=String,Description="Gene ID">
##INFO=<ID=EXON,Number=0,Type=Flag,Description="Exon Region">
##INFO=<ID=FC,Number=.,Type=String,Description="Functional Consequence">
printf (STDOUT "##FILTER=<ID=q20,Description=\"Quality below 20\">\n");
printf (STDOUT "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
printf (STDOUT "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");

printf (STDOUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t****\n");

open (my $ifs, $calls) or die $!;
while (my $line = <$ifs>) {
  chomp ($line);
  if (substr ($line, 0, 4) eq "#Sex") {
    if (substr ($line, 5, 1) eq "M") {
      $sex = 1;
    }
  }
  next if (substr ($line, 0, 1) eq "#");
  my @t = split ("\t", $line);
  my $gt = $t[1];
  my $prob = $t[2];
  my @tt = split (":", $t[0]);
  my $chr = $tt[0];
  my $pos = $tt[1];
  my $id = $tt[2];
  @tt = split ("/", $tt[3]);
  my $ref = $tt[0];
  my $alt = $tt[1];
  my $rc = $t[3];
  my $ac = $t[4];

  my $a0 = 0;
  my $a1 = 0;
  if (($sex == 0) || (($chr ne "Y") && ($chr ne "X"))) {
    if ($gt eq "AA") {
      # 
    } elsif ($gt eq "AB") {
      $a1 = 1;
    } elsif ($gt eq "BB") {
      $a0 = 1;
      $a1 = 1;
    } else {
      #
    }
  } else {
    if ($gt eq "A") {
      #
    } elsif ($gt eq "B") {
      $a0 = 1;
      $a1 = 1;
    } else {
      #
    }
  }

  my $qual = "*";
  my $filter = "*";
  my $info = "*";
  my $format = "GT:GQ";
  printf (STDOUT "%s\t%s\t%s\t%s\t%s", $chr, $pos, $id, $ref, $alt);
  printf (STDOUT "\t%s\t%s\t%s\t%s", $qual, $filter, $info, $format);
  printf (STDOUT "\t%s/%s:%s", $a0, $a1, $rc + $ac);

  printf (STDOUT "\t%s\n", $gt);
}
