#!/usr/bin/perl

use strict;

#
# Find regions with overrepresented repeats
#
# Copyright (C) Lauris Kaplinski and University of Tartu
# You can distribute and modify this script freely
#
# Usage:
# find_regions.pl OVERREPFILE FASTAFILE MINLEN MINMOVAVG [MAXLEN]
#
# OVERREPFILE - a list of k-mers that are overrepresented in given sequence (glistquery output)
# FASTAFILE - the sequence to be analyzed in FastA format
# MINLEN - minimum repeat length
# MINMOVAVG - overrepresentation moving average cutoff
# MAXLEN - the maximum length of repeat
#

my $overrep_file = $ARGV[0];
my $fasta_file = $ARGV[1];
my $min_len = $ARGV[2];
my $min_movavg = $ARGV[3];
my $max_len = $ARGV[4];
if ($max_len eq "") {
  $max_len = 10000;
}

my $wordlen = 16;

my %overrep;
printf (STDERR "Loading oligo file (%s)\n", $overrep_file);
open (my $ifs, $overrep_file) or die ("Cannot open $overrep_file");
while (my $line = <$ifs>) {
  chomp($line);
  $line =~ s/\r//;
  my @tokenz = split ("\t", $line);
  my $oligo = $tokenz[0];
  my $count = $tokenz[1];
  #printf (STDERR "%s\t%s\n", $oligo, $count);
  $overrep{$oligo} = $count;
}
close ($ifs);
printf (STDERR "Done\n");

printf (STDERR "Loading FastA file (%s)\n", $fasta_file);
open (my $ifs, $fasta_file) or die ("Cannot open $fasta_file");
my $line = <$ifs>;
my $seq = "";
while ($line = <$ifs>) {
  chomp ($line);
  $seq .= $line;
}
close ($ifs);
printf (STDERR "Done\n");

my $idx = 1;
my $nwords = length ($seq) - $wordlen;
printf (STDERR "Sequence contains %d words\n", $nwords);
my $start = -1;
my $end = -1;
my $best = 0;
my $sum = 0;
for (my $i = 0; $i < $nwords; $i++) {
  my $word = substr ($seq, $i, $wordlen);
  my $count = $overrep{$word};
  if ($count >= $min_movavg) {
    # This is overrepresented
    $sum += $count;
    if ($start < 0) {
      # Start region
      $start = $i;
      $end = $i + 32;
      printf (STDERR "Starting region at %d", $i);
    } else {
      # Continue region
      $end = $i + 32;
    }
  } else {
    # This is not
    if ($start >= 0) {
      # Test if region ends
      my $len = $i + 1 - $start;
      my $movavg = $sum / $len;
      if ($movavg < $min_movavg) {
        # Region ends
        $len = $end - $start;
        printf (STDERR " ending at %d length %d\n", $i, $len);
        $movavg = $sum / ($len - 31);
        if (($len >= $min_len) && ($len <= $max_len)) {
          my $reg = substr ($seq, $start, $len);
          printf (STDOUT ">Repeat_%d %d-%d length %d avg %.2f\n", $idx, $i, $i + $len, $len, $movavg);
          printf (STDOUT "%s\n", $reg);
          $idx += 1;
        }
        $sum = 0;
        $start = -1;
      }
    } else {
      $sum = 0;
    }
  }
}

