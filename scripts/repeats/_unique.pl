#!/usr/bin/perl

use strict;

# Collate unique repeats according to BLAST table
#
# unique.pl FASTAFILE BLASTFILE

my $fasta_file = $ARGV[0];
my $blast_file = $ARGV[1];

my @ids;
my @seqs;
my %incl;
open (my $ifs, $fasta_file) or die $!;
while (my $line = <$ifs>) {
  chomp($line);
  $line =~ s/\r//;
  my @t = split (" ", $line);
  my $id = substr ($t[0], 1);
  $line = <$ifs>;
  chomp($line);
  $line =~ s/\r//;
  my $seq = $line;
  next if (length ($seq) > 2000);
  push (@ids, $id);
  push (@seqs, $seq);
  #printf (STDOUT ">%s\n%s\n", $id, $seq);
  $incl{$id} = 1;
}
close ($ifs);

open ($ifs, $blast_file) or die $!;
while (my $line = <$ifs>) {
  chomp($line);
  $line =~ s/\r//;
  my @t = split ("\t", $line);
  my $id0 = $t[0];
  my $len0 = $t[1];
  my $id1 = $t[2];
  my $len1 = $t[3];
  my $ident = $t[4];
  my $alen = $t[5];
  next if ($id0 eq $id1);
  next if ($id0 gt $id1);
  next if ($incl{$id0} == 0);
  next if ($incl{$id1} == 0);
  next if ($ident < 90);
  next if (abs (($alen - $len0) / $alen) > 0.1);
  next if (abs (($alen - $len1) / $alen) > 0.1);
  #printf (STDOUT "%s is identical to %s\n", $id0, $id1);
  #printf (STDOUT "%s\n", $line);
  $incl{$id1} = 0;
}
close ($ifs);

for (my $i = 0; $i < @ids; $i++) {
  my $id = $ids[$i];
  my $seq = $seqs[$i];
  next if ($incl{$id} == 0);
  printf (STDOUT ">%s\n%s\n", $id, $seq);
}
