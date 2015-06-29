#!/usr/bin/perl

use strict;

# Parse final list of BLAST vs all chromosomes
#
# unique.pl FASTAFILE BLASTFILE TARGET

my $fasta_file = $ARGV[0];
my $blast_file = $ARGV[1];
my $tgt = $ARGV[2];

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
  push (@ids, $id);
  push (@seqs, $seq);
  $incl{$id} = 1;
}
close ($ifs);

my %tcount;
my %ocount;
open ($ifs, $blast_file) or die $!;
while (my $line = <$ifs>) {
  chomp($line);
  $line =~ s/\r//;
  my @t = split ("\t", $line);
  my $qid = $t[0];
  my $dbid = $t[1];
  if ($dbid eq $tgt) {
    $tcount{$qid} += 1;
  } else {
    $ocount{$qid} += 1;
  }
}
close ($ifs);

for (my $i = 0; $i < @ids; $i++) {
  my $id = $ids[$i];
  my $seq = $seqs[$i];
  #printf (STDOUT "%s T:%s O:%s\n", $id, $tcount{$id}, $ocount{$id});
  if (($tcount{$id} > 0) && ($ocount{$id} == 0)) {
    printf (STDOUT ">%s %s:%s\n%s\n", $id, $tgt, $tcount{$id}, $seq);
  }
}
