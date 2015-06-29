#!/usr/bin/perl

use strict;

#
# Collate (semi)identical repeats
#
# Copyright (C) Lauris Kaplinski and University of Tartu
# You can distribute and modify this script freely
#
# Usage:
# collate_repeats.pl BLASTFILE FASTAFILE
#
# BLASTFILE - BLAST tabulated output file, aligning overrepresented regions against themselves
# FASTAFILE - the original FastA file
#

my $blast_file = $ARGV[0];
my $fasta_file = $ARGV[1];

my %iddict;
printf (STDERR "Loading BLAST file (%s)\n", $blast_file);
open (my $ifs, $blast_file) or die ("Cannot open $blast_file");
my @lines = <$ifs>;
close ($ifs);
printf (STDERR "Done\n");

foreach (@lines) {
  my $line = $_;
  chomp($line);
  $line =~ s/\r//;
  my @tokenz = split ("\t", $line);
  my $dbid = $tokenz[0];
  $iddict{$dbid} = 1;
}

my @ids;
my %names;
my %seqs;
printf (STDERR "Loading FastA file (%s)\n", $fasta_file);
open (my $ifs, $fasta_file) or die ("Cannot open $fasta_file");
while (my $line = <$ifs>) {
  chomp($line);
  $line =~ s/\r//;
  my $name = substr ($line, 1);;
  my @tokenz = split (" ", $name);
  my $id = $tokenz[0];
  $line = <$ifs>;
  chomp($line);
  $line =~ s/\r//;
  my $seq = $line;
  printf (STDERR "Adding %s\n", $id);
  push (@ids, $id);
  $names{$id} = $name;
  $seqs{$id} = $seq;
}
close ($ifs);
printf (STDERR "Done\n");

@ids = sort(@ids);

my %id;
foreach (@ids) {
  my $dbid = $_;
  foreach (@lines) {
    my $line = $_;
    chomp($line);
    $line =~ s/\r//;
    my @tokenz = split ("\t", $line);
    if ($tokenz[0] eq $dbid) {
      my $dblen = $tokenz[1];
      my $qid = $tokenz[2];
      next if (!exists ($names{$qid}));
      my $qlen = $tokenz[3];
      my $ident = $tokenz[4];
      my $len = $tokenz[5];
      if (($ident > 90) && (abs ($dblen/$qlen - 1) < 0.05) && (abs ($dblen/$len - 1) < 0.05)) {
        #printf (STDOUT "%s and %s are the same (%d %d %d %.2f)\n", $dbid, $qid, $dblen, $qlen, $len, $ident);
        #printf (STDOUT "%s %s\n", $id{$dbid}, $id{$qid});
        if ($id{$qid} eq "") {
          if ($id{$dbid} eq "") {
            $id{$dbid} = $dbid;
          }
          if ($qid ne $dbid) {
            $id{$qid} = $dbid;
            #printf (STDOUT "%s->%s\n", $qid, $dbid);
          }
        }
      }
    }
  }
}

# Print unique
#foreach (@ids) {
#  my $dbid = $_;
#  if ($id{$dbid} eq $dbid) {
#    printf (STDOUT ">%s\n%s\n", $names{$dbid}, $seqs{$dbid});
#  }
#}

# Print groups
foreach (@ids) {
  my $dbid = $_;
  if ($id{$dbid} eq $dbid) {
    printf (STDOUT "\nGroup %s\n\n", $dbid);
    printf (STDOUT ">%s\n%s\n\n", $names{$dbid}, $seqs{$dbid});
    foreach (@lines) {
      my $line = $_;
      chomp($line);
      $line =~ s/\r//;
      my @tokenz = split ("\t", $line);
      if ($tokenz[0] eq $dbid) {
        my $dblen = $tokenz[1];
        my $qid = $tokenz[2];
        next if (!exists ($names{$qid}));
        my $qlen = $tokenz[3];
        my $ident = $tokenz[4];
        my $len = $tokenz[5];
        if ($id{$qid} ne $dbid) {
          printf (STDOUT ">%s\n%s\n", $names{$qid}, $seqs{$qid});
          #printf (STDOUT "%s\n", $line);
        }
      }
    }
  }
}
