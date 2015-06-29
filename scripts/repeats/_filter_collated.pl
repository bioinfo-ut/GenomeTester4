#!/usr/bin/perl

use strict;

# Keep only groups with minimum number of members
#
# Copyright (C) Lauris Kaplinski and University of Tartu
# You can distribute and modify this script freely
#
# Usage:
# filter_collated.pl GROUP_FILE  MIN_NUM_MATCHES
#


my $group_file = $ARGV[0];
my $min_num = $ARGV[1];

open (my $ifs, $group_file) or die ("Cannot open $group_file");

my $gidx = 0;
my $block = "";
my $num_members = 0;
while (my $line = <$ifs>) {
  #chomp($line);
  #$line =~ s/\r//;
  if (substr ($line, 0, 5) eq "Group") {
    if (($gidx > 0) && ($num_members > $min_num)) {
      printf (STDOUT "%s", $block);
    }
    $num_members = 0;
    $block = "";
    $gidx += 1;
  } elsif (substr ($line, 0, 1) eq ">") {
    $num_members += 1;
  }
  $block .= $line;
}
if (($gidx > 0) && ($num_members > $min_num)) {
  printf (STDOUT "%s", $block);
}
close ($ifs);
