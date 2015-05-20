#!/usr/bin/perl
use strict;

# Front-end for making intersection of more than two list files
# Maarja Lepamets

use File::Basename;

my $nlists = scalar @ARGV;
my $list1;
my $list2;
my $out = "intrs";



my $k = 1;
while ($nlists != 1) {
    
  if ($nlists == scalar @ARGV) {
  
    print STDERR "mkdir -p ${out}_${k}\n";
    system "mkdir -p ${out}_${k}";
    
    my $i = 0;
    while ($i < $nlists) {
      $list1 = $ARGV[$i];
      if ($i == $nlists - 1) {
        print STDERR "cp $list1 ${out}_${k}/copy_".basename($list1) . "\n";
        system "cp $list1 ${out}_${k}/copy_".basename($list1);
        last;
      }
      $list2 = $ARGV[$i + 1];
  
      print STDERR "/mambakodu/maarjalepamets/programs/glistcompare $list1 $list2 -o ${out}_${k}/$i"."_".($i + 1)." -i\n";
      system "/mambakodu/maarjalepamets/programs/glistcompare $list1 $list2 -o ${out}_${k}/$i"."_".($i + 1)." -i";  
  
      $i += 2;
      
    }
      
  } else {
    
    #system "rm -r -f ${out}_".($k - 2);
    
    my $i = 0;
    my $dir = "${out}_".($k - 1); 
    my @files = <$dir/*>;
    $nlists = scalar @files;
    
    if ($nlists == 2) {
      $list1 = $files[$i];
      $list2 = $files[$i + 1];
      print STDERR "/mambakodu/maarjalepamets/programs/glistcompare $list1 $list2 -o intrs -i\n";
      system "/mambakodu/maarjalepamets/programs/glistcompare $list1 $list2 -o intrs -i";
      #system "rm -r -f ${out}_".($k - 1);
      last;
    }
    
    print STDERR "mkdir -p ${out}_${k}\n";
    system "mkdir -p ${out}_${k}";
    while ($i < $nlists) {
      $list1 = $files[$i];
      if ($i == $nlists - 1) {
        
        print STDERR "cp $list1 ${out}_${k}/copy_".basename($list1) . "\n";
        system "cp $list1 ${out}_${k}/copy_".basename($list1);
        last;
      }
      $list2 = $files[$i + 1];
      
      print STDERR "/mambakodu/maarjalepamets/programs/glistcompare $list1 $list2 -o ${out}_${k}/$i"."_".($i + 1)." -i\n";
      system "/mambakodu/maarjalepamets/programs/glistcompare $list1 $list2 -o ${out}_${k}/$i"."_".($i + 1)." -i";
      $i += 2;
      
    }
    
  }
  
  $nlists = int($nlists / 2 + 0.5);
  $k += 1;
}
