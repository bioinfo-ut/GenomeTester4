GenomeTester4 package

GenomeTester4 is a toolkit for creating and manipulating k-mer lists. It
contains 3 programs: glistmaker, glistcompare and glistquery. It is
developed by Department of Bioinformatics, University of Tartu and
distributed under GPL version 3.0 (or later).


1. Quick usage for the impatient



  glistmaker|glistcompare|glistquery -h

Prints out quick description of command line arguments



  glistmaker FASTA_FILES -w WORD_LENGTH -o OUTPUT_NAME

Generates a list of all unique k-mers in input files with their frequencies. The output file is
named OUTPUT_NAME_WORD_LENGTH_0_0.list



  glistquery LIST_FILE

Prints out all k-mers and their frequencies in list file



  glistquery LIST_FILE -q WORD

Prints out the frequency of given word in list file



  glistquery LIST_FILE -f WORD_FILE

Prints out the frequencies of all words in WORD_FILE in list file



  glistcompare LIST_FILE_1 LIST_FILE_2 --union -o OUTPUT_NAME

Generates union of all words in both input lists. The frequencies in final
list are sums of both input frequencies.



  glistcompare LIST_FILE_1 LIST_FILE_2 --intersection -o OUTPUT_NAME

Generates intersection of words in both input lists. The frequencies in
final list are smaller frequencies of either list



  glistcompare LIST_FILE_1 LIST_FILE_2 --difference -o OUTPUT_FILE

Generates complement (words in the first file and NOT in second file) of
words in both input lists. The frequencies in final list are the original
frequencies in the first list.



2. Compiling

GenomeTester4 is written in standard C. The only external dependency should
be pthreads library that is standard in all Linux systems.
Binaries compiled with full optimization are included in directory "bin".
If you for whatever reason have to compile these manually, just enter into
src subdirectory and type:

make clean
make
