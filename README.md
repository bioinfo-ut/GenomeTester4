# GenomeTester4
A toolkit for performing set operations - union, intersection and complement - on k-mer lists.  
Copyright (C) University of Tartu 2015-2016

Building GenomeTester4 binaries:  
Change into subdirectory 'src' and type:
```
make clean
make
```
# FastGT 
An alignment-free genotype caller.  
Copyright (C) University of Tartu 2015-2017  
  
FastGT is implemented as a sub-project of GenomeTester4, sharing some source files.  
Source files are inside the GenomeTester4 distribution.  
FastGT has two binaries - gmer_counter and gmer_caller.  
  
Building FastGT binaries:  
Change into subdirectory 'src\' and type:  
```
make gmer_counter
make gmer_caller
```
* k-mer lists and user manual for FastGT are available from http://bioinfo.ut.ee/FastGT/
