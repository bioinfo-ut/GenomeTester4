# GenomeTester4
A toolkit for performing set operations - union, intersection and complement - on k-mer lists.  
Copyright (C) University of Tartu 2015-2016.  
Please cite: Kaplinski L, Lepamets M, Remm M. (2015). GenomeTester4: a toolkit for performing basic set operations – union, intersection and complement on k-mer lists. GigaScience, 4:58  

Building GenomeTester4 binaries:  
Change into subdirectory 'src' and type:
```
make clean
make
```
# FastGT 
An alignment-free genotype caller.  
Copyright (C) University of Tartu 2015-2017  
Please cite: Pajuste F-D, Kalpinski L, Möls M, Puurand T, Lepamets M, Remm M. (2017). FastGT: an alignment-free method for calling common SNVs directly from raw sequencing reads. Scientific Reports, 7:2537  

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
