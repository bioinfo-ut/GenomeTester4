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
FastGT utilizes two binaries from the GenomeTester4 package: gmer_counter and gmer_caller.  

Building FastGT binaries:  
```
git clone https://github.com/bioinfo-ut/GenomeTester4.git
cd GenomeTester4/
cd src
make gmer_counter
make gmer_caller
```
* k-mer lists and user manual for FastGT are available from http://bioinfo.ut.ee/FastGT/
See also README.FastGT.md


# KATK
A genotype caller for rare and de novo variants.  
Copyright (C) University of Tartu 2016-2021

KATK is implemented as a sub-project of GenomeTester4, sharing some source files.  
Source files are inside the GenomeTester4 distribution.  
KATK utilizes two binaries from the GenomeTester4 package: gmer_counter and gassembler.  

Building KATK binaries:  
```
git clone https://github.com/bioinfo-ut/GenomeTester4.git
cd GenomeTester4/
cd src
make gmer_counter
make gassembler
```
* KATK data files required for calling variants from the entire exome are available from http://bioinfo.ut.ee/KATK/
See also README.KATK.md
