# FastGT genotype caller
Copyright (C) University of Tartu 2015-2017  
  
FastGT is a fast and accurate genotype caller from sequencing data.  
It uses Empirical Bayes classifier for calling genotypes from the counts of polymorphism-specific k-mers.  
FastGT can detect both canonical (diploid for autosomes and female X, haploid for male X and Y) and
non-canonical (0-4 alleles) genotypes.  
  
* Binaries  
FastGT has two binaries - gmer_counter and gmer_caller. Pre-compiled binaries are in directory 'bin'.  
  
* Compilation  
Change int src subdirectory and type:
'''
make gmer_counter
make gmer_caller
'''

* Usage

First one has to prepare the database of specific k-mers for each allele of polymorphism of interest.
k-mer databases are available from http://bioinfo.ut.ee/FastGT/  
  
K-mers are counted from raw reads in FASTQ file using program gmer_counter:  
'''
gmer_counter -db DATABASE FASTQ_FILE(S) > COUNTS_FILE.txt
'''
  
Then the genotypes are called using program gmer_caller:  
'''
gmer_caller COUNTS_FILE.txt > GENOTYPE_FILE.txt
'''

Genotype file can be converted to VCF format:
'''
generate_vcf.pl GENOTYPE_FILE.txt > GENOTYPE_FILE.vcf
'''