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
```
make gmer_counter
make gmer_caller
```

* Usage   
First one has to prepare the database of specific k-mers for each allele of polymorphism of interest.  
k-mer databases are available from http://bioinfo.ut.ee/FastGT/  
  
K-mers are counted from raw reads in FASTQ file using program gmer_counter:  
```
gmer_counter -db DATABASE FASTQ_FILE(S) > COUNTS_FILE.txt
```
  
Then the genotypes are called using program gmer_caller:  
```
gmer_caller COUNTS_FILE.txt > GENOTYPE_FILE.txt
```

Genotype file can be converted to VCF format:  
```
generate_vcf.pl GENOTYPE_FILE.txt > GENOTYPE_FILE.vcf
```
   
Additional options for gmer_counter:   
```
    -db DATABASE     - SNP/KMER database file
    -dbb DBBINARY    - binary database file
    -w FILENAME      - write binary database to file
    -32              - use 32-bit integeres for counts (default 16-bit)
    --max_kmers NUM  - maximum number of kmers per node
    --header         - print header row
    --total          - print the total number of kmers per node
    --unique         - print the number of nonzero kmers per node
    --kmers          - print individual kmer counts (default if no other output)
    --compile_index  - Add read index to database and write it to file
    --distribution NUM  - print kmer distribution (up to given number)
    --num_threads    - number of worker threads (default 24)
    --low_memory     - optimize for low memory usage
    -D               - increase debug level
```

Additional options for gmer_caller:   
```
    --training_size NUM - Use NUM markers for training (default 100000)
    --runs NUMBER       - Perfom NUMBER runs of model training (use 0 for no training)
    --num_threads NUM   - Use NUM threads (min 1, max 32, default 16)
    --header            - Print table header
    --non_canonical     - Output non-canonical genotypes
    --prob_cutoff       - probability cutoff for calling genotype (default 0)
    --alternatives      - Print probabilities of all alternative genotypes
    --info              - Print information about individual
    --no_genotypes      - Print only summary information, not actual genotypes
    --model TYPE        - Model type (full, diploid, haploid)
    --params PARAMS     - Model parameters (error, p0, p1, p2, coverage, size, size2)
    --coverage NUM      - Average coverage of reads
    -D                  - increase debug level
```