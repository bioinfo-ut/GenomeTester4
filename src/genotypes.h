#ifndef __GENOTYPES_H__
#define __GENOTYPES_H__

enum {
  X, A, B, AA, AB, BB, AAA, AAB, BBA, BBB, AAAA, AAAB, BBBA, AABB, BBBB, NUM_GENOTYPES
};

void genotype_probabilities (double a[], float avg_maf, unsigned int count_a, unsigned int count_b, double l_error, double p_0, double p_1, double p_2, double lambda, double size, double size2);

#endif
