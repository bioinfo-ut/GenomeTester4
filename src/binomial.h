#ifndef __BINOMIAL_H__
#define __BINOMIAL_H__

/* Find all factors of numbers 0...n-1 and write to table */
void find_factos (unsigned int n, unsigned int *factors, unsigned int stride, unsigned int *nfactors);

unsigned int combination (unsigned int n, unsigned int k);
double combinations_d (unsigned int n, unsigned int k);
double log_combinations_d (unsigned int n, unsigned int k);
double combination_k_r_1 (unsigned int k, double r);

double poisson (unsigned int k, double lambda);

double dbinom (unsigned int x, unsigned int n, double p);
double log_dbinom (unsigned int x, unsigned int n, double p);

double dnbinom (unsigned int x, double k, double p);
double dnbinom_mu (unsigned int x, double k, double mu);

double PDF (double x, double mu, double sigma);
double CDF (double x, double mu, double sigma);

#endif
