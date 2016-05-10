#ifndef __BINOMIAL_H__
#define __BINOMIAL_H__

unsigned int combination (unsigned int n, unsigned int k);
double combination_k_r_1 (unsigned int k, double r);

double poisson (unsigned int k, unsigned int lambda);

double dbinom (unsigned int x, unsigned int n, double p);

double dnbinom (unsigned int x, double k, double p);
double dnbinom_mu (unsigned int x, double k, double mu);

#endif
