#ifndef __BINOMIAL_H__
#define __BINOMIAL_H__

/* Prepare lookup tables for small combination numbers */
/* Needed if simultaneous querys from multiple threads may occur */
void init_combination_tables (void);

unsigned int combination (unsigned int n, unsigned int k);
double combinations_d (unsigned int n, unsigned int k);
double log_combinations_d (unsigned int n, unsigned int k);
float log_combinations_f (unsigned int n, unsigned int k);
double combination_k_r_1 (unsigned int k, double r);
double log_combination_k_r_1 (unsigned int k, double r);
float log_combination_k_r_f (unsigned int k, float r);

double poisson (unsigned int k, double lambda);

double dbinom (unsigned int x, unsigned int n, double p);
double log_dbinom (unsigned int x, unsigned int n, double p);

double dnbinom (unsigned int x, double k, double p);
double dnbinom_mu (unsigned int x, double size, double mu);

/* Calculates negative binomial distribution with precalculated log number of combinations */
double dnbinom_mu_precalc (unsigned int x, double size, double mu, double log_n_combinations);
float dnbinom_mu_precalc_f (unsigned int x, float size, float mu, float log_n_combinations);

double PDF (double x, double mu, double sigma);
double CDF (double x, double mu, double sigma);

#endif
