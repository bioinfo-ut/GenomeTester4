#define __BINOMIAL_C__

#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "binomial.h"

#define MAX_PROD 16384

static double *t_log_factorial_d = NULL;
static float *t_log_factorial_f = NULL;
static double *t_log_sums_d = NULL;
static float *t_log_sums_f = NULL;

void
init_combination_tables (void)
{
  /* Double log factorial */
  if (!t_log_factorial_d) {
    unsigned int i;
    t_log_factorial_d = (double *) malloc (MAX_PROD * sizeof (double));
    t_log_factorial_d[0] = 0;
    for (i = 1; i < MAX_PROD; i++) {
      t_log_factorial_d[i] = t_log_factorial_d[i - 1] + log (i);
    }
  }
  /* Float log factorial */
  if (!t_log_factorial_f) {
    unsigned int i;
    t_log_factorial_f = (float *) malloc (MAX_PROD * sizeof (float));
    t_log_factorial_f[0] = 0;
    for (i = 1; i < MAX_PROD; i++) {
      t_log_factorial_f[i] = t_log_factorial_f[i - 1] + logf (i);
    }
  }
  /* Double sum of logarithms */
  if (!t_log_sums_d) {
    unsigned int i, j;
    double *logs = (double *) malloc (MAX_PROD * sizeof (double));
    t_log_sums_d = (double *) malloc (MAX_PROD * sizeof (double));
    for (i = 1; i < MAX_PROD; i++) {
      logs[i] = log (i);
      t_log_sums_d[i] = logs[i];
      for (j = 2; j < i; j++) {
        t_log_sums_d[i] += logs[j];
      }
    }
  }
  /* Float sum of logarithms */
  if (!t_log_sums_f) {
    unsigned int i, j;
    float *logs = (float *) malloc (MAX_PROD * sizeof (float));
    t_log_sums_f = (float *) malloc (MAX_PROD * sizeof (float));
    for (i = 1; i < MAX_PROD; i++) {
      logs[i] = logf (i);
      t_log_sums_f[i] = logs[i];
      for (j = 2; j < i; j++) {
        t_log_sums_f[i] += logs[j];
      }
    }
  }
}

static double
log_factorial (unsigned int v)
{
  double dv = v;
  double val = 0;
  if (!t_log_factorial_d) init_combination_tables ();
  while (v >= MAX_PROD) {
    val += log (dv);
    dv -= 1;
    v -= 1;
  }
  val += t_log_factorial_d[v];
  return val;
}

static float
log_factorial_f (unsigned int v)
{
  float dv = v;
  float val = 0;
  if (!t_log_factorial_f) init_combination_tables ();
  while (v >= MAX_PROD) {
    val += log (dv);
    dv -= 1;
    v -= 1;
  }
  val += t_log_factorial_f[v];
  return val;
}

double
combinations_d (unsigned int n, unsigned int k)
{
  return exp (log_combinations_d (n, k));
}

double
log_combinations_d (unsigned int n, unsigned int k)
{
  if (!k || (k == n)) return 0;
  if (k == 1) return log (n);
  if (n >= MAX_PROD) {
    double num = log (n), den = log (k);
    while (k > 1) {
      n -= 1;
      k -= 1;
      num += log (n);
      den -= log (k);
    }
    return num - den;
  } else {
    double s;
    if (!t_log_sums_d)  init_combination_tables ();
    s = t_log_sums_d[n] - t_log_sums_d[n - k] - t_log_sums_d[k];
    return s;
  }
}

float
log_combinations_f (unsigned int n, unsigned int k)
{
  if (!k || (k == n)) return 0;
  if (k == 1) return logf (n);
  if (n >= MAX_PROD) {
    float num = logf (n), den = logf (k);
    while (k > 1) {
      n -= 1;
      k -= 1;
      num += logf (n);
      den -= logf (k);
    }
    return num - den;
  } else {
    float s;
    if (!t_log_sums_f)  init_combination_tables ();
    s = t_log_sums_f[n] - t_log_sums_f[n - k] - t_log_sums_f[k];
    return s;
  }
}

unsigned int
combination (unsigned int n, unsigned int k)
{
  if (!k || (k == n)) return 1;
  double num = n, den = k;
  while (k > 1) {
    n -= 1;
    k -= 1;
    num *= n;
    den *= k;
  }
  return (unsigned int) (num / den + 0.5);
}

double
poisson (unsigned int k, double lambda)
{
  double dk = k;
  double dl = lambda;
  double p = exp (-dl);
  while (k > 0) {
    p *= dl;
    p /= dk;
    dk -= 1;
    k -= 1;
  }
  return p;
}

double
dbinom (unsigned int x, unsigned int n, double p)
{
  if ((x == 0) && (p == 0)) return 1;
  if ((x == n) && (p == 1)) return 1;
  double c = combinations_d (n, x);
  double p0 =  pow (p, x);
  double p1 =  pow (1 - p, n - x);
  return c * p0 * p1;
}

double
log_dbinom (unsigned int x, unsigned int n, double p)
{
  if ((x == 0) && (p == 0)) return 0;
  if ((x == n) && (p == 1)) return 0;
  double c = log_combinations_d (n, x);
  assert (!isnan (c));
  double p0 = log (p) * x;
  assert (!isnan (p0));
  double p1 = log (1 - p) * (n - x);
  assert (!isnan (p1));
  return c + p0 + p1;
}

double
combination_k_r (unsigned int k, double r)
{
  return exp (log_combination_k_r (k, r));
}

double
log_combination_k_r (unsigned int k, double r)
{
  if (!k) return 0;
  double val = lgamma (k + r) - lgamma (r) - log_factorial (k);
  return val;
}

float
log_combination_k_r_f (unsigned int k, float r)
{
  if (!k) return 0;
  float val = lgammaf (k + r) - lgammaf (r) - log_factorial_f (k);
  return val;
}

double
dnbinom (unsigned int x, double size, double p)
{
  double c, p0, p1, val;
  c = log_combination_k_r (x, size);
  assert (!isnan (c));
  p0 = log (p) * x;
  p1 = log (1 - p) * size;
  val = c + p0 + p1;
  assert (!isnan (val));
  return exp (val);
}

double
dnbinom_mu (unsigned int x, double size, double mu)
{
  double p;
  if (size <= 0) return 0;
  if (mu <= 0) return 0;
  p = mu / (size + mu);
  return dnbinom (x, size, p);
}

double
dnbinom_precalc (unsigned int x, double size, double p, double log_n_combinations)
{
  double c, p0, p1, val;
  c = log_n_combinations;
  p0 =  log (p) * x;
  p1 =  log (1 - p) * size;
  val = exp (c + p0 + p1);
  assert (!isnan (val));
  return val;
}

float
dnbinom_precalc_f (unsigned int x, float size, float p, float log_n_combinations)
{
  float c, p0, p1, val;
  c = log_n_combinations;
  p0 =  logf (p) * x;
  p1 =  logf (1 - p) * size;
  val = expf (c + p0 + p1);
  assert (!isnan (val));
  return val;
}

double
dnbinom_mu_precalc (unsigned int x, double size, double mu, double log_n_combinations)
{
  double p;
  if (size <= 0) return 0;
  if (mu <= 0) return 0;
  p = mu / (size + mu);
  return dnbinom_precalc (x, size, p, log_n_combinations);
}

float
dnbinom_mu_precalc_f (unsigned int x, float size, float mu, float log_n_combinations)
{
  float p;
  if (size <= 0) return 0;
  if (mu <= 0) return 0;
  p = mu / (size + mu);
  return dnbinom_precalc_f (x, size, p, log_n_combinations);
}

#define pi 3.14159265359

double
PDF (double x, double mu, double sigma)
{
  return 1 / (sigma * sqrt (2 * pi)) * exp (-((x - mu) * (x - mu)) / (2 * sigma * sigma));
}

static double
PDF_0 (double x)
{
  return (1 / sqrt (2 * pi)) * exp (-(x * x) / 2);
}

#define b0 0.2316419
#define b1 0.319381530
#define b2 -0.356563782
#define b3 1.781477937
#define b4 -1.821255978
#define b5 1.330274429

double
CDF (double x, double mu, double sigma)
{
  double t;
  x -= mu;
  x /= sigma;
  if (x < 0) return 1 - CDF (-x, 0, 1);
  t = 1 / (1 + b0 * x);
  return 1 - PDF_0 (x) * (b1 * t + b2 * t * t + b3 * t * t * t + b4 * t * t * t * t + b5 * t * t * t * t * t);
}
