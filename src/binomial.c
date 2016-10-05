#define __BINOMIAL_C__

#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "binomial.h"

void
find_factos (unsigned int n, unsigned int *factors, unsigned int stride, unsigned int *nfactors)
{
  unsigned int i, j;
  memset (factors, 0, n * stride * sizeof (unsigned int));
  nfactors[0] = 0;
  nfactors[1] = 0;
  for (i = 2; i < n; i++) nfactors[i] = 1;
  factors[2 * stride] = 2;
  for (i = 2; i < n; i++) {
    if (nfactors[i] == 1) {
      /* i is prime */
      factors[i * stride] = i;
      for (j = 2; j < n; j++) {
        unsigned int k = j * i;
        if (k < n) {
          unsigned int l;
          /* Fill in d[k] */
          for (l = 0; l < nfactors[j]; l++) {
            factors[k * stride + l] = factors[j * stride + l];
          }
          factors[k * stride + l] = i;
          nfactors[k] = l + 1;
        }
      }
    }
  }
}

#define MAX_PROD 4096

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
    static double *logs = NULL;
    static double *log_sums = NULL;
    double s;
    if (!logs) {
      unsigned int i, j;
      logs = (double *) malloc (MAX_PROD * sizeof (double));
      log_sums = (double *) malloc (MAX_PROD * sizeof (double));
      for (i = 1; i < MAX_PROD; i++) {
        logs[i] = log (i);
        log_sums[i] = logs[i];
        for (j = 2; j < i; j++) {
          log_sums[i] += logs[j];
        }
      }
    }
    s = log_sums[n] - log_sums[n - k] - log_sums[k];
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
  double c = combinations (n, x);
  double p0 =  pow (p, x);
  double p1 =  pow (1 - p, n - x);
  return c * p0 * p1;
}

double
log_dbinom (unsigned int x, unsigned int n, double p)
{
  double c = log_combinations (n, x);
  double p0 = log (p) * x;
  double p1 = log (1 - p) * (n - x);
  return c + p0 + p1;
}

double
combination_k_r_1 (unsigned int k, double r)
{
  if (!k) return 1;
  double val = (k + r - 1) / k;
  while (k > 1) {
    k -= 1;
    val *= (k + r - 1);
    val /= k;
  }
  return val;
}

double
dnbinom (unsigned int x, double r, double p)
{
  /* unsigned int c = combination (x + r - 1, x); */
  double c = combination_k_r_1 (x, r);
  double p0 =  pow (p, x);
  double p1 =  pow (1 - p, r);
  if ((p0 == 0) || (p1 == 0)) return 0;
  double val = c * p0 * p1;
  assert (!isnan (c));
  assert (!isnan (val));
  return val;
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
