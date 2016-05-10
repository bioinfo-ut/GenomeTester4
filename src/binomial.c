#define __BINOMIAL_C__

#include <math.h>
#include <assert.h>

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
poisson (unsigned int k, unsigned int lambda)
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
  unsigned int c = combination (n, x);
  double p0 =  pow (p, x);
  double p1 =  pow (1 - p, n - x);
  return c * p0 * p1;
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
