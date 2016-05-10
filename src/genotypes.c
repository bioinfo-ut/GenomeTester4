#define __GENOTYPES_C__

#include <math.h>

#include "binomial.h"

#include "genotypes.h"

void
genotype_probabilities (double a[], float pB, unsigned int var1, unsigned int var2, double l_viga, double p_0, double p_1, double p_2, double lambda, double size, double size2)
{
  double q0, q1;
  double p[NUM_GENOTYPES];
  double p_A_alleel, p_B_alleel;
  double p_lisa, p_lisa1, p_lisa2;

  p_B_alleel = pB;
  p_A_alleel = 1 - p_B_alleel;

  p[X] = p_0;
  p[A] = p_A_alleel * p_1;
  p[B] = p_B_alleel * p_1;
  p[AA] = p_A_alleel * p_A_alleel * p_2;
  p[AB] = 2 * p_A_alleel * p_B_alleel * p_2;
  p[BB] = p_B_alleel * p_B_alleel * p_2;

  /* tõenäosus näha midagi enam kui 2-s korduses */
  p_lisa = 1 - p_0 - p_1 - p_2;
  if (p_lisa >= 0) {
    /* 3-se genotüübi (näteks AAB) tõenäosus */
    p_lisa1 = (-1 + sqrtf (1 + 4 * p_lisa)) / 2;
    /* 4-se genotüübi (näiteks AAAB) tõenäosus */
    p_lisa2 = p_lisa1 * p_lisa1;
  } else {
    p_lisa1 = 0;
    p_lisa2 = 0;
  }
  
  /* fixme: mu or not? */
  p[AAA] = dbinom (3, 3, p_A_alleel) * p_lisa1;
  p[BBB] = dbinom (0, 3, p_A_alleel) * p_lisa1;
  p[AAB] = dbinom (2, 3, p_A_alleel) * p_lisa1;
  p[BBA] = dbinom (1, 3, p_A_alleel) * p_lisa1;
  p[AAAA] = dbinom (4, 4, p_A_alleel) * p_lisa2;
  p[BBBB] = dbinom (0, 4, p_A_alleel) * p_lisa2;
  p[AAAB] = dbinom (3, 4, p_A_alleel) * p_lisa2;
  p[AABB] = dbinom (2, 4, p_A_alleel) * p_lisa2;
  p[BBBA] = dbinom (1, 4, p_A_alleel) * p_lisa2;

  /* a_0=dnbinom(var2, mu=l_viga, size=size+size2*l_viga)*dnbinom(var1, mu=l_viga, size=size+size2*l_viga)*p_0 */
  q0 = dnbinom_mu (var1, size + size2 * l_viga, l_viga);
  q1 = dnbinom_mu (var2, size + size2 * l_viga, l_viga);
  a[X] = q0 * q1 * p[0];
    
  /* a_A=dnbinom(var1, mu=lambda/2, size=size+size2*lambda/2)*dnbinom(var2, mu=l_viga, size=size+size2*l_viga)*p_A */
  q0 = dnbinom_mu (var1, size + size2 * lambda / 2, lambda / 2);
  q1 = dnbinom_mu (var2, size + size2 * l_viga, l_viga);
  a[A] = q0 * q1 * p[A];
    
  /* a_B=dnbinom(var2, mu=lambda/2, size=size+size2*lambda/2)*dnbinom(var1, mu=l_viga, size=size+size2*l_viga)*p_B */
  q0 = dnbinom_mu (var2, size + size2 * lambda / 2, lambda / 2);
  q1 = dnbinom_mu (var1, size + size2 * l_viga, l_viga);
  a[B] = q0 * q1 * p[B];

  /* a_AA=dnbinom(var1, mu=lambda, size=size+size2*lambda*1)*dnbinom(var2, mu=l_viga, size=size+size2*l_viga)*p_AA */
  q0 = dnbinom_mu (var1, size + size2 * lambda, lambda);
  q1 = dnbinom_mu (var2, size + size2 * l_viga, l_viga);
  a[AA] = q0 * q1 * p[AA];

  /* a_AB=dnbinom(var1, mu=lambda/2, size=size+size2*lambda/2)*dnbinom(var2, mu=lambda/2, size=size+size2*lambda/2)*p_AB */
  q0 = dnbinom_mu (var1, size + size2 * lambda / 2, lambda / 2);
  q1 = dnbinom_mu (var2, size + size2 * lambda / 2, lambda / 2);
  a[AB] = q0 * q1 * p[AB];

  /* a_BB=dnbinom(var2, mu=lambda, size=size+size2*lambda)*dnbinom(var1, mu=l_viga, size=size+size2*l_viga)*p_BB */
  q0 = dnbinom_mu (var2, size + size2 * lambda, lambda);
  q1 = dnbinom_mu (var1, size + size2 * l_viga, l_viga);
  a[BB] = q0 * q1 * p[BB];

  /* a_AAA=dnbinom(var1, mu=lambda*1.5, size=size+size2*lambda*1.5)*dnbinom(var2, mu=l_viga, size=size+size2*l_viga)*p_AAA */
  q0 = dnbinom_mu (var1, size + size2 * lambda * 1.5, lambda * 1.5);
  q1 = dnbinom_mu (var2, size + size2 * l_viga, l_viga);
  a[AAA] = q0 * q1 * p[AAA];
  
  /* a_AAB=dnbinom(var1, mu=lambda, size=size+size2*lambda)*dnbinom(var2, mu=lambda/2, size=size+size2*lambda/2)*p_AAB */
  q0 = dnbinom_mu (var1, size + size2 * lambda, lambda);
  q1 = dnbinom_mu (var2, size + size2 * lambda / 2, lambda / 2);
  a[AAB] = q0 * q1 * p[AAB];
  
  /* a_BBA=dnbinom(var1, mu=lambda/2, size=size+size2*lambda/2)*dnbinom(var2, mu=lambda, size=size+size2*lambda)*p_BBA */
  q0 = dnbinom_mu (var1, size + size2 * lambda / 2, lambda / 2);
  q1 = dnbinom_mu (var2, size + size2 * lambda, lambda);
  a[BBA] = q0 * q1 * p[BBA];

  /* a_BBB=dnbinom(var1, mu=l_viga, size=size+size2*l_viga)*dnbinom(var2, mu=lambda*1.5, size=size+size2*lambda*1.5)*p_BBB */
  q0 = dnbinom_mu (var1, size + size2 * l_viga, l_viga);
  q1 = dnbinom_mu (var2, size + size2 * lambda * 1.5, lambda * 1.5);
  a[BBB] = q0 * q1 * p[BBB];

  /* a_AAAA=dnbinom(var1, mu=lambda*2, size=size+size2*lambda*2)*dnbinom(var2, mu=l_viga, size=size+size2*l_viga)*p_AAAA */
  q0 = dnbinom_mu (var1, size + size2 * lambda * 2, lambda * 2);
  q1 = dnbinom_mu (var2, size + size2 * l_viga, l_viga);
  a[AAAA] = q0 * q1 * p[AAAA];

  /* a_AAAB=dnbinom(var1, mu=lambda*1.5, size=size+size2*lambda*1.5)*dnbinom(var2, mu=lambda/2, size=size+size2*lambda/2)*p_AAAB */
  q0 = dnbinom_mu (var1, size + size2 * lambda * 1.5, lambda * 1.5);
  q1 = dnbinom_mu (var2, size + size2 * lambda / 2, lambda / 2);
  a[AAAB] = q0 * q1 * p[AAAB];

  /* a_BBBA=dnbinom(var1, mu=lambda/2, size=size+size2*lambda/2)*dnbinom(var2, mu=lambda*1.5, size=size+size2*lambda*1.5)*p_BBBA */
  q0 = dnbinom_mu (var1, size + size2 * lambda / 2, lambda / 2);
  q1 = dnbinom_mu (var2, size + size2 * lambda * 1.5, lambda * 1.5);
  a[BBBA] = q0 * q1 * p[BBBA];
  
  /* a_AABB=dnbinom(var1, mu=lambda, size=size+size2*lambda)*dnbinom(var2, mu=lambda, size=size+size2*lambda)*p_AABB */
  q0 = dnbinom_mu (var1, size + size2 * lambda, lambda);
  q1 = dnbinom_mu (var2, size + size2 * lambda, lambda);
  a[AABB] = q0 * q1 * p[AABB];
  
  /* a_BBBB=dnbinom(var1, mu=l_viga, size=size+size2*l_viga)*dnbinom(var2, mu=lambda*2, size=size+size2*lambda*2)*p_BBBB */
  q0 = dnbinom_mu (var1, size + size2 * l_viga, l_viga);
  q1 = dnbinom_mu (var2, size + size2 * lambda * 2, lambda * 2);
  a[BBBB] = q0 * q1 * p[BBBB];
}

