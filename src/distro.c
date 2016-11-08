#define __DISTRO_C__

#include <assert.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "binomial.h"
#include "sequence.h"
#include "simplex.h"
#include "utils.h"
#include "wordmap.h"
#include "wordtable.h"

static unsigned int debug = 2;

typedef struct _List List;
typedef struct _DistroData DistroData;

static float find_coeffs (float values[], const DistroData *dd, unsigned int ncounts);

static void print_coverages (unsigned int nlists, const char *lists[]);
static void get_distro (wordmap *map, unsigned int *dist, unsigned int min, unsigned int max);
static unsigned int find_median_coverage (wordmap *map);

#define EXPECTED_KMER_COV_BITS 14

#define MAX_LISTS 1000
#define MAP_BLOCK 1000000
#define WORDLENGTH 25

#define DISTRO_SIZE 60
#define HAPLOID_VALUE 10
#define NCOEFFS 6

struct _List {
  const char *filename;
  
}

struct _DistroData {
  /* Total number of samples */
  unsigned int n_samples;
  /* Observed values */
  unsigned int o[DISTRO_SIZE];
  /* Values >= DISTRO_SIZE */
  unsigned int high;
  /* Sum of observed values */
  unsigned int sum_o;
};

static float
logit (float p)
{
  return logf (p / (1 - p));
}

static float
logit_1 (float a)
{
  return 1 / (1 + expf (-a));
}

static double lambdas[NCOEFFS + 1];
static double *log_combos;

int
main (int argc, const char *argv[])
{
  const char *lfile;
  const char **lists;
  unsigned int *coverages;
  unsigned int nlists;
  const unsigned char *cdata;
  size_t csize, cpos;
  wordmap **maps;
  unsigned int i, kmer_idx;
  unsigned long long *lptr;
  unsigned long long current, next;

  lambdas[0] = 0.1;
  for (i = 1; i <= NCOEFFS; i++) {
    lambdas[i] = i * HAPLOID_VALUE;
  }
  
  lists = (const char **) malloc (1000 * 8);
  coverages = (unsigned int *) malloc (1000 * 4);
  nlists = 0;
    
  lfile = argv[1];
  cdata = (const unsigned char *) mmap_by_filename (lfile, &csize);
  cpos = 0;
  while (cpos < csize) {
    unsigned long long s, e;
    char b[1024];
    unsigned int size;
    s = cpos;
    while ((cdata[s] <= ' ') && (s < csize)) s += 1;
    if (s >= csize) break;
    e = s;
    while ((cdata[e] > ' ') && (e < csize)) e += 1;
    memcpy (b, cdata + s, e - s);
    b[e - s] = 0;
    s = e;
    while ((cdata[s] <= ' ') && (s < csize)) s += 1;
    if (s >= csize) break;
    e = s;
    while ((cdata[e] > ' ') && (e < csize)) e += 1;
    size = strtol ((const char *) cdata + s, NULL, 10);
    lists[nlists] = strdup (b);
    coverages[nlists] = size;
    nlists += 1;
    if (nlists > MAX_LISTS) break;
    cpos = e;
  }

  /* Precaulculate combinations */
  log_combos = (double *) malloc ((nlists + 1) * sizeof (double));
  for (i = 0; i <= nlists; i++) {
    log_combos[i] = log_combinations_d (nlists, i);
    fprintf (stderr, "%u %u %g\n", nlists, i, log_combos[i]);
  } 
  
  maps = (wordmap **) malloc (nlists * 8);
  for (i = 0; i < nlists; i++) {
    maps[i] = wordmap_new (lists[i], 0);
  }
  lptr = (unsigned long long *) malloc (nlists * 8);
  memset (lptr, 0, nlists * 8);
  /* Pick first current */
  current = WORDMAP_WORD (maps[0], 0);
  for (i = 0; i < nlists; i++) {
    if ((lptr[i] < maps[i]->header->nwords) && (WORDMAP_WORD (maps[i], lptr[i]) < current)) {
      current = WORDMAP_WORD (maps[i], lptr[i]);
    }
  }
  next = WORDMAP_WORD (maps[0], 0 + MAP_BLOCK);
  kmer_idx = 0;
  while (current < next) {
    DistroData dd = { 0 };
    char c[33];
    float coeffs[NCOEFFS];
    float result;
    unsigned int add_zero;
    
    dd.n_samples = nlists;
    /* Calculate counts */
    for (i = 0; i < nlists; i++) {
      if ((lptr[i] < maps[i]->header->nwords) && (WORDMAP_WORD (maps[i], lptr[i]) == current)) {
        unsigned int freq = WORDMAP_FREQ (maps[i], lptr[i]);
        freq = (freq * 2 * HAPLOID_VALUE + coverages[i] / 2) / coverages[i];
        if (freq < DISTRO_SIZE) {
          dd.o[freq] += 1;
          dd.sum_o += 1;
        } else {
          dd.high += 1;
        }
        lptr[i] += 1;
      }
    }
    add_zero = dd.n_samples - dd.sum_o - dd.high;
    dd.o[0] += add_zero;
    dd.sum_o += add_zero;
    /* fprintf (stdout, "Sum %g\n", sum); */
    if ((kmer_idx >= 0) && (dd.sum_o >= 400)) {
      float dsum;

      result = find_coeffs (coeffs, &dd, DISTRO_SIZE);
      word2string (c, current, WORDLENGTH);
      fprintf (stdout, "%s :", c);
      
      for (i = 1; i < DISTRO_SIZE; i++) {
        fprintf (stdout, " %u", dd.o[i]);
      }
      fprintf (stdout, " (%u)", dd.high);
      fprintf (stdout, ": %u\n", dd.sum_o);
      fprintf (stdout, " %g :", result);
      dsum = 0;
      for (i = 0; i < NCOEFFS; i++) {
        dsum += logit_1 (coeffs[i]);
        fprintf (stdout, " %.3f", logit_1 (coeffs[i]));
      }
      fprintf (stdout, " : %.2f\n", dsum);
    }
    kmer_idx += 1;
    /* Pick new current */
    current = 0xffffffffffffffff;
    for (i = 0; i < nlists; i++) {
      if ((lptr[i] < maps[i]->header->nwords) && (WORDMAP_WORD (maps[i], lptr[i]) < current)) {
        current = WORDMAP_WORD (maps[i], lptr[i]);
      }
    }
  }
  
  /*print_coverages (argc, argv);*/
  
  return 1;
}

/* Generates expected distrbution 0..DISTRO_SIZE, last element is sum of all above */

static double
generate_distribution (int nvalues, const float values[], double *e)
{
  double sum_values, sum_q;
  unsigned int i, j;

  sum_values = 0;
  for (j = 0; j < nvalues; j++) sum_values += logit_1 (values[j]);
  sum_q = 0;
  for (i = 0; i < DISTRO_SIZE; i++) {
    double pp = 0;
    for (j = 0; j < nvalues; j++) {
      pp += logit_1 (values[j]) * poisson (i, lambdas[j]);
    }
    if (pp < 1e-60) pp = 1e-60;
    e[i] = pp;
    sum_q += pp;
  }
  e[DISTRO_SIZE] = sum_values - sum_q;
  return sum_q;
}

static float
distance (int nvalues, const float values[], void *data)
{
  const DistroData *dd = (const DistroData *) data;
  unsigned int i;
  double e[DISTRO_SIZE + 1];
  double S_Q, S_C, D, delta;

  /* Generate theoretical distribution */
  S_Q = generate_distribution (nvalues, values, e);
  /* Divergence */
  D = 1;
  for (i = 1; i < DISTRO_SIZE; i++) {
    double p = dbinom (dd->o[i], dd->n_samples, e[i]);
    if (p < 1e-60) p = 1e-60;
    D -= log (p);
    assert (isfinite (D));
  }

  S_C = 0;
  for (i = 0; i < nvalues; i++) {
    S_C += logit_1 (values[i]);
  }
  delta = S_Q * dd->n_samples - dd->sum_o;
  D *= (1 + abs (delta) / 100);
  
  if (debug > 2) {
    for (i = 1; i < DISTRO_SIZE; i++) {
      fprintf (stdout, "%.3f : ", e[i]);
    }
    fprintf (stdout, "\n");
    for (i = 0; i < NCOEFFS; i++) {
      fprintf (stdout, "%.2f ", logit_1 (values[i]));
    }
    fprintf (stdout, "%g\n", D);
  }
  return D;
}

static float
find_coeffs (float values[], const DistroData *dd, unsigned int ncounts)
{
  float deltas[NCOEFFS];
  unsigned int i;
  float result;
  double q[DISTRO_SIZE + 1];

  for (i = 0; i < NCOEFFS; i++) {
    values[i] = logit (0.5f / NCOEFFS);
    deltas[i] = values[i] / 40;
    /* values[i] = logit (1e-8f); */
  }
  result = downhill_simplex (NCOEFFS, values, deltas, 1e-5f, 5, 500, distance, (void *) dd);
  
  if (debug > 1) {
    fprintf (stdout, "\n");
    for (i = 0; i < NCOEFFS; i++) {
      fprintf (stdout, "%.2f ", lambdas[i]);
    }
    fprintf (stdout, "\n");
    for (i = 0; i < NCOEFFS; i++) {
      fprintf (stdout, "%.2f ", logit_1 (values[i]));
    }
    fprintf (stdout, "\n");
    generate_distribution (NCOEFFS, values, q);
    fprintf (stdout, "EXP:");
    for (i = 0; i < DISTRO_SIZE + 1; i++) {
      fprintf (stdout, " %2u", (unsigned int) (q[i] * dd->n_samples + 0.5));
    }
    fprintf (stdout, "\n");
    fprintf (stdout, "OBS:");
    for (i = 0; i < DISTRO_SIZE + 1; i++) {
      fprintf (stdout, " %2u", dd->o[i]);
    }
    fprintf (stdout, "\n");
  }
  return result;
}

static void
print_coverages (unsigned int nlists, const char *lists[])
{
  unsigned int i;
  for (i = 1; i < nlists; i++) {
    wordmap *map;
    unsigned int cov, j;
    unsigned int dist[30];
    unsigned int max;
    
    map = wordmap_new (lists[i], 1);
    if (!map) continue;

    get_distro (map, dist, 11, 40);
    if (debug) {
      for (j = 12; j <= 40; j++) {
        fprintf (stderr, "%u\t%u\n", j, dist[j - 11]);
      }
    }
    max = 11;
    for (j = 12; j <= 40; j++) {
      if (dist[j - 11] > dist[max - 11]) max = j;
    }
    if ((max > 11) && (max < 40)) {
      fprintf (stdout, "%s\t%u\n", lists[i], max);
    }
    fprintf (stderr, "%s max = %u freq = %u\n", lists[i], max, dist[max - 11]);
#if 0
    cov = find_median_coverage (lists[i]);
    fprintf (stdout, "%s\t%u\n", lists[i], cov);
#endif
    wordmap_release (map);
  }
}

static void
get_distro (wordmap *map, unsigned int *dist, unsigned int min, unsigned int max)
{
  unsigned long long i;
  memset (dist, 0, (max - min + 1) * 4);
  for (i = 0; i < map->header->nwords; i++) {
    unsigned long long freq = WORDMAP_FREQ (map, i);
    if ((freq >= min) && (freq <= max)) {
      dist[freq - min] += 1;
    }
  }
}

static unsigned int
find_median_coverage (wordmap *map)
{
  unsigned long long i;
  unsigned int min, max, med;
  min = 1000000000;
  max = 0;
  for (i = 0; i < map->header->nwords; i++) {
    unsigned long long freq = WORDMAP_FREQ (map, i);
    if (freq < min) min = freq;
    if (freq > max) max = freq;
  }
  med = (min + max) / 2;
  while (max > min) {
    unsigned long long above = 0, below = 0, equal;
    for (i = 0; i < map->header->nwords; i++) {
      unsigned long long freq = WORDMAP_FREQ (map, i);
      if (freq > med) above += 1;
      if (freq < med) below += 1;
    }
    equal = map->header->nwords - above - below;
    if (debug > 1) fprintf (stderr, "Min %u med %u max %u, below %llu, equal %llu, above %llu\n", min, med, max, below, equal, above);
    /* Special case: min == med, max == med + 1 */
    if (max == (min + 1)) {
      if (above > (below + equal)) {
        /* Max is true median */
        med = max;
      }
      break;
    }
    if (above > below) {
      if ((above - below) < equal) break;
      min = med;
    } else if (below > above) {
      if ((below - above) < equal) break;
      max = med;
    } else {
      break;
    }
    med = (min + max) / 2;
  }
  return med;
}

