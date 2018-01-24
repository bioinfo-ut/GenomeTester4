
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "binomial.h"
#include "genotypes.h"
#include "utils.h"
#include "simplex.h"
#include "thread-pool.h"

static unsigned int debug = 0;

#define MAX_THREADS 32

#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))

AosoraThreadPool *pool = NULL;

static const char *gt[] = { "-", "A", "B", "AA", "AB", "BB", "AAA", "AAB", "BBA", "BBB", "AAAA", "AAAB", "BBBA", "AABB", "BBBB" };

enum { L_VIGA, P_0, P_1, P_2, LAMBDA, SIZE, SIZE2 };

/* SNV k-mer calls */
typedef struct _SNPCall SNPCall;
struct _SNPCall {
  unsigned int line;
  unsigned short counts[2];
};

static float distanceL3 (int ndim, const float params[], void *data);

static unsigned int get_pair_median (const unsigned char *cdata, unsigned long long csize, const unsigned char *lines[], unsigned int autos[], unsigned int nautos);

typedef struct _L3Data L3Data;
typedef struct _L3Optim L3Optim;

struct _L3Optim {
  L3Data *l3;
  unsigned int first_call;
  unsigned int n_calls;
  double sum;
};

struct _L3Data {
  float *params;
  unsigned int n_calls;
  unsigned int *var1;
  unsigned int *var2;
  float pB;
  unsigned int n_threads;
  AosoraThreadPool *pool;
  L3Optim optims[MAX_THREADS];
};

#define MIN_P (1.0f / 8192)

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

static float
logit_clamped (float p, float min, float max)
{
  if (p <= min) {
    p = min;
  } else if (p >= max) {
    p = max;
  } else {
    p = (p - min) / (max - min);
  }
  return logit (p);
}

static float
logit_1_clamped (float a, float min, float max)
{
  a = logit_1 (a);
  a = min + (max - min) * a;
  return a;
}

static void
print_params (const float params[], const char *prefix, FILE *ofs)
{
  float l_viga, p_0, p_1, p_2, lambda, size, size2;

  l_viga = logit_1_clamped (params[0], MIN_P, 0.1f);
  p_0 = logit_1_clamped (params[1], MIN_P, 1 - MIN_P);
  p_1 = logit_1_clamped (params[2], MIN_P, 1 - MIN_P);
  p_2 = logit_1_clamped (params[3], MIN_P, 1 - MIN_P);
  lambda = expf (params[4]);
  size = params[5];
  size2 = -expf (params[6]);
  fprintf (ofs, "%s %g %g %g %g %g %g %g\n", prefix, l_viga, p_0, p_1, p_2, lambda, size, size2);
}

static unsigned int
count_lines (const unsigned char *cdata, unsigned long long csize)
{
  unsigned long long cpos;
  unsigned int nlines;
  nlines = 0;
  cpos = 0;
  while (cpos < csize) {
    while ((cpos < csize) && (cdata[cpos] != '\n')) cpos += 1;
    if (cdata[cpos] == '\n') nlines += 1;
    cpos += 1;
  }
  return nlines;
}

static unsigned int
parse_lines (const unsigned char *lines[], const unsigned char *cdata, unsigned long long csize)
{
  unsigned long long cpos;
  unsigned int nlines;
  nlines = 0;
  cpos = 0;
  while (cpos < csize) {
    lines[nlines] = cdata + cpos;
    while ((cpos < csize) && (cdata[cpos] != '\n')) cpos += 1;
    if (cdata[cpos] == '\n') nlines += 1;
    cpos += 1;
  }
  return nlines;
}

SNPCall *
parse_calls (const unsigned char *lines[], const unsigned int indices[], unsigned int nautos, unsigned int pair_median)
{
  unsigned int i;
  SNPCall *calls = (SNPCall *) malloc (nautos * sizeof (SNPCall));
  for (i = 0; i < nautos; i++) {
    unsigned int ntokenz;
    const unsigned char *tokenz[16];
    unsigned int lengths[16];
    unsigned int npairs, j;
    int best_a, best_b, best_delta;
    ntokenz = split_line (lines[indices[i]], lines[indices[i] + 1] - lines[indices[i]], tokenz, lengths, 8);
    if (ntokenz < 4) continue;
    npairs = (ntokenz - 2) / 2;
    best_delta = 0x7fffffff;
    for (j = 0; j < npairs; j++) {
      int a = strtol ((const char *) tokenz[2 + 2 * j], NULL, 10);
      int b = strtol ((const char *) tokenz[2 + 2 * j + 1], NULL, 10);
      int delta = (a + b) - (int) pair_median;
      if (delta < 0) delta = - delta;
      if (delta < best_delta) {
        best_a = a;
        best_b = b;
        best_delta = delta;
      }
    }
    calls[i].line = indices[i];
    calls[i].counts[0] = best_a;
    calls[i].counts[1] = best_b;
  }
  return calls;
}

/* Build list of unique random integers between 0...size */

unsigned int *
build_training_set (unsigned int set_size, unsigned int subset_size)
{
  unsigned int *train, i;
  train = (unsigned int *) malloc (set_size * sizeof (unsigned int));
  /* for (i = 0; i < nautos; i++) train[i] = autos[i]; */
  for (i = 0; i < set_size; i++) train[i] = i;
  for (i = 0; i < subset_size; i++) {
    unsigned int p;
    unsigned int t;
    p = (unsigned int) rand_long_long (0, set_size - 1);
    t = train[i];
    train[i] = train[p];
    train[p] = t;
  }
  return train;
}

float
calculate_allele_freq (SNPCall *calls, unsigned int set[], unsigned int set_size)
{
  double ppB, npB;
  unsigned int i;

  ppB = npB = 0;
  for (i = 0; i < set_size; i++) {
    unsigned int c0, c1;
    if (set) {
      c0 = calls[set[i]].counts[0];
      c1 = calls[set[i]].counts[1];
    } else {
      c0 = calls[i].counts[0];
      c1 = calls[i].counts[1];
    }
    if (c0 + c1) {
      ppB += (1.0f * c1) / (c0 + c1);
      npB += 1;
    }
  }
  if (npB) {
    return (float) (ppB / npB);
  } else {
    return 0;
  }
}

static void
train_model (SNPCall *calls, unsigned int ncalls, unsigned int max_training, unsigned int nruns, float v[], float *pB, unsigned int mul, unsigned int nthreads)
{
  unsigned int *train;
  unsigned int ntrain;
  double s0, s1, ppB, npB;
  double keskmine;
  float params[7], deltas[7];
  L3Data l3;
  unsigned int i;
  unsigned int chunk_size;

  /* Train model */
  if (debug) fprintf (stderr, "Building training set...");
  ntrain = MIN (ncalls, max_training);
  train = build_training_set (ncalls, ntrain);
  if (debug) fprintf (stderr, "done\n");

  if (debug) fprintf (stderr, "Calculating mean...");
  s0 = s1 = ppB = npB = 0;
  for (i = 0; i < ntrain; i++) {
    unsigned int c0 = calls[train[i]].counts[0];
    unsigned int c1 = calls[train[i]].counts[1];
    s0 += c0;
    s1 += c1;
    if (c0 + c1) {
      ppB += (1.0f * c1) / (c0 + c1);
      npB += 1;
    }
  }
  *pB = calculate_allele_freq (calls, train, ntrain);
  keskmine = (s0 + s1) / ntrain;
  if (debug) fprintf (stderr, "done\n");
  if (debug) {
    fprintf (stderr, "A %g B %g\n", s0, s1);
    fprintf (stderr, "Training size %u mean %.1f\n", ntrain, keskmine);
    fprintf (stderr, "pB %.3f\n", *pB);
  }
  if (keskmine == 0) {
    fprintf (stderr, "No calls in training sample, aborting model optimization\n");
    return;
  }
  if (pB == 0) {
    fprintf (stderr, "No B allele in training sample, aborting model optimization\n");
    return;
  }

  /* Genoomis mitteesineva k-meeri keskmine kohtamiste arv readide seas (keskmine vigade arv) - l_viga */
  /* OPTIMIZING_PARAM = logit (ACTUAL_PARAM) */
  params[0] = logit_clamped (v[0], MIN_P, 0.1f);
  /* N-täheliste genotüüpide tõenäosused */
  /* OPTIMIZING_PARAM = logit (ACTUAL_PARAM) */
  params[1] = logit_clamped (v[1], MIN_P, 1 - MIN_P); /* p0 */
  params[2] = logit_clamped (v[2], MIN_P, 1 - MIN_P); /* p1 */
  params[3] = logit_clamped (v[3], MIN_P, 1 - MIN_P); /* p2 */
  /* Keskmine katvus */
  /* OPTIMIZING_PARAM = log (ACTUAL_PARAM) */
  params[4] = logf (mul * keskmine);
  params[5] = v[5]; /* size */
  params[6] = logf (-v[6]); /* size2 */

  for (i = 0; i < 7; i++) deltas[i] = params[i] / 10;

  l3.params = params;
  l3.n_calls = ntrain;
  l3.pB = *pB;
  l3.var1 = malloc (ntrain * sizeof (float));
  l3.var2 = malloc (ntrain * sizeof (float));
  for (i = 0; i < ntrain; i++) {
    l3.var1[i] = calls[train[i]].counts[0];
    l3.var2[i] = calls[train[i]].counts[1];
  }

  chunk_size = (ntrain + nthreads - 1) / nthreads;
  if (chunk_size < 2000) chunk_size = 2000;
  nthreads = (ntrain + chunk_size - 1) / chunk_size;
  l3.n_threads = nthreads;
  l3.pool = pool;
  for (i = 0; i < l3.n_threads; i++) {
    l3.optims[i].l3 = &l3;
    l3.optims[i].first_call = i * chunk_size;
    l3.optims[i].n_calls = chunk_size;
    if (l3.optims[i].first_call + l3.optims[i].n_calls > l3.n_calls) {
      l3.optims[i].n_calls = l3.n_calls - l3.optims[i].first_call;
    }
    l3.optims[i].sum = 0;
  }

  downhill_simplex (7, params, deltas, 1e-6, nruns, 100, distanceL3, &l3);

  if (debug) {
    float dist = distanceL3 (7, params, &l3);
    print_params (params, "Best", stderr);
    fprintf (stderr, "Best distance %.6f\n", dist);
  }

  free (l3.var1);
  free (l3.var2);

  /* Kui mitu sellist k-meeri-lugemit näeme siis, kui antud k-meeri tegelikult genoomis ei esinenud... */
  v[0] = logit_1_clamped (params[0], MIN_P, 0.1f);
  /* Tõenäosus omada midagi 0-korduses, 1-s korduses, 2-s korduses (nagu autosoomis) */
  /* Tõenäosus nÃ¤ha midagi enam kui kahes korduses (p_lisa) leitakse valemiga 1-p_0-p_1-p_2 */
  v[1] = logit_1_clamped (params[1], MIN_P, 1 - MIN_P);
  v[2] = logit_1_clamped (params[2], MIN_P, 1 - MIN_P);
  v[3] = logit_1_clamped (params[3], MIN_P, 1 - MIN_P);
  /* Keskmine katvus (autosoomide korral, sugukromosoomide korral katvus/2) */
  v[4] = expf (params[4]);
  /* Ülehajuvusparameeter (Kui palju laiem on readide jaotus */
  v[5] = params[5];
  v[6] = -expf (params[6]);

  free (train);
}

typedef struct _PData PData;
struct _PData {
  double a[NUM_GENOTYPES];
  double sum;
  unsigned int best;
};
typedef struct _POptim POptim;
struct _POptim {
  SNPCall *calls;
  unsigned int first;
  unsigned int ncalls;
  float pB;
  float *params;
  PData *pdata;
};

static void
calc_run (void *data)
{
  unsigned int i;
  POptim *optim = (POptim *) data;
  /* fprintf (stderr, "Run %u-%u\n", optim->first, optim->first + optim->ncalls); */
  for (i = 0; i < optim->ncalls; i++) {
    unsigned int var1, var2, j;
    double best;
    var1 = optim->calls[i].counts[0];
    var2 = optim->calls[i].counts[1];
    genotype_probabilities (optim->pdata[i].a, optim->pB, var1, var2, optim->params[L_VIGA], optim->params[P_0], optim->params[P_1], optim->params[P_2], optim->params[LAMBDA], optim->params[SIZE], optim->params[SIZE2]);
    optim->pdata[i].sum = optim->pdata[i].a[0];
    optim->pdata[i].best = 0;
    best = optim->pdata[i].a[0];
    for (j = 1; j < NUM_GENOTYPES; j++) {
      optim->pdata[i].sum += optim->pdata[i].a[j];
      if (optim->pdata[i].a[j] > best) {
        optim->pdata[i].best = j;
        best = optim->pdata[i].a[j];
      }
    }
  }
}

static void
print_genotypes (const unsigned char *lines[], SNPCall *calls, unsigned int ncalls, float params[], float pB, unsigned int nalleles, float pc, unsigned int alt, unsigned int nthreads)
{
  POptim optim[32];
  unsigned int chunk_size = 5000;
  unsigned int pos = 0;
  unsigned int iter = 0;
  PData *pdata = (PData *) malloc (nthreads * chunk_size * sizeof (PData));

  memset (optim, 0, sizeof (optim));

  while (pos < ncalls) {
    unsigned int cpos = pos;
    unsigned int tidx, j, t;
    if (debug > 1) fprintf (stderr, "Printing iteration %u\n", iter++);
    for (tidx = 0; tidx < nthreads; tidx++) {
      unsigned int cend;
      if (cpos >= ncalls) break;
      optim[tidx].calls = calls + cpos;
      optim[tidx].first = cpos;
      optim[tidx].pdata = pdata + tidx * chunk_size;
      cend = cpos + chunk_size;
      if (cend >= ncalls) cend = ncalls;
      optim[tidx].ncalls = cend - cpos;
      optim[tidx].pB = pB;
      optim[tidx].params = params;
      aosora_thread_pool_submit (pool, NULL, calc_run, NULL, NULL, &optim[tidx]);
      cpos = cend;
    }
    aosora_thread_pool_run (pool);
    for (t = 0; t < tidx; t++) {
      unsigned int k;
      for (k = 0; k < optim[t].ncalls; k++) {
        double *a;
        double summa;
        unsigned int best_gt;
        char c[64];
        const unsigned char *p;
        unsigned int cancall;
        SNPCall *call = calls + optim[t].first + k;
        a = optim[t].pdata[k].a;
        summa = optim[t].pdata[k].sum;
        best_gt = optim[t].pdata[k].best;

        p = lines[call->line];
        j = 0;
        while (p[j] != '\t') j += 1;
        memcpy (c, p, j);
        c[j] = 0;
        fprintf (stdout, "%s", c);
        cancall = 0;
        if (nalleles == 0) {
          cancall = 1;
        } else if (nalleles == 1) {
          if ((best_gt == A) || (best_gt == B)) cancall = 1;
        } else if (nalleles == 2) {
          if ((best_gt == AA) || (best_gt == AB) || (best_gt == BB)) cancall = 1;
        }
        if (a[best_gt] < pc) cancall = 0;
        if (!optim[t].calls[k].counts[0] && !optim[t].calls[k].counts[1]) cancall = 0;
        if (cancall) {
          fprintf (stdout, "\t%s\t%.2f", gt[best_gt], a[best_gt] / summa);
        } else {
          fprintf (stdout, "\tNC\t");
        }
        fprintf (stdout, "\t%u\t%u", call->counts[0], call->counts[1]);
        if (alt) {
          for (j = 0; j < NUM_GENOTYPES; j++) {
            fprintf (stdout, "\t%.2f", a[j] / summa);
          }
        }
        fprintf (stdout, "\n");
      }
    }
    pos = cpos;
    if (debug > 1) fprintf (stderr, "Writing, pos %u\n", pos);
  }
  free (pdata);
}

static void
print_usage (FILE *ofs)
{
  fprintf (ofs, "Usage:\n");
  fprintf (ofs, "  gmer_caller ARGUMENTS COUNTS_FILE\n");
  fprintf (ofs, "Arguments:\n");
  fprintf (ofs, "    --training_size NUM - Use NUM markers for training (default 100000)\n");
  fprintf (ofs, "    --runs NUMBER       - Perfom NUMBER runs of model training (use 0 for no training)\n");
  fprintf (ofs, "    --num_threads NUM   - Use NUM threads (min 1, max %u, default %u)\n", MAX_THREADS, MAX_THREADS / 2);
  fprintf (ofs, "    --header            - Print table header\n");
  fprintf (ofs, "    --non_canonical     - Output non-canonical genotypes\n");
  fprintf (ofs, "    --prob_cutoff       - probability cutoff for calling genotype (default 0)\n");
  fprintf (ofs, "    --alternatives      - Print probabilities of all alternative genotypes\n");
  fprintf (ofs, "    --info              - Print information about individual\n");
  fprintf (ofs, "    --no_genotypes      - Print only summary information, not actual genotypes\n");
  fprintf (ofs, "    --model TYPE        - Model type (full, diploid, haploid)\n");
  fprintf (ofs, "    --params PARAMS     - Model parameters (error, p0, p1, p2, coverage, size, size2)\n");
  fprintf (ofs, "    --coverage NUM      - Average coverage of reads\n");
  fprintf (ofs, "    -D                  - increase debug level\n");
}

enum { MODEL_FULL, MODEL_DIPLOID, MODEL_HAPLOID };

int
main (int argc, const char *argv[])
{
  const char *call_fn = NULL;
  unsigned int nruns = 5;
  unsigned int max_training = 100000;
  unsigned int nthreads = MAX_THREADS / 2;
  unsigned int header = 0;
  unsigned int non_canonical = 0;
  float prob_cutoff = 0;
  unsigned int alternatives = 0;
  unsigned int info = 0;
  unsigned int print_gt = 1;
  unsigned int model = MODEL_FULL;
  unsigned int params_specified = 0;
  unsigned int coverage_specified = 0;

  unsigned int i;
  float pB;
  int aidx;
  float params[7], x_params[7];
  const unsigned char *cdata;
  unsigned long long csize;
  unsigned int nlines;
  const unsigned char **lines;
  unsigned int na, nx, ny;
  unsigned int *a, *x, *y;
  unsigned int a_med, x_med, y_med;
  double p_XX, p_X, p_Y, p_1;
  SNPCall *calls_a, *calls_x, *calls_y;

  /* Initial parameters (diploid) */
  params[L_VIGA] = 0.0547219f;
  params[P_0] = 4.2603e-05f;
  params[P_1] = 0.014934f;
  params[P_2] = 0.985023f;
  params[LAMBDA] = 0;
  params[SIZE] = 65.48f;
  params[SIZE2] = -0.6792684f;

  srand (1);

  /* Initialize combination table for multithreaded use */
  init_combination_tables ();

  aidx = 1;
  while (aidx < argc) {
    if (!strcmp (argv[aidx], "-D")) {
      debug += 1;
    } else if (!strcmp (argv[aidx], "--runs")) {
      aidx += 1;
      if (aidx >= argc) {
        print_usage (stderr);
        exit (1);
      }
      nruns = strtol (argv[aidx], NULL, 10);
    } else if (!strcmp (argv[aidx], "--training_size")) {
      aidx += 1;
      if (aidx >= argc) {
        print_usage (stderr);
        exit (1);
      }
      max_training = strtol (argv[aidx], NULL, 10);
    } else if (!strcmp (argv[aidx], "--num_threads")) {
      aidx += 1;
      if (aidx >= argc) {
        print_usage (stderr);
        exit (1);
      }
      nthreads = strtol (argv[aidx], NULL, 10);
    } else if (!strcmp (argv[aidx], "--header")) {
      header = 1;
    } else if (!strcmp (argv[aidx], "--non_canonical")) {
      non_canonical = 1;
    } else if (!strcmp (argv[aidx], "--prob_cutoff")) {
      aidx += 1;
      if (aidx >= argc) {
        print_usage (stderr);
        exit (1);
      }
      prob_cutoff = atof (argv[aidx]);
    } else if (!strcmp (argv[aidx], "--model")) {
      aidx += 1;
      if (aidx >= argc) {
        print_usage (stderr);
        exit (1);
      }
      if (!strcmp (argv[aidx], "full")) {
        model = MODEL_FULL;
      } else if (!strcmp (argv[aidx], "diploid")) {
        model = MODEL_DIPLOID;
      } else if (!strcmp (argv[aidx], "haploid")) {
        model = MODEL_HAPLOID;
      } else {
        print_usage (stderr);
        exit (1);
      }
      prob_cutoff = atof (argv[aidx]);
    } else if (!strcmp (argv[aidx], "--params")) {
      aidx += 1;
      if ((aidx + 6) >= argc) {
        print_usage (stderr);
        exit (1);
      }
      for (i = 0; i < 7; i++) {
        params[i] = atof (argv[aidx + i]);
      }
      params_specified = 1;
      aidx += 6;
    } else if (!strcmp (argv[aidx], "--coverage")) {
      aidx += 1;
      if (aidx >= argc) {
        print_usage (stderr);
        exit (1);
      }
      params[LAMBDA] = atof (argv[aidx]);
      coverage_specified = 1;
    } else if (!strcmp (argv[aidx], "--alternatives")) {
      alternatives = 1;
    } else if (!strcmp (argv[aidx], "--info")) {
      info = 1;
    } else if (!strcmp (argv[aidx], "--no_genotypes")) {
      print_gt = 0;
    } else {
      if (call_fn) {
        print_usage (stderr);
        exit (1);
      }
      call_fn = argv[aidx];
    }
    aidx += 1;
  }

  if (!call_fn) {
    fprintf (stderr, "No input file specified\n");
    print_usage (stderr);
  }
  if ((nthreads < 1) || (nthreads > MAX_THREADS)) {
    fprintf (stderr, "Invalid number of threads %u - should be 1-%u\n", nthreads, MAX_THREADS);
    print_usage (stderr);
  }

  /* If no user-specified params and haploid model change default probabilities */
  if ((model == MODEL_HAPLOID) && !params_specified) {
    params[P_1] = 0.985023f;
    params[P_2] = 0.014934f;
  }

  pool = aosora_thread_pool_new (nthreads);

  /* Read calls */
  if (debug) fprintf (stderr, "Reading %s...", call_fn);
  cdata = (const unsigned char *) gt4_mmap (call_fn, &csize);
  if (!cdata) {
    fprintf (stderr, "Cannot read %s\n", call_fn);
    exit (1);
  }
  nlines = count_lines (cdata, csize);
  if (nlines < 1) {
    fprintf (stderr, "File contains no lines\n");
    exit (1);
  }
  if (debug) fprintf (stderr, "done (%u lines)\n", nlines);
  if (debug) fprintf (stderr, "Building line table...");
  lines = (const unsigned char **) malloc ((nlines + 1) * sizeof (unsigned char *));
  nlines = parse_lines (lines, cdata, csize);
  lines[nlines] = cdata + csize;
  if (debug) fprintf (stderr, "done\n");

  /* Count chromosomes */
  /* Full counts are only needed for FULL model, otherwise treat all as autosomes */
  if (debug) fprintf (stderr, "Counting chromosomes...");
  na = nx = ny = 0;
  for (i = 0; i < nlines; i++) {
    if ((model != MODEL_FULL) || ((lines[i][0] > '0') && (lines[i][0] <= '9'))) {
      na += 1;
    } else if (lines[i][0] == 'X') {
      nx += 1;
    } else if (lines[i][0] == 'Y') {
      ny += 1;
    }
  }
  a = (unsigned int *) malloc (na * sizeof (unsigned int));
  x = (unsigned int *) malloc (nx * sizeof (unsigned int));
  y = (unsigned int *) malloc (ny * sizeof (unsigned int));
  na = nx = ny = 0;
  for (i = 0; i < nlines; i++) {
    if ((model != MODEL_FULL) || ((lines[i][0] > '0') && (lines[i][0] <= '9'))) {
      a[na++] = i;
    } else if (lines[i][0] == 'X') {
      x[nx++] = i;
    } else if (lines[i][0] == 'Y') {
      y[ny++] = i;
    }
  }
  if (debug) fprintf (stderr, "done\n");
  if (debug) fprintf (stderr, "Autosomes %u X %u Y %u\n", na, nx, ny);

  /* Pair medians */
  a_med = x_med = y_med = 0;
  if (debug) fprintf (stderr, "Calculating medians...");
  a_med = get_pair_median (cdata, csize, lines, a, na);
  if (model == MODEL_FULL) {
    x_med = get_pair_median (cdata, csize, lines, x, nx);
    y_med = get_pair_median (cdata, csize, lines, y, ny);
  }
  if (debug) fprintf (stderr, "done\n");
  if (debug) fprintf (stderr, "Autosomes/unspecified %u X %u Y %u\n", a_med, x_med, y_med);

  /* Determine sex */
  /* Only for full model */
  p_XX = p_X = p_Y = p_1 = 0;
  if (model == MODEL_FULL) {
    p_XX = poisson (x_med, a_med);
    p_X = poisson (x_med, a_med / 2);
    p_Y = poisson (y_med, a_med / 2);
    p_1 = poisson (y_med, 1);
    if (debug) {
      fprintf (stderr, "XX %g X %g Y %g 0 %g\n", p_XX, p_X, p_Y, p_1);
      if (p_XX > p_X) {
        fprintf (stderr, "Probably female\n");
      } else {
        fprintf (stderr, "Probably male\n");
      }
    }
    if (p_XX > p_X) {
      if (p_Y > p_1) {
        fprintf (stderr, "Y inconsistency: p_1 %g p_Y %g p_X %g p_XX %g\n", p_1, p_Y, p_X, p_XX);
      }
    } else {
      if (p_Y < p_1) {
        fprintf (stderr, "Y inconsistency: p_1 %g p_Y %g p_X %g p_XX %g\n", p_1, p_Y, p_X, p_XX);
      }
    }
  }

  /* Parse autosome/unspecified calls calls */
  if (debug) fprintf (stderr, "Reading autosome/unspecified calls...");
  calls_a = parse_calls (lines, a, na, a_med);
  if (debug) fprintf (stderr, "done\n");

  /* Train autosome model */
  /* Only if na > 0 && nruns > 0 */
  if (nruns && (na > 0)) {
    if (debug) fprintf (stderr, "Training autosome/unspecified model\n");
    if ((model == MODEL_FULL) || (model == MODEL_DIPLOID)) {
      /* Train full model */
      train_model (calls_a, na, max_training, nruns, params, &pB, 1, nthreads);
    } else if (model == MODEL_HAPLOID) {
      /* Haploid model */
      train_model (calls_a, na, max_training, nruns, params, &pB, 2, nthreads);
    }
  } else {
    /* Estimate MAF */
    pB = calculate_allele_freq (calls_a, NULL, na);
  }


  if (info) {
    if (model == MODEL_FULL) fprintf (stdout, "#Sex\t%s\n", (p_XX > p_X) ? "F" : "M");
    fprintf (stdout, "#EstimatedCoverage\t%g\n", params[LAMBDA]);
    fprintf (stdout, "#AverageMAF\t%g\n", pB);
    fprintf (stdout, "#AutosomeModel\t%g %g %g %g %g %g %g\n", params[L_VIGA], params[P_0], params[P_1], params[P_2], params[LAMBDA], params[SIZE], params[SIZE2]);
  }

  /* Default X model is diploid */
  memcpy (x_params, params, sizeof (params));
  /* If model is FULL, parse X calls */
  if (model == MODEL_FULL) {
    if (debug) fprintf (stderr, "Reading X calls...");
    calls_x = parse_calls (lines, x, nx, x_med);
    if (debug) fprintf (stderr, "done\n");
    /* Train separate haploid model */
    if ((nx > 0) && nruns && (p_XX <= p_X)) {
      /* Male */
      x_params[P_1] = 0.98f;
      x_params[P_2] = 0.01f;
      if (debug) fprintf (stderr, "Training X model\n");
      train_model (calls_x, nx, max_training, nruns, x_params, &pB, 2, nthreads);
      if (info) {
        fprintf (stdout, "#XModel\t%g %g %g %g %g %g %g\n", x_params[L_VIGA], x_params[P_0], x_params[P_1], x_params[P_2], x_params[LAMBDA], x_params[SIZE], x_params[SIZE2]);
      }
    }
  }

  /* Print genotypes */
  if (print_gt) {
    if (header) {
      fprintf (stdout, "#ID\tGT\tPROB\tA_KMERS\tB_KMERS");
      if (header) {
        unsigned int j;
        for (j = 0; j < NUM_GENOTYPES; j++) {
          fprintf (stdout, "\t%s", gt[j]);
        }
      }
      fprintf (stdout, "\n");
    }
    if (model != MODEL_HAPLOID) {
      print_genotypes (lines, calls_a, na, params, pB, (non_canonical) ? 0 : 2, prob_cutoff, alternatives, nthreads);
    } else {
      print_genotypes (lines, calls_a, na, params, pB, (non_canonical) ? 0 : 1, prob_cutoff, alternatives, nthreads);
    }
    if (model == MODEL_FULL) {
      if (p_XX > p_X) {
        /* Female */
        print_genotypes (lines, calls_x, nx, params, pB, (non_canonical) ? 0 : 2, prob_cutoff, alternatives, nthreads);
      } else {
        /* Male */
        print_genotypes (lines, calls_x, nx, x_params, pB, (non_canonical) ? 0 : 1, prob_cutoff, alternatives, nthreads);
        if (debug) fprintf (stderr, "Reading Y calls...");
        calls_y = parse_calls (lines, y, ny, y_med);
        if (debug) fprintf (stderr, "done\n");
        print_genotypes (lines, calls_y, ny, x_params, pB, (non_canonical) ? 0 : 1, prob_cutoff, alternatives, nthreads);
      }
    }
  }

  aosora_thread_pool_delete (pool);

  return 0;
}

static double
mlogL3 (float l_viga, float p_0, float p_1, float p_2, float lambda, float size, float size2, unsigned int n_calls, float pB, const unsigned int var1[], const unsigned int var2[])
{
  unsigned int i;
  double sum;
  sum = 0;
  for (i = 0; i < n_calls; i++) {
    double a[NUM_GENOTYPES];
    double abi;
    int j;

    genotype_probabilities (a, pB, var1[i], var2[i], l_viga, p_0, p_1, p_2, lambda, size, size2);

    abi = 0;
    for (j = 0; j < NUM_GENOTYPES; j++) {
      assert (!isnan (a[j]));
      abi += a[j];
    }
    assert (!isnan (abi));

    /* abi[abi<1e-30]=1e-30 */
    if (abi < 1e-30) abi = 1e-30;
    /* l = sum(log(abi))+3000000 */
    sum += log (abi);
  }

  return -sum;
}

static void
optim_run (void *data)
{
  L3Optim *optim = (L3Optim *) data;
  float l_viga, p_0, p_1, p_2, lambda, size, size2;
  l_viga = logit_1_clamped (optim->l3->params[0], MIN_P, 0.1f);
  p_0 = logit_1_clamped (optim->l3->params[1], MIN_P, 1 - MIN_P);
  p_1 = logit_1_clamped (optim->l3->params[2], MIN_P, 1 - MIN_P);
  p_2 = logit_1_clamped (optim->l3->params[3], MIN_P, 1 - MIN_P);
  lambda = expf (optim->l3->params[4]);
  size = optim->l3->params[5];
  size2 = -expf (optim->l3->params[6]);
  if (debug > 2) fprintf (stderr, "Optimization: %u-%u\n", optim->first_call, optim->first_call + optim->n_calls);
  optim->sum = mlogL3 (l_viga, p_0, p_1, p_2, lambda, size, size2, optim->n_calls, optim->l3->pB, optim->l3->var1 + optim->first_call, optim->l3->var2 + optim->first_call);
}

static float
distanceL3 (int ndim, const float params[], void *data)
{
  static int iter = 1;
  L3Data *l3;
  unsigned int i;
  float l_viga, p_0, p_1, p_2, lambda, size, size2;
  double delta0, delta1;
  double result;

  l3 = (L3Data *) data;

  if (debug > 1) print_params (params, "Params", stderr);

  l_viga = logit_1_clamped (l3->params[0], MIN_P, 0.1f);
  assert (!isnan (l_viga));
  p_0 = logit_1_clamped (l3->params[1], MIN_P, 1 - MIN_P);
  assert (!isnan (p_0));
  p_1 = logit_1_clamped (l3->params[2], MIN_P, 1 - MIN_P);
  assert (!isnan (p_1));
  p_2 = logit_1_clamped (l3->params[3], MIN_P, 1 - MIN_P);
  assert (!isnan (p_1));
  lambda = expf (l3->params[4]);
  assert (!isnan (lambda));
  size = l3->params[5];
  size2 = -expf (l3->params[6]);

  for (i = 0; i < l3->n_threads; i++) {
    l3->optims[i].sum = 0;
    aosora_thread_pool_submit (l3->pool, NULL, optim_run, NULL, NULL, &l3->optims[i]);
  }
  aosora_thread_pool_run (l3->pool);
  result = 0;
  for (i = 0; i < l3->n_threads; i++) {
    result += l3->optims[i].sum;
  }

  /* fixme: How to implement that? */
  if (p_0 + p_1 + p_2 > 1) {
    result = result + 10000 - 100000 * (1 - p_0 - p_1 - p_2);
  }
  delta0 = size + size2 * lambda / 2;
  if (delta0 < 0) {
    result = result + 10000 + 100 * delta0;
  }

  delta1 = size + size2 * l_viga;
  if (delta1 < 0) {
    result = result + 10000 + 100 * delta1;
  }

  if (debug > 1) fprintf (stderr, "Iteration %d delta %.6f\n", iter++, result);
  return (float) result;
}

static double
logL2 (const float arg[], unsigned int count, unsigned int tulem2[], unsigned int katvus2[])
{
  /* AAA, BBB */
  double p1 = arg[0];
  /* AAB, BBA */
  double p2 = arg[1];
  /* ABA, ABB, BAA, BAB */
  double p3 = p2;
  /* Fraction of other DNA */
  double x = arg[2];
#if 0
  if (x < 0) x = 0;
  if (x > 1) x = 1;
#endif
  /* Number of errors */
  double z = arg[3];
#if 0
  if (z < 0) z = 0;
  if (z > 1) z = 1;
#endif
  double logL;
  unsigned int i;

  logL = 0;
  for (i = 0; i < count; i++) {
    double v;

    v  = p1 * dbinom (tulem2[i], katvus2[i], z / 2);
    v += p1 * dbinom (tulem2[i], katvus2[i], 1 - z / 2);
    v += p2 * dbinom (tulem2[i], katvus2[i], (1 - x / 2) * (1 - z) + z / 2);
    v += p2 * dbinom (tulem2[i], katvus2[i], (x / 2) * (1 - z) + z / 2);
    v += 2 * p3 * dbinom (tulem2[i], katvus2[i], 0.5 * (1 - z) + z / 2);
    v += (0.5 - p1 - p2 - p3) * dbinom (tulem2[i], katvus2[i], (0.5 + x / 2) * (1 - z) + z / 2);
    v += (0.5 - p1 - p2 - p3) * dbinom (tulem2[i], katvus2[i], (0.5 - x / 2) * (1 - z) + z / 2);

    if (isnan (v)) {
      fprintf (stdout, "NaN %u %u\n", tulem2[i], katvus2[i]);
    }

    if (v <= 0) {
      /* fprintf (stdout, "Negative %u %u %g p1 = %g p2 = %g x = %g z = %g dbinom = %g\n", tulem2[i], katvus2[i], v, p1, p2, x, z, dbinom (tulem2[i], katvus2[i], 0.9)); */
      logL += 1000000;
    } else {
      logL -= log (v);
    }
  }
  return logL;
}

static unsigned int
get_pair_median (const unsigned char *cdata, unsigned long long csize, const unsigned char *lines[], unsigned int autos[], unsigned int nautos)
{
  unsigned int pair_median, ntokenz, i;
  const unsigned char *tokenz[16];
  unsigned int lengths[16];
  unsigned int max, min, med;
  unsigned int *medians_6 = (unsigned int *) malloc (nautos * sizeof (unsigned int));
  for (i = 0; i < nautos; i++) {
    unsigned int npairs, sum, j;
    ntokenz = split_line (lines[autos[i]], lines[autos[i] + 1] - lines[autos[i]], tokenz, lengths, 8);
    if (ntokenz < 4) continue;
    npairs = (ntokenz - 2) / 2;
    sum = 0;
    for (j = 0; j < npairs; j ++) {
      sum += strtol ((const char *) tokenz[2 + 2 * j], NULL, 10);
      sum += strtol ((const char *) tokenz[2 + 2 * j + 1], NULL, 10);
    }
    sum = sum * 6 / npairs;
    medians_6[i] = sum;
  }
  max = 0;
  min = 0xffffffff;
  for (i = 0; i < nautos; i++) {
    if (medians_6[i] > max) max = medians_6[i];
    if (medians_6[i] < min) min = medians_6[i];
  }
  med = (min + max) / 2;
  while (max > min) {
    unsigned int above = 0, below = 0, equal;
    for (i = 0; i < nautos; i++) {
      if (medians_6[i] > med) above += 1;
      if (medians_6[i] < med) below += 1;
    }
    equal = nautos - above - below;
    if (debug > 1) fprintf (stderr, "Trying median %u (%u) - equal %u, below %u, above %u\n", med / 6, med, equal, below, above);
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
  pair_median = med / 6;
  free (medians_6);

  return pair_median;
}
