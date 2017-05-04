#define __DISTRO_C__

#include <math.h>
#include <assert.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "binomial.h"
#include "queue.h"
#include "sequence.h"
#include "simplex.h"
#include "utils.h"
#include "wordmap.h"
#include "wordtable.h"

static unsigned int debug = 0;

typedef struct _Model Model;
typedef struct _List List;
typedef struct _Population Population;
typedef struct _Params Params;
typedef struct _KMer KMer;
typedef struct _Task Task;
typedef struct _DistroQueue DistroQueue;

static unsigned int read_population_data (Population *pop, const char *lfile, unsigned int max_lists, unsigned int min_cov);

static float find_coeffs (KMer *kmer);

static unsigned long long get_this_or_next (GT4WordMap *map, unsigned long long query);
static void print_coverages (unsigned int nlists, const char *lists[]);
static void get_distro (GT4WordMap *map, unsigned int *dist, unsigned int min, unsigned int max);
static unsigned int find_median_coverage (GT4WordMap *map);
static void generate_expected_distribution (float e[], const Params *params, unsigned int normalize);
static void generate_observed_distribution (float o[], const KMer *kmer, unsigned int haploid_coverage);
static void print_params (FILE *ofs, const Params *params, unsigned int header);
static void process_task (Queue *queue, unsigned int idx, Task *task);
static void process (Queue *queue, unsigned int idx, void *data);

static float default_coeffs[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

#define EXPECTED_KMER_COV_BITS 14

#define MAX_LISTS 1000
#define MAP_BLOCK 1000000
#define KMER_BLOCK_SIZE 4
#define WORDLENGTH 25
#define MIN_SIZE 4
#define MAX_SIZE 5000
#define MIN_ESIZE 0.05f
#define DEFAULT_NUM_THREADS 16

/*
 * Number of distinct possible coverages
 * The high distribution has coefficient 1 - SUM(coeffs)
 */
#define MAX_AMPLITUDE 6
#define N_AMPLITUDES (MAX_AMPLITUDE + 1)
enum {
  A0, A1, A2, A3, A4, A5, A6, P_MANY,
  P_SIZE,
  P_HCOV,
  P_ESIZE,
  N_PARAMS
};

/* Linear values */
enum {
  V0, V1, V2, V3, V4, V5, V6,
  V_SIZE,
  V_HCOV,
  V_ESIZE,
  N_VALUES
};

/*
 * Number of normalized k-mer counts in distribution
 * o[DISTRO_SIZE] is the sum of all higher counts
 */
#define MIN_COVERAGE 16
#define MAX_COVERAGE 40
#define HAPLOID_COVERAGE_NORM (MIN_COVERAGE / 2)
#define DIPLOID_COVERAGE_NORM (2 * HAPLOID_COVERAGE_NORM)
#define DISTRO_SIZE (MAX_AMPLITUDE * HAPLOID_COVERAGE_NORM)
#define ERROR_COVERAGE 0.1f

struct _Model {
  const char *id;
  /* Number of actual parameters */
  unsigned int n_params;
  /* Number of linear optimization values */
  unsigned int n_values;
  void (*p2v) (float d[], const float s[]);
  void (*v2p) (float d[], const float s[]);
  float (*distance) (int nvalues, const float values[], void *data);
};

static void initialize_models (Model models[]);

struct _List {
  const char *filename;
  unsigned int coverage;
  float x_coverage;
  float y_coverage;
  GT4WordMap *map;
  unsigned long long lptr;
};

struct _Population {
  unsigned long long current_kmer;
  unsigned int nlists;
  List *lists;
  unsigned int min_coverage;
  unsigned int max_coverage;
  float f_fraction;
  float m_fraction;
};

struct _Params {
  float values[N_PARAMS];
};

struct _KMer {
  Task *task;
  unsigned long long value;
  unsigned int obs[4096];
  Params params;
  float result;
};

struct _Task {
  DistroQueue *distro;
  unsigned int idx;
  unsigned int nkmers;
  KMer kmers[KMER_BLOCK_SIZE];
  Task *next;
  /* Temporaries for distance calculation */
  /* Last negative binomial distribution size */
  float last_size;
  /* Floating point log combination */
  float combinations[N_AMPLITUDES * MAX_COVERAGE];
  float last_esize;
  float e_logc[N_AMPLITUDES * MAX_COVERAGE];
  /* Precalculated per-coverage remainders */
  float sums[MAX_COVERAGE + 1];
};

struct _DistroQueue {
  Queue queue;
  Task *free_tasks;
  Task *allocated_tasks;
  Task *completed_tasks;
  Population *population;
  unsigned int processing;
  unsigned int finished;
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

/* Expand linear parameters to N_PARAMS actual parameters */
/* C0_logit, C1_logit, C2_logit, C3_logit, C4_logit, C5_logit, SIZE_log, HCOV, ESIZE */
static void
lin2params (float d[], const float s[])
{
  unsigned int i;
  float rem = 1;
  for (i = 0; i < N_AMPLITUDES; i++) {
    float q = logit_1 (s[i]) * rem;
    d[i] = q;
    assert (!isnan (d[i]));
    rem -= d[i];
    if (rem < 1e-6f) rem = 1e-6f;
  }
  d[P_MANY] = rem;
  d[P_SIZE] = logit_1 (s[V_SIZE]) * (MAX_SIZE - MIN_SIZE) + MIN_SIZE;
  assert (!isnan (d[P_SIZE]));
  d[P_HCOV] = 0.5f + logit_1 (s[V_HCOV]);
  assert (!isnan (d[P_HCOV]));
  d[P_ESIZE] = logit_1 (s[V_ESIZE]) * (MAX_SIZE - MIN_ESIZE) + MIN_ESIZE;
  assert (!isnan (d[P_ESIZE]));
}

/* Expand actual parameters to N_PARAMS - 1 linear parameters */
/* C0_logit, C1_logit, C2_logit, C3_logit, C4_logit, C5_logit, SIZE_log, HCOV, ESIZE */
static void
params2lin (float d[], const float s[])
{
  unsigned int i;
  float q, rem = 1;
  /* fixme: expand by sum */
  for (i = 0; i < N_AMPLITUDES; i++) {
    q = s[i];
    q /= rem;
    if (q < 1e-6f) q = 1e-6f;
    if (q > (1 - 1e-6f)) q = 1 - 1e-6f;
    d[i] = logit (q);
    assert (!isnan (d[i]));
    rem -= q;
    if (rem < 1e-6) rem = 1e-6;
  }
  q = (s[P_SIZE] - MIN_SIZE) / (MAX_SIZE - MIN_SIZE);
  if (q > (1 - 1e-6f)) q = 1 - 1e-6f;
  d[V_SIZE] = logit (q);
  assert (!isnan (d[V_SIZE]));
  d[V_HCOV] = logit (s[P_HCOV] - 0.5f);
  assert (!isnan (d[V_HCOV]));
  q = (s[P_ESIZE] - MIN_ESIZE) / (MAX_SIZE - MIN_ESIZE);
  if (q > (1 - 1e-6f)) q = 1 - 1e-6f;
  d[V_ESIZE] = logit (q);
  assert (!isnan (d[V_ESIZE]));
}

static double
rand_d (double min, double max)
{
  return min + (max - min) * rand () / (RAND_MAX + 1.0);
}

Model models[2];

int
main (int argc, const char *argv[])
{
  unsigned int argi, i, kmer_idx;

  unsigned long long current = 0;
  unsigned int max_lists = 10000;
  unsigned int max_kmers = 1000000000;
  unsigned int find_coverages = 0, print_oe = 0, print_distributions = 0;
  const char *lfile = NULL;
  Population pop;
  unsigned int nthreads = DEFAULT_NUM_THREADS;
  DistroQueue dq;
  Task *tasks;
  Task *current_task = NULL;

  initialize_models (models);

  argi = 1;
  while (argi < argc) {
    if (!strcmp (argv[argi], "--kmer")) {
      argi += 1;
      if (argi >= argc) exit (1);
      current = string_to_word (argv[argi], strlen (argv[argi]));
      max_kmers = 1;
    } else if (!strcmp (argv[argi], "--max_kmers")) {
      argi += 1;
      if (argi >= argc) exit (1);
      max_kmers = atoi (argv[argi]);
    } else if (!strcmp (argv[argi], "--max_lists")) {
      argi += 1;
      if (argi >= argc) exit (1);
      max_lists = atoi (argv[argi]);
    } else if (!strcmp (argv[argi], "--num_threads")) {
      argi += 1;
      if (argi >= argc) exit (1);
      nthreads = atoi (argv[argi]);
    } else if (!strcmp (argv[argi], "--print_oe")) {
      print_oe = 1;
    } else if (!strcmp (argv[argi], "--print_distributions")) {
      print_distributions = 1;
    } else if (!strcmp (argv[argi], "--find_coverages")) {
      find_coverages = 1;
    } else if (!strcmp (argv[argi], "-D")) {
      debug += 1;
    } else if (!strcmp (argv[argi], "--nbinom")) {
      fprintf (stdout, "%g\n", dnbinom_mu (atof (argv[argi + 1]), atof (argv[argi + 2]), atof (argv[argi + 3])));
      exit (0);
    } else {
      lfile = argv[argi];
    }
    argi += 1;
  }

  if (find_coverages) {
    unsigned int d[50];
    read_population_data (&pop, lfile, max_lists, 0);
    if (debug) fprintf (stderr, "Number of usable lists: %u\n", pop.nlists);
    for (i = 0; i < pop.nlists; i++) {
      unsigned int j;
      fprintf (stdout, "%s", pop.lists[i].filename);
      GT4WordMap *map = gt4_wordmap_new (pop.lists[i].filename, 0);
      get_distro (map, d, 2, 51);
      j = 1;
      while (j < 50) {
        if (d[j - 1] < d[j]) break;
        j += 1;
      }
      unsigned int max = j;
      while (j < 50) {
        if (d[j] > d[max]) max = j;
        j += 1;
      }
      fprintf (stdout, "\t%u\n", max + 2);
      if (debug) fprintf (stderr, "%u %s %u\n", i, pop.lists[i].filename, max + 2);
      gt4_wordmap_delete (map);
    }
    exit (0);
  }

  if (nthreads < 1) nthreads = 1;
  if (nthreads > 32) nthreads = 32;

  /* Initialize static tables */
  log_combinations_d (4, 2);
  log_combinations_f (4, 2);
  log_combination_k_r_1 (4, 2);
  log_combination_k_r_f (4, 2);
  
  /* Population */
  read_population_data (&pop, lfile, max_lists, MIN_COVERAGE);
  if (debug) {
    fprintf (stderr, "Number of usable lists: %u\n", pop.nlists);
    fprintf (stderr, "Females %.2f Males %.2f\n", pop.f_fraction, pop.m_fraction);
  }

  for (i = 0; i < pop.nlists; i++) {
    pop.lists[i].map = gt4_wordmap_new (pop.lists[i].filename, 0);
  }

#if 1
  /* Set list pointers */
  if (debug) fprintf (stderr, "Searching starting k-mer\n");
  if (!current) current = WORDMAP_WORD (pop.lists[0].map, (unsigned long long) (0.25 * pop.lists[0].map->header->nwords));
  for (i = 0; i < pop.nlists; i++) {
    pop.lists[i].lptr = get_this_or_next (pop.lists[i].map, current);
    if (debug) fprintf (stderr, ".");
  }
  if (debug) fprintf (stderr, " done\n");
#endif

  /* Queue */
  memset (&dq, 0, sizeof (DistroQueue));
  queue_init (&dq.queue, nthreads);
  dq.population = &pop;
  tasks = (Task *) malloc (256 * sizeof (Task));
  memset (tasks, 0, 256 * sizeof (Task));
  for (i = 0; i < 256; i++) {
    tasks[i].distro = &dq;
    tasks[i].next = dq.free_tasks;
    dq.free_tasks = &tasks[i];
  }
  queue_create_threads (&dq.queue, process, &dq);

  unsigned int job_idx = 0;
  kmer_idx = 0;
  unsigned int reading_finished = 0;
  unsigned int finished = 0;
  unsigned int next_to_print = 0;
  /* fixme: This is wrong */
  while (!finished) {
    unsigned int n_lists, n_observed;

    if (!reading_finished) {
      if (!current_task) {
        /* Pick new task */
        queue_lock (&dq.queue);
        if (!dq.free_tasks && !dq.completed_tasks && dq.allocated_tasks) {
          Task *task = dq.allocated_tasks;
          dq.allocated_tasks = task->next;
          queue_unlock (&dq.queue);
          process_task (&dq.queue, 0, task);
          queue_lock (&dq.queue);
          task->next = dq.completed_tasks;
          dq.completed_tasks = task;
        }
        if (dq.free_tasks) {
          current_task = dq.free_tasks;
          dq.free_tasks = current_task->next;
          /* Initialize task */
          memset (current_task, 0, sizeof (Task));
          current_task->distro = &dq;
          current_task->idx = job_idx++;
        }
        queue_unlock (&dq.queue);
      }

      if (current_task) {
      KMer *kmer = &current_task->kmers[current_task->nkmers];
      kmer->task = current_task;
      kmer->value = current;

      /* Calculate counts */
      memset (kmer->obs, 0, sizeof (kmer->obs));
      n_lists = 0;
      n_observed = 0;
      for (i = 0; i < pop.nlists; i++) {
        if (pop.lists[i].lptr < pop.lists[i].map->header->nwords) {
          n_lists += 1;
          if (WORDMAP_WORD (pop.lists[i].map, pop.lists[i].lptr) == current) {
            unsigned int freq = WORDMAP_FREQ (pop.lists[i].map, pop.lists[i].lptr);
            current_task->kmers[current_task->nkmers].obs[i] = freq;
            pop.lists[i].lptr += 1;
            n_observed += 1;
          }
        }
      }

      if (n_observed > 10) {
        current_task->nkmers += 1;
        kmer_idx += 1;
      }
      
      if (!n_lists || (kmer_idx >= max_kmers)) {
        reading_finished = 1;
        if (debug) fprintf (stderr, "Reading finished\n");
      } else {
        /* Pick new current */
        current = 0xffffffffffffffff;
        for (i = 0; i < pop.nlists; i++) {
          if ((pop.lists[i].lptr < pop.lists[i].map->header->nwords) && (WORDMAP_WORD (pop.lists[i].map, pop.lists[i].lptr) < current)) {
            current = WORDMAP_WORD (pop.lists[i].map, pop.lists[i].lptr);
          }
        }
      }

      if ((current_task->nkmers >= KMER_BLOCK_SIZE) || reading_finished) {
        /* Submit */
        if (debug) fprintf (stderr, "Submitting task %u\n", current_task->idx);
        queue_lock (&dq.queue);
        if (!dq.allocated_tasks) {
          current_task->next = dq.allocated_tasks;
          dq.allocated_tasks = current_task;
        } else {
          Task *t = dq.allocated_tasks;
          while (t->next) t = t->next;
          t->next = current_task;
        }
        current_task = NULL;
        queue_broadcast (&dq.queue);
        queue_unlock (&dq.queue);
      }
    }
    }

    Task *print_task = NULL;
    queue_lock (&dq.queue);
    if (dq.completed_tasks) {
      Task *p, *t;
      p = NULL;
      t = dq.completed_tasks;
      while (t) {
        if (t->idx == next_to_print) {
          if (p) {
            p->next = t->next;
          } else {
            dq.completed_tasks = t->next;
          }
          print_task = t;
          break;
        }
        p = t;
        t = t->next;
      }
    }
    queue_unlock (&dq.queue);

    if (print_task) {
      if (debug) fprintf (stderr, "Printing task %u with %u kmers\n", print_task->idx, print_task->nkmers);
      for (i = 0; i < print_task->nkmers; i++) {
        KMer *kmer = &print_task->kmers[i];
        char c[33];

        word2string (c, kmer->value, WORDLENGTH);
        fprintf (stdout, "KMer:\t%s\n", c);

        /* Print results */
        if (debug) fprintf (stderr, "Result = %g\n", kmer->result);
        print_params (stdout, &kmer->params, 0);
        if (print_oe) {
          float e[DISTRO_SIZE + 1];
          unsigned int j;
          /* Print ruler */
          fprintf (stdout, "R: ");
          for (j = 0; j <= DISTRO_SIZE; j++) {
            if ((fmod (j, kmer->params.values[P_HCOV] * HAPLOID_COVERAGE_NORM)) < 1) {
              fprintf (stdout, "_|_");
            } else {
              fprintf (stdout, "___");
            }
          }
          fprintf (stdout, "\n");
          /* Print observed */
          generate_observed_distribution (e, kmer, HAPLOID_COVERAGE_NORM);
          fprintf (stdout, "O:");
          for (j = 0; j <= DISTRO_SIZE; j++) {
            fprintf (stdout, " %2u", (unsigned int) (e[j] * kmer->task->distro->population->nlists + 0.5f));
          }
          fprintf (stdout, "\n");
          /* Print expected */
          generate_expected_distribution (e, &kmer->params, 1);
          fprintf (stdout, "E:");
          for (j = 0; j <= DISTRO_SIZE; j++) {
            fprintf (stdout, " %2u", (unsigned int) (e[j] * kmer->task->distro->population->nlists + 0.5f));
          }
          fprintf (stdout, "\n");
          if (print_distributions) {
            Params p;
            unsigned int j, k, l;
            for (j = 0; j < N_AMPLITUDES; j++) {
              p = kmer->params;
              for (k = 0; k < N_AMPLITUDES; k++) if (k != j) p.values[k] = 0;
              generate_expected_distribution (e, &p, 0);
              fprintf (stdout, "%u:", j);
              for (l = 0; l <= DISTRO_SIZE; l++) {
                fprintf (stdout, " %2u", (unsigned int) (e[l] * kmer->task->distro->population->nlists + 0.5f));
              }
              fprintf (stdout, "\n");
            }
          }
        }
      }
      queue_lock (&dq.queue);
      next_to_print = print_task->idx + 1;
      print_task->next = dq.free_tasks;
      dq.free_tasks = print_task;
      queue_broadcast (&dq.queue);
      queue_unlock (&dq.queue);
    }

    queue_lock (&dq.queue);
    if (reading_finished && !dq.allocated_tasks && !dq.completed_tasks && !dq.processing) {
      dq.finished = 1;
      if (dq.queue.nthreads_running < 2) finished = 1;
      queue_broadcast (&dq.queue);
    }
    queue_unlock (&dq.queue);
  }

  queue_finalize (&dq.queue);
  
  return 1;
}

/* Generates expected distrbution 0..DISTRO_SIZE, last element is remainder (sum == 1) */

static void
generate_expected_distribution (float e[], const Params *params, unsigned int normalize)
{
  unsigned int i, j;
  double sum = 0;
  memset (e, 0, (DISTRO_SIZE + 1) * sizeof (float));
  for (j = 0; j <= DISTRO_SIZE; j++) {
    double p = params->values[0] * dnbinom_mu (j, params->values[P_ESIZE], ERROR_COVERAGE);
    e[(j != 1) ? j : 0] += (float) p;
    if (j < DISTRO_SIZE) sum += p;
  }
  for (i = 1; i < N_AMPLITUDES; i++) {
    for (j = 0; j <= DISTRO_SIZE; j++) {
      double lambda = i * params->values[P_HCOV] * HAPLOID_COVERAGE_NORM;
      double p = params->values[i] * dnbinom_mu (j, params->values[P_SIZE], lambda);
      e[(j != 1) ? j : 0] += (float) p;
      if (j < DISTRO_SIZE) sum += p;
    }
  }
  if (normalize) e[DISTRO_SIZE] = 1 - sum;
}

/* Generates normalized observed distribution 0..DISTRO_SIZE, last element is remainder (sum == 1) */

static void
generate_observed_distribution (float o[], const KMer *kmer, unsigned int haploid_coverage)
{
  unsigned int i, j;
  memset (o, 0, (DISTRO_SIZE + 1) * sizeof (float));
  /* Calculate counts */
  for (i = 0; i < kmer->task->distro->population->nlists; i++) {
    if (debug > 3) fprintf (stderr, "List %s\n", kmer->task->distro->population->lists[i].filename);
    double p = (double) 2 * haploid_coverage / kmer->task->distro->population->lists[i].coverage;
    if (debug > 3) fprintf (stderr, "p = %g\n", p);
    unsigned int freq = kmer->obs[i];
    if (debug > 3) fprintf (stderr, "freq = %u\n", freq);
    /* Normalized freq */
    unsigned int fnorm = 0;
    if (freq >= 40) {
      fnorm = (unsigned int) (freq * p + 0.5);
    } else {
      for (j = 0; j < freq; j++) if (rand_d (0, 1) < p) fnorm += 1;
    }
    if (fnorm == 1) fnorm = 0;
    if (debug > 3) fprintf (stderr, "fnorm = %u\n", fnorm);
    if (fnorm < DISTRO_SIZE) {
      o[fnorm] += 1;
    } else {
      o[DISTRO_SIZE] += 1;
    }
  }
  for (i = 0; i <= DISTRO_SIZE; i++) o[i] /= kmer->task->distro->population->nlists;
}

static float
find_coeffs_1 (KMer *kmer, float haploid_coverage)
{
  float a[N_PARAMS - 1], d[N_PARAMS - 1];
  unsigned int i;
  float sum, res;

  /* Initialize coefficients */
  sum = 0;
  for (i = 0; i <= N_AMPLITUDES; i++) sum += default_coeffs[i];
  for (i = 0; i <= N_AMPLITUDES; i++) kmer->params.values[i] = default_coeffs[i] / sum;
  kmer->params.values[P_HCOV] = haploid_coverage;
  kmer->params.values[P_SIZE] = 100;
  kmer->params.values[P_ESIZE] = 1;

  models[1].p2v (a, kmer->params.values);
  for (i = 0; i < models[1].n_values; i++) d[i] = a[i] / 40;
  res = downhill_simplex (models[1].n_values, a, d, 1e-5f, 2, 200, models[1].distance, (void *) kmer);
  models[1].v2p (kmer->params.values, a);
  
  return res - 0.9f * (kmer->params.values[0] + kmer->params.values[2] + kmer->params.values[4] - kmer->params.values[1] - kmer->params.values[3] - kmer->params.values[5]);;
}

#define SPAN 12

static float
find_coeffs (KMer *kmer)
{
  float a[N_PARAMS - 1], d[N_PARAMS - 1];
  float result, best;
  Params params;
  int i;

  result = find_coeffs_1 (kmer, 1);
  best = result;
  params = kmer->params;
  if (debug > 1) fprintf (stderr, "C %.1f %g", (float) HAPLOID_COVERAGE_NORM, result);
  for (i = 1; i < SPAN; i++) {
    float cov = (HAPLOID_COVERAGE_NORM - i / 4.0f) / HAPLOID_COVERAGE_NORM;
    result = find_coeffs_1 (kmer, cov);
    if (debug > 1) fprintf (stderr, " - %.2f %g", cov * HAPLOID_COVERAGE_NORM, result);
    if (result < best) {
      best = result;
      params = kmer->params;
    } else if (result > best) {
      break;
    }
  }
  for (i = 1; i < SPAN; i++) {
    float cov = (HAPLOID_COVERAGE_NORM + i / 4.0f) / HAPLOID_COVERAGE_NORM;
    result = find_coeffs_1 (kmer, cov);
    if (debug > 1) fprintf (stderr, " - %.1f %g", cov * HAPLOID_COVERAGE_NORM, result);
    if (result < best) {
      best = result;
      params = kmer->params;
    } else if (result > best) {
      break;
    }
  }
  if (debug > 1) fprintf (stderr, "\n");
  
  kmer->params = params;

  models[0].p2v (a, kmer->params.values);
  for (i = 0; i < models[0].n_values; i++) d[i] = a[i] / 40;
  kmer->result = downhill_simplex (models[0].n_values, a, d, 1e-5f, 4, 400, models[0].distance, (void *) kmer);
  models[0].v2p (kmer->params.values, a);
  
  return kmer->result;
}

/* Models */

enum {
  MODEL_ORIG_FULL,
  MODEL_ORIG_AMPL
};

static float
distance (const float a[], float h_cov, float size, float e_size, void *data)
{
  KMer *kmer = (KMer *) data;
  Task *task = kmer->task;
  Population *pop = task->distro->population;
  unsigned int i;
  float D = 0;

  /* Precalculate coefficients */
  if (size != task->last_size) {
    for (i = 0; i < (MAX_AMPLITUDE / 2) * pop->max_coverage; i++) {
      task->combinations[i] = log_combination_k_r_f (i, size);
    }
    task->last_size = size;
  }
  if (e_size != task->last_esize) {
    for (i = 0; i < (MAX_AMPLITUDE / 2) * pop->max_coverage; i++) {
      task->e_logc[i] = log_combination_k_r_f (i, e_size);
    }
    task->last_esize = e_size;
  }
  memset (task->sums, 0, sizeof (task->sums));
  
  for (i = 0; i < pop->nlists; i++) {
    unsigned int j;
    float s;
    if (kmer->obs[i] < 2) {
      s = 0;
      for (j = 0; j < 2; j++) {
        s += a[0] * dnbinom_mu_precalc_f (j, e_size, ERROR_COVERAGE, task->e_logc[j]);
        assert (!isnan (s));
      }
      for (j = 1; j < N_AMPLITUDES; j++) {
        float lambda = j * h_cov * pop->lists[i].coverage / 2;
        s += a[j] * dnbinom_mu_precalc_f (0, size, lambda, task->combinations[0]);
        assert (!isnan (s));
        s += a[j] * dnbinom_mu_precalc_f (1, size, lambda, task->combinations[1]);
        assert (!isnan (s));
      }
    } else if (kmer->obs[i] < ((MAX_AMPLITUDE / 2) * pop->lists[i].coverage)) {
      s = 0;
      s += a[0] * dnbinom_mu_precalc_f (kmer->obs[i], e_size, ERROR_COVERAGE, task->e_logc[kmer->obs[i]]);
      assert (!isnan (s));
      for (j = 1; j < N_AMPLITUDES; j++) {
        float lambda = j * h_cov * pop->lists[i].coverage / 2;
        s += a[j] * dnbinom_mu_precalc_f (kmer->obs[i], size, lambda, task->combinations[kmer->obs[i]]);
        assert (!isnan (s));
      }
    } else {
      if (!task->sums[pop->lists[i].coverage]) {
        unsigned int o;
        s = 1;
        for (o = 0; o < ((MAX_AMPLITUDE / 2) * pop->lists[i].coverage); o++) {
          s -= a[0] * dnbinom_mu_precalc_f (o, e_size, ERROR_COVERAGE, task->e_logc[o]);
          assert (!isnan (s));
        }
        for (j = 1; j < N_AMPLITUDES; j++) {
          float lambda = j * h_cov * pop->lists[i].coverage / 2;
          for (o = 0; o < ((MAX_AMPLITUDE / 2) * pop->lists[i].coverage); o++) {
            s -= a[j] * dnbinom_mu_precalc_f (o, size, lambda, task->combinations[o]);
            assert (!isnan (s));
          }
        }
        if (s < 1e-8f) s = 1e-8f;
        task->sums[pop->lists[i].coverage] = s;
      } else {
        s = task->sums[pop->lists[i].coverage];
      }
    }
    if (s < 1e-8f) s = 1e-8f;
    D -= logf (s);
    assert (!isnan (D));
  }
  float x = a[P_HCOV];
  float s = 0.37f;
  D -= (-(x * x) / (2 * s * s));
  D /= pop->nlists;
  if (D < -1e20) D = -1e20;
  if (D > 1e20) D = 1e20;
  return D;
}

/* Original full model */

static float
distance_full (int nvalues, const float values[], void *data)
{
  float a[32], D;
  lin2params (a, values);
  D = distance (a, a[P_HCOV], a[P_SIZE], a[P_ESIZE], data);
  return D;
}

/* Original amplitude-only model */

static float
distance_partial (int nvalues, const float values[], void *data)
{
  float a[32], D;
  lin2params (a, values);
  D = distance (a, a[P_HCOV], a[P_SIZE], a[P_ESIZE], data);
  return D;
}

static void
initialize_models (Model models[])
{
  /* Full model */
  models[0].id = "Full";
  models[0].n_params = N_PARAMS;
  models[0].n_values = N_VALUES;
  models[0].p2v = params2lin;
  models[0].v2p = lin2params;
  models[0].distance = distance_full;
  /* Partial amplitude model */
  models[1].id = "Amplitudes";
  models[1].n_params = N_PARAMS;
  models[1].n_values = N_AMPLITUDES;
  models[1].p2v = params2lin;
  models[1].v2p = lin2params;
  models[1].distance = distance_partial;
}

static float
distance_by_model (Model *model, void *data)
{
}

static float
distance_weighted (int nvalues, const float values[], void *data)
{
  float a[N_PARAMS];
  lin2params (a, values);
  float result = distance (values, a[P_HCOV], a[P_SIZE], a[P_ESIZE], data);
  return result - 0.9f * (a[0] + a[2] + a[4] + a[6] - a[1] - a[3] - a[5]);
}

static void
process_task (Queue *queue, unsigned int idx, Task *task)
{
  unsigned int i;
  if (debug) fprintf (stderr, "Thread %u analyzing task %u\n", idx, task->idx);
  for (i = 0; i < task->nkmers; i++) {
    if (debug > 1) fprintf (stderr, "Thread %d: Analyzing kmer %u\n", idx, i);
    task->kmers[i].result = find_coeffs (&task->kmers[i]);
  }
}

static void
process (Queue *queue, unsigned int idx, void *data)
{
  DistroQueue *dq = (DistroQueue *) data;
  unsigned int finished = 0;

  if (debug > 1) {
    queue_lock (&dq->queue);
    fprintf (stderr, "Thread %d started (total %d)\n", idx, dq->queue.nthreads_running);
    queue_unlock (&dq->queue);
  }

  while (!finished) {
    queue_lock (&dq->queue);
    if (dq->allocated_tasks) {
      dq->processing += 1;
      /* Pick task */
      Task *task = dq->allocated_tasks;
      dq->allocated_tasks = task->next;
      queue_unlock (&dq->queue);
      process_task (queue, idx, task);
      queue_lock (&dq->queue);
      task->next = dq->completed_tasks;
      dq->completed_tasks = task;
      dq->processing -= 1;
      queue_broadcast (&dq->queue);
      queue_unlock (&dq->queue);
    } else if (!dq->finished) {
      /* Wait for main thread */
      if (debug > 1) fprintf (stderr, "Thread %d: Waiting\n", idx);
      queue_wait (&dq->queue);
      if (debug > 1) fprintf (stderr, "Thread %d: Woke up\n", idx);
      queue_unlock (&dq->queue);
    } else {
      /* Finished */
      queue_unlock (&dq->queue);
      finished = 1;
    }
  }
}

static void
print_params (FILE *ofs, const Params *params, unsigned int header)
{
  unsigned int j;
  if (header) {
    fprintf (ofs, "A0\tA1\tA2\tA3\tA4\tA5\tA+\tSIZE\tHCOV\tESIZE\n");
  }
  for (j = 0; j < N_PARAMS; j++) {
    if (j > 0) fprintf (ofs, "\t");
    fprintf (ofs, "%.4f", params->values[j]);
  }
  fprintf (ofs, "\n");
}

static unsigned long long
get_this_or_next (GT4WordMap *map, unsigned long long query)
{
  unsigned long long low, high, mid, word;
  low = 0;
  high = map->header->nwords - 1;
  while (low <= high) {
    mid = (low + high) / 2;
    word = WORDMAP_WORD (map, mid);
    if (word < query) {
      low = mid + 1;
    } else if (word > query) {
      if (mid == low) return high;
      high = mid;
    } else {
      return mid;
    }
  }
  return mid;
}

static void
print_coverages (unsigned int nlists, const char *lists[])
{
  unsigned int i;
  for (i = 1; i < nlists; i++) {
    GT4WordMap *map;
    unsigned int j;
    unsigned int dist[30];
    unsigned int max;
    
    map = gt4_wordmap_new (lists[i], 1);
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
    gt4_wordmap_delete (map);
  }
}

static void
get_distro (GT4WordMap *map, unsigned int *dist, unsigned int min, unsigned int max)
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
find_median_coverage (GT4WordMap *map)
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

static unsigned int
read_population_data (Population *pop, const char *lfile, unsigned int max_lists, unsigned int min_cov)
{
  const unsigned char *cdata;
  unsigned long long csize, cpos;
  unsigned int num_f = 0, num_m = 0;

  memset (pop, 0, sizeof (Population));
  pop->lists = (List *) malloc (max_lists * sizeof (List));
  pop->nlists = 0;

  cdata = (const unsigned char *) gt4_mmap (lfile, &csize);
  if (!cdata) return 0;
  cpos = 0;
  while (cpos < csize) {
    unsigned long long s, e;
    char b[1024];
    unsigned int size = 0;
    float xcov = 1;
    float ycov = 1;
    /* Skip initial whitespace */
    s = cpos;
    while ((cdata[s] <= ' ') && (s < csize)) s += 1;
    if (s >= csize) break;
    /* Get filename */
    e = s;
    while ((cdata[e] > ' ') && (e < csize)) e += 1;
    memcpy (b, cdata + s, e - s);
    b[e - s] = 0;
    /* Skip whitespace */
    s = e;
    while ((cdata[s] <= ' ') && (cdata[s] != '\n') && (s < csize)) s += 1;
    if (s >= csize) break;
    /* Read coverage */
    e = s;
    if (cdata[s] != '\n') {
      size = strtol ((const char *) cdata + s, NULL, 10);
      while ((cdata[e] > ' ') && (e < csize)) e += 1;
    }
    /* Skip whitespace */
    s = e;
    while ((cdata[s] <= ' ') && (cdata[s] != '\n') && (s < csize)) s += 1;
    if (s >= csize) break;
    /* Read X coverage */
    e = s;
    if (cdata[s] != '\n') {
      xcov = atof ((const char *) cdata + s);
      while ((cdata[e] > ' ') && (e < csize)) e += 1;
    }
    /* Skip whitespace */
    s = e;
    while ((cdata[s] <= ' ') && (cdata[s] != '\n') && (s < csize)) s += 1;
    if (s >= csize) break;
    /* Read Y coverage */
    e = s;
    if (cdata[s] != '\n') {
      ycov = atof ((const char *) cdata + s);
      while ((cdata[e] > ' ') && (e < csize)) e += 1;
    }
    
    if ((size >= min_cov) && (size <= MAX_COVERAGE)) {
      if (!pop->min_coverage || (size < pop->min_coverage)) {
        pop->min_coverage = size;
      }
      if (size > pop->max_coverage) {
        pop->max_coverage = size;
      }
      pop->lists[pop->nlists].filename = strdup (b);
      pop->lists[pop->nlists].coverage = size;
      pop->lists[pop->nlists].x_coverage = xcov;
      pop->lists[pop->nlists].y_coverage = ycov;
      if (xcov > 10 * ycov) {
        num_f += 1;
      } else {
        num_m += 1;
      }
      pop->nlists += 1;
    }
    if (pop->nlists >= max_lists) break;
    cpos = e;
  }
  pop->f_fraction = (float) num_f / (num_f + num_m);
  pop->m_fraction = (float) num_m / (num_f + num_m);
  return pop->nlists;
}

