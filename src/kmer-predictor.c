/*
 * Predicts trait by kmer counts
 */

int debug = 0;

#define MAX_LISTS 1024

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "version.h"
#include "word-list-stream.h"

static void linear_regression (double *x, double *y, unsigned int n_samples, double *a, double *b, double *r);
static void print_help (FILE *ofs, int exitvalue);

#define DELTA 20

static unsigned long long max_kmers = 1000000000;

int
main (int argc, const char *argv[]) {
  const char *lists_name = NULL;
  const char *kmers_name = NULL;
  const char *write_coeffs_name = NULL;

  const unsigned char *cdata;
  long long unsigned int csize, cpos;
  char *sample_names[MAX_LISTS];
  char *list_names[MAX_LISTS];
  double ffs[MAX_LISTS];
  unsigned int n_lists = 0;
  double avg_ff;
  double *a, *b, *scale;
  double pred_ffs[MAX_LISTS];
  unsigned int i;

  GT4WordListStream *all_list, *lists[MAX_LISTS];
  GT4WordSArrayImplementation *all_impl, *impls[MAX_LISTS];
  GT4WordSArrayInstance *all_inst, *insts[MAX_LISTS];
  
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--version")) {
      fprintf (stdout, "kmer_predictor version %u.%u.%u (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_MICRO, VERSION_QUALIFIER);
      return 0;
    } else if (!strcmp (argv[i], "-h") || !strcmp (argv[i], "--help") || !strcmp (argv[i], "-?")) {
      print_help (stdout, 0);
    } else if (!strcmp (argv[i], "--kmers")) {
      i += 1;
      if (i >= (unsigned int) argc) {
        print_help (stderr, 1);
      }
      kmers_name = argv[i];
    } else if (!strcmp (argv[i], "--lists")) {
      i += 1;
      if (i >= (unsigned int) argc) {
        print_help (stderr, 1);
      }
      lists_name = argv[i];
    } else if (!strcmp (argv[i], "--write_coefficients")) {
      i += 1;
      if (i >= (unsigned int) argc) {
        print_help (stderr, 1);
      }
      write_coeffs_name = argv[i];
    } else if (!strcmp (argv[i], "--max_kmers")) {
      i += 1;
      if (i >= (unsigned int) argc) {
        print_help (stderr, 1);
      }
      max_kmers = strtoll (argv[i], NULL, 10);
    } else if (!strcmp (argv[i], "-D")) {
      debug += 1;
    } else {
      fprintf(stderr, "Unknown argument: %s!\n", argv[i]);
      print_help (stderr, 1);
    }
  }

  if (!kmers_name || !lists_name) {
    print_help (stderr, 1);
  }

  cdata = gt4_mmap (lists_name, &csize);
  if (!cdata) {
    fprintf (stderr, "Cannot mmap %s\n", lists_name);
    exit (1);
  }
  cpos = 0;
  while (cpos < csize) {
    const unsigned char *tokenz[4];
    unsigned int lengths[4], n_tokenz;
    uint64_t lend = cpos;
    while ((lend < csize) && (cdata[lend] != '\n')) lend += 1;
    n_tokenz = split_line (cdata + cpos, lend - cpos, tokenz, lengths, 4);
    if (n_tokenz == 3) {
      char n[32];
      sample_names[n_lists] = strndup ((const char *) tokenz[0], lengths[0]);
      list_names[n_lists] = strndup ((const char *) tokenz[1], lengths[1]);
      strncpy (n, (const char *) tokenz[2], lengths[2]);
      ffs[n_lists] = atof (n);
      n_lists += 1;
      if (n_lists >= MAX_LISTS) break;
    }
    cpos = lend + 1;
  }
  gt4_munmap (cdata, csize);
  if (debug) fprintf (stderr, "Num lists %u\n", n_lists);
  /* Normalize FF */
  avg_ff = 0;
  for (i = 0; i < n_lists; i++) avg_ff += ffs[i];
  avg_ff /= n_lists;
  for (i = 0; i < n_lists; i++) {
    ffs[i] -= avg_ff;
    if (debug > 1) fprintf (stderr, "%s\t%s\t%.2f\n", sample_names[i], list_names[i], ffs[i]);
  }

  /* Find regression */
  all_list = gt4_word_list_stream_new (kmers_name, VERSION_MAJOR);
  all_impl = (GT4WordSArrayImplementation *) az_object_get_interface (AZ_OBJECT(all_list), GT4_TYPE_WORD_SARRAY, (void **) &all_inst);
  if (debug) fprintf (stdout, "Num words: %llu\n", all_inst->num_words);
  gt4_word_sarray_get_first_word (all_impl, all_inst);
  for (i = 0; i < n_lists; i++) {
    lists[i] = gt4_word_list_stream_new (list_names[i], VERSION_MAJOR);
    if (!lists[i]) {
      fprintf (stderr, "Cannot open list %s\n", list_names[i]);
      exit (1);
    }
    impls[i] = (GT4WordSArrayImplementation *) az_object_get_interface (AZ_OBJECT(lists[i]), GT4_TYPE_WORD_SARRAY, (void **) &insts[i]);
    gt4_word_sarray_get_first_word (impls[i], insts[i]);
  }
  a = (double *) malloc (all_inst->num_words * sizeof (double));
  b = (double *) malloc (all_inst->num_words * sizeof (double));
  scale = (double *) malloc (all_inst->num_words * sizeof (double));
  if (debug) fprintf (stderr, "Calculating correlation coefficients ");
  while (all_inst->idx < all_inst->num_words) {
    double counts[MAX_LISTS];
    double r, t, n_zero = 0;
    double avg_0 = 0, avg_1 = 0;
    double count_0 = 0, count_1 = 0;
    if (debug && !(all_inst->idx % 10000)) fprintf (stderr, ".");
    if (debug > 2) fprintf (stdout, "%llu", all_inst->idx);
    for (i = 0; i < (n_lists - DELTA); i++) {
      while ((insts[i]->idx < insts[i]->num_words) && (insts[i]->word < all_inst->word)) {
        gt4_word_sarray_get_next_word (impls[i], insts[i]);
      }
      if ((insts[i]->idx < insts[i]->num_words) && (insts[i]->word == all_inst->word)) {
        if (debug > 2) fprintf (stdout, "\t%u", insts[i]->count);
        counts[i] = insts[i]->count;
        avg_1 += insts[i]->count * ffs[i];
        count_1 += insts[i]->count;
      } else {
        if (debug > 2) fprintf (stdout, "\t0");
        counts[i] = 0;
        avg_0 += ffs[i];
        count_0 += 1;
        n_zero += 1;
      }
    }
    if (debug > 2) fprintf (stdout, "\n");
    /* linear_regression (counts, ffs, n_lists - DELTA, &a[all_inst->idx], &b[all_inst->idx], &r);
    t = r * sqrtf ((n_lists - 2) / (1 - r * r));*/
    if (count_0) {
      a[all_inst->idx] = avg_0 / count_0;
    } else {
      a[all_inst->idx] = 0;
    }
    if (count_1) {
      b[all_inst->idx] = avg_1 / count_1;
    } else {
      b[all_inst->idx] = 0;
    }
    /* fprintf (stderr, "%llu %g %g %g %g\n", all_inst->idx, avg_0, avg_1, a[all_inst->idx], b[all_inst->idx]); */
    scale[all_inst->idx] = (n_zero * (n_lists - n_zero)) / (n_lists * n_lists);
    /* fprintf (stderr, "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", a[all_inst->idx], b[all_inst->idx], r, t, scale[all_inst->idx]); */
    gt4_word_sarray_get_next_word (all_impl, all_inst);
    if (all_inst->idx >= max_kmers) break;
  }
  if (debug) fprintf (stderr, "\n");
  /* Predict ffs */
  if (debug) fprintf (stderr, "Predicting ");
  for (i = 0; i < n_lists; i++) {
    if (debug) fprintf (stderr, ".");
    pred_ffs[i] = 0;
    gt4_word_sarray_get_first_word (all_impl, all_inst);
    gt4_word_sarray_get_first_word (impls[i], insts[i]);
    while (all_inst->idx < all_inst->num_words) {
      double count = 0;
      while ((insts[i]->idx < insts[i]->num_words) && (insts[i]->word < all_inst->word)) {
        gt4_word_sarray_get_next_word (impls[i], insts[i]);
      }
      if ((insts[i]->idx < insts[i]->num_words) && (insts[i]->word == all_inst->word)) {
        count = insts[i]->count;
      }
      /* pred_ffs[i] += scale[all_inst->idx] * (a[all_inst->idx] + b[all_inst->idx] * count); */
      if (!count) {
        pred_ffs[i] += a[all_inst->idx];
      } else {
        pred_ffs[i] += b[all_inst->idx];
      }
      gt4_word_sarray_get_next_word (all_impl, all_inst);
      if (all_inst->idx >= max_kmers) break;
    }
    /* fprintf (stderr, "%s\t%.3f\t%.3f\n", sample_names[i], ffs[i], pred); */
  }
  if (debug) fprintf (stderr, "\n");
  double pa, pb, pr;
  linear_regression (pred_ffs, ffs, n_lists - DELTA, &pa, &pb, &pr);
  if (write_coeffs_name) {
    FILE *ofs = fopen (write_coeffs_name, "w+");
    fprintf (ofs, "AVG_FF\t%.3g\n", avg_ff);
    fprintf (ofs, "SCALE\t%g\t%g\t%g\n", pa, pb, pr);
    for (i = 0; i < all_inst->num_words; i++) {
      if (i >= max_kmers) break;
      fprintf (ofs, "%g\t%g\n", a[i], b[i]);
    }
    fclose (ofs);
  }
  for (i = 0; i < n_lists; i++) {
    double pred;
    pred = pa + pb * pred_ffs[i];
    fprintf (stderr, "%s\t%.3f\t%.3f\n", sample_names[i], ffs[i] + avg_ff, pred + avg_ff);
  }

  gt4_word_list_stream_delete ((GT4WordListStream *) all_list);
  for (i = 0; i < n_lists; i++) {
    gt4_word_list_stream_delete ((GT4WordListStream *) lists[i]);
  }

  return 0;
}

static void
linear_regression (double *x, double *y, unsigned int n, double *a, double *b, double *r)
{
  double sx = 0, sy = 0, sx2 = 0, sy2 = 0, sxy = 0;
  double d;
  unsigned int i;
  for (i = 0; i < n; i++) {
    sx += x[i];
    sy += y[i];
    sx2 += x[i] * x[i];
    sy2 += y[i] * y[i];
    sxy += x[i] * y[i];
  }
  d = n * sx2 - sx * sx;
  if (d == 0) {
    *a = *b = *r = 0;
    return;
  }
  *a = (sy * sx2 - sx * sxy) / d;
  *b = (n * sxy - sx * sy) / d;
  d = (n * sx2 - sx * sx) * (n * sy2 - sy * sy);
  if (d <= 0) {
    *a = *b = *r = 0;
    return;
  }
}

static void
print_help (FILE *ofs, int exit_value)
{
  fprintf (ofs, "kmer_predictor version %u.%u.%u (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_MICRO, VERSION_QUALIFIER);
  fprintf (ofs, "Usage: kmer_predictor OPTIONS\n");
  fprintf (ofs, "Options:\n");
  fprintf (ofs, "    -v, --version            - print version information and exit\n");
  fprintf (ofs, "    -h, --help               - print this usage screen and exit\n");
  fprintf (ofs, "    -D                       - increase debug level\n");
  exit (exit_value);
}

