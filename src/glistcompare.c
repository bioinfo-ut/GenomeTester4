/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 * 
 * Cpyright (C) 2014 University of Tartu
 *
 * Authors: Maarja Lepamets and Lauris Kaplinski
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
               
#include "common.h"
#include "index-map.h"
#include "utils.h"
#include "sequence.h"
#include "listmaker-queue.h"
#include "version.h"
#include "word-array-sorted.h"
#include "word-list-stream.h"
#include "word-map.h"
#include "word-table.h"

#define USE_SCOUTS use_scouts

enum Rules {
  RULE_DEFAULT,
  RULE_ADD,
  RULE_SUBTRACT,
  RULE_MIN,
  RULE_MAX,
  RULE_FIRST,
  RULE_SECOND,
  RULE_ONE,
  RULE_TWO
};

enum SubsetMethods {
  /* Choose n k-mers from all kmers, collate counts */
  RAND_ALL,
  /* Choose n unique k-mers, copy counts */
  RAND_UNIQUE,
  /* Choose n unique k-mers weighted by counts, copy counts */
  RAND_WEIGHTED_UNIQUE
};

/* Wordmap correctness is not tested */
static int compare_wordmaps (AZObject *list1, AZObject *list2, int find_union, int find_intrsec, int find_diff, int find_ddiff, int subtract, int countonly, const char *out, unsigned int cutoff, int rule);
static int compare_wordmaps_mm (AZObject *list1, AZObject *list2, int find_diff, int find_ddiff, int subtract, int countonly, const char *out, unsigned int cutoff, unsigned int nmm, int rule);
/* No actual writing will be done if ofile is 0 */
/* Upon completion list header is filled regardless of whether actual writing was performed */
static unsigned int union_multi (AZObject *m[], unsigned int nmaps, unsigned int cutoff, int ofile, GT4ListHeader *header);
static unsigned int subset (GT4WordSListImplementation *impl, GT4WordSListInstance *inst, unsigned int subset_method, unsigned long long subset_size, const char *filename);
static unsigned long long fetch_relevant_words (GT4WordTable *table, AZObject *map, AZObject *querymap, unsigned int cutoff, unsigned int nmm, FILE *f, int subtract, int countonly, unsigned long long *totalfreq);
static void print_help (int exitvalue);

#define MAX_FILES 1024

int debug = 0;

unsigned int use_scouts = 1;
unsigned int stream = 0;

int main (int argc, const char *argv[])
{
  int arg_idx, v, i;
  unsigned int nfiles = 0;
  const char *fnames[MAX_FILES];
  AZObject *objs[MAX_FILES];
  char *end;
  int rule = RULE_DEFAULT;
  long seed = -1;
  unsigned int wlen = 0, err = 0;

  /* default values */
  unsigned int cutoff = 1, nmm = 0;
  int find_union = 0, find_intrsec = 0, find_diff = 0, find_ddiff = 0, subtraction = 0, countonly = 0, print_operation = 0;
  int find_subset = 0;
  int subset_method = RAND_ALL;
  unsigned long long subset_size = 0;
  const char *outputname = "out";

  if (argc <= 1) {
    print_help (1);
  }

  for (arg_idx = 1; arg_idx < argc; arg_idx++) {
    if (argv[arg_idx][0] != '-') {
      /* File name */
      if (nfiles >= MAX_FILES) {
        fprintf (stderr, "Too many file arguments (max %d)\n", MAX_FILES);
        print_help (1);
      }
      fnames[nfiles++] = argv[arg_idx];
      continue;
    }
    if (!strcmp(argv[arg_idx], "-v") || !strcmp(argv[arg_idx], "--version")) {
      fprintf (stdout, "glistcompare version %u.%u.%u (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_MICRO, VERSION_QUALIFIER);
      return 0;
    } else if (!strcmp (argv[arg_idx], "-h") || !strcmp (argv[arg_idx], "--help") || !strcmp (argv[arg_idx], "-?")) {
      print_help (0);
    } else if (!strcmp (argv[arg_idx], "-o") || !strcmp (argv[arg_idx], "--outputname")) {
      if (!argv[arg_idx + 1] || argv[arg_idx + 1][0] == '-') {
        fprintf (stderr, "Warning: No output name specified!\n");
        arg_idx += 1;
        continue;
      }
      outputname = argv[arg_idx + 1];
      arg_idx += 1;
    } else if (!strcmp (argv[arg_idx], "-c") || !strcmp (argv[arg_idx], "--cutoff")) {
      if (!argv[arg_idx + 1]) {
        fprintf (stderr, "Warning: No frequency cut-off specified! Using the default value: %d.\n", cutoff);
        continue;
      }
      cutoff = strtol (argv[arg_idx + 1], &end, 10);
      if (*end != 0) {
        fprintf(stderr, "Error: Invalid frequency cut-off: %s! Must be an integer.\n", argv[arg_idx + 1]);
        print_help (1);
      }
      arg_idx += 1;
    } else if (!strcmp (argv[arg_idx], "-mm")|| !strcmp (argv[arg_idx], "--mismatch")) {
      if (!argv[arg_idx + 1]) {
        fprintf (stderr, "Warning: No number of mismatches specified!");
        continue;
      }
      nmm = strtol (argv[arg_idx + 1], &end, 10);
      if (*end != 0) {
        fprintf (stderr, "Error: Invalid number of mismatches: %s! Must be an integer.\n", argv[arg_idx + 1]);
        print_help (1);
      }
      arg_idx += 1;
    } else if (!strcmp (argv[arg_idx], "-u") || !strcmp (argv[arg_idx], "--union")) {
      find_union = 1;
    } else if (!strcmp (argv[arg_idx], "-i") || !strcmp (argv[arg_idx], "--intersection")) {
      find_intrsec = 1;
    } else if (!strcmp (argv[arg_idx], "-d") || !strcmp (argv[arg_idx], "--difference")) {
      find_diff = 1;
    } else if (!strcmp (argv[arg_idx], "-dd") || !strcmp (argv[arg_idx], "--double_difference")) {
      find_ddiff = 1;
    } else if (!strcmp (argv[arg_idx], "-du") || !strcmp (argv[arg_idx], "--diff_union")) {
            find_diff = 1;
      subtraction = 1;
    } else if (!strcmp (argv[arg_idx], "--count_only")) {
      countonly = 1;
    } else if (!strcmp (argv[arg_idx], "-r")|| !strcmp (argv[arg_idx], "--rule")) {
      arg_idx += 1;
      if (arg_idx >= argc) {
        print_help (1);
      }
      if (!strcmp (argv[arg_idx], "default")) {
        rule = RULE_DEFAULT;
                        } else if (!strcmp (argv[arg_idx], "add")) {
        rule = RULE_ADD;
                        } else if (!strcmp (argv[arg_idx], "subtract")) {
        rule = RULE_SUBTRACT;
                        } else if (!strcmp (argv[arg_idx], "min")) {
        rule = RULE_MIN;
                        } else if (!strcmp (argv[arg_idx], "max")) {
        rule = RULE_MAX;
                        } else if (!strcmp (argv[arg_idx], "first")) {
        rule = RULE_FIRST;
                        } else if (!strcmp (argv[arg_idx], "second")) {
        rule = RULE_SECOND;        
                        } else if (!strcmp (argv[arg_idx], "1")) {
        rule = RULE_ONE;        
                        } else if (!strcmp (argv[arg_idx], "2")) {
        rule = RULE_TWO;        
                        }
    } else if (!strcmp (argv[arg_idx], "-ss") || !strcmp (argv[arg_idx], "--subset")) {
      find_subset = 1;
      arg_idx += 1;
      if (arg_idx >= argc) {
        print_help (1);
      }
      if (!strcmp (argv[arg_idx], "rand")) {
        subset_method = RAND_ALL;
      } else if (!strcmp (argv[arg_idx], "rand_unique")) {
        subset_method = RAND_UNIQUE;
      } else if (!strcmp (argv[arg_idx], "rand_weighted_unique")) {
        subset_method = RAND_WEIGHTED_UNIQUE;
      } else {
        print_help (1);
      }
      arg_idx += 1;
      if (arg_idx >= argc) {
        print_help (1);
      }
      subset_size = strtoll (argv[arg_idx], &end, 10);
      if (*end != 0) {
        fprintf (stderr, "Error: Invalid subset size: %s! Must be an integer.\n", argv[arg_idx]);
        print_help (1);
      }
    } else if (!strcmp (argv[arg_idx], "--seed")) {
      arg_idx += 1;
      if (arg_idx >= argc) {
        print_help (1);
      }
      seed = strtoll (argv[arg_idx], &end, 10);
    } else if (!strcmp (argv[arg_idx], "--print_operation")) {
      print_operation = 1;
    } else if (!strcmp (argv[arg_idx], "--disable_scouts")) {
      use_scouts = 0;
    } else if (!strcmp (argv[arg_idx], "--stream")) {
      stream = 1;
    } else if (!strcmp (argv[arg_idx], "-D")) {
      debug += 1;
    } else {
      fprintf(stderr, "Unknown argument: %s!\n", argv[arg_idx]);
      print_help (1);
    }
  }
  if (debug) fprintf (stderr, "Rule: %d\n", rule);

  debug_wordmap = debug;
  
  if (debug) fprintf (stderr, "Num files: %d\n", nfiles);

  if (seed == -1) {
    srand48 ((unsigned int) time(NULL));
  } else {
    srand48 (seed);
  }

  /* Subset is incompatible with --stream */
  if (nmm || find_subset) {
    if (stream) fprintf (stderr, "Warning: Subset and mismatches are incompatible with streaming, using mapping\n");
    stream = 0;
  }

  /* Build list of objects */
  for (i = 0; i < nfiles; i++) {
    GT4WordSListImplementation *impl;
    GT4WordSListInstance *inst;
    FILE *ifs;
    uint32_t code;

    ifs = fopen (fnames[i], "r");
    if (!ifs) {
      fprintf (stderr, "Error: Cannot open %s\n", fnames[i]);
      err = 1;
    }
    fread (&code, 4, 1, ifs);
    fclose (ifs);
    if (code == GT4_LIST_CODE) {
      if (stream) {
        objs[i] = (AZObject *) gt4_word_list_stream_new (fnames[0], VERSION_MAJOR);
      } else {
        objs[i] = (AZObject *) gt4_word_map_new (fnames[i], VERSION_MAJOR, use_scouts, 0);
      }
    } else if (code == GT4_INDEX_CODE) {
      objs[i] = (AZObject *) gt4_index_map_new (fnames[i], VERSION_MAJOR, 0);
    } else {
      fprintf (stderr, "Error: File %s has unknown format\n", fnames[i]);
      err = 1;
    }
    impl = (GT4WordSListImplementation *) az_object_get_interface (objs[i], GT4_TYPE_WORD_SLIST, (void **) &inst);
    if (!impl) {
      fprintf (stderr, "Error: File %s is invalid or corrupted\n", fnames[i]);
      err = 1;
    }
    if (!wlen) {
      wlen = inst->word_length;
    } else if (inst->word_length != wlen) {
      fprintf (stderr, "Error: File %s has different word length (%u != %u)\n", fnames[i], inst->word_length, wlen);
      err = 1;
    }
  }
  if (err) {
    fprintf (stderr, "Stopping...\n");
    exit (1);
  }
  
  /* Subset */
  if (find_subset) {
    GT4WordSListImplementation *impl;
    GT4WordSListInstance *inst;
    char out_name[2048], tmp_name[2048];
    if (nfiles != 1) {
      fprintf (stderr, "Error: Subsetting multiple files is not supported\n");
      exit (1);
    }
    impl = (GT4WordSListImplementation *) az_object_get_interface (objs[0], GT4_TYPE_WORD_SLIST, (void **) &inst);
    if (((subset_method == RAND_UNIQUE) || (subset_method == RAND_WEIGHTED_UNIQUE)) && (subset_size > inst->num_words)) {
      fprintf (stderr, "Error: Unique subset size (%llu) is bigger than number of unique kmers (%llu)\n", subset_size, inst->num_words);
      exit (1);
    }
    snprintf (tmp_name, 2048, "%s_subset_%u.list.tmp", outputname, inst->word_length);
    tmp_name[2047] = 0;
    subset (impl, inst, subset_method, subset_size, tmp_name);
    snprintf (out_name, 2048, "%s_subset_%u.list", outputname, inst->word_length);
    out_name[2047] = 0;
    if (rename (tmp_name, out_name)) {
      fprintf (stderr, "Error: cannot rename %s to %s\n", tmp_name, out_name);
    }
    return 0;
  }

  if (nfiles < 2) {
    fprintf (stderr, "Error: At least 2 list/index files are needed\n");
    exit (1);
  }

  if (nfiles > 2) {
    if (!find_union || find_intrsec || find_diff || find_ddiff) {
      fprintf(stderr, "Error: Algorithm incompatible with multiple files!\n");
      print_help (1);
    }
    if (nmm) {
      fprintf(stderr, "Error: Multiple files are not compatible with mismatches!\n");
      print_help (1);
    }
    if (rule != RULE_DEFAULT) {
      fprintf(stderr, "Error: Explicit rule incompatible with multiple files!\n");
      print_help (1);
    }
  }

  /* both differences */
  if (find_ddiff) find_diff = 1;

  /* checking parameter values */
  if (!find_diff && nmm) fprintf(stderr, "Warning: Number of mismatches are not used!\n");
  if (!find_diff && subtraction) fprintf(stderr, "Warning: Subtraction is not used!\n");
  if (strlen (outputname) > 200) {
    fprintf (stderr, "Error: Output name exceeds the 200 character limit.\n");
    exit (1);
  }
  
  if (!find_intrsec && (rule == RULE_SUBTRACT || rule == RULE_MIN || rule == RULE_FIRST || rule == RULE_SECOND)) {
    fprintf (stderr, "Error: Rules min, subtract, fist and second can only be used with finding the intersection.\n");
    exit (1);
  }

  if (print_operation) {
    unsigned int i;
    fprintf (stdout, "Operation\t%s%s%s%s\nFiles\t%u\n", (find_union) ? "U" : "", (find_intrsec) ? "I" : "", (find_diff) ? "D" : "", (find_ddiff) ? "X" : "", nfiles);
    for (i = 0; i < nfiles; i++) {
      fprintf (stdout, "%u\t%s\n", i, fnames[i]);
    }
  }

  if (nmm) {
    v = compare_wordmaps_mm (objs[0], objs[1], find_diff, find_ddiff, subtraction, countonly, outputname, cutoff, nmm, rule);
  } else if (nfiles == 2) {
    v = compare_wordmaps (objs[0], objs[1], find_union, find_intrsec, find_diff, find_ddiff, subtraction, countonly, outputname, cutoff, rule);
  } else {
    GT4ListHeader header;
    char out_name[2048], tmp_name[2048];
    int ofile = 0;
    if (!countonly) {
      snprintf (tmp_name, 2048, "%s_%d_union.list.tmp", outputname, wlen);
      tmp_name[2047] = 0;
      ofile = creat (tmp_name, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
      if (ofile < 0) {
        fprintf (stderr, "Error: Cannot create output file %s\n", tmp_name);
        exit (1);
      }
    }
    v = union_multi (objs, nfiles, cutoff, ofile, &header);
    if (ofile > 0) {
      close (ofile);
      snprintf (out_name, 2048, "%s_%d_union.list", outputname, wlen);
      out_name[2047] = 0;
      if (rename (tmp_name, out_name)) {
        fprintf (stderr, "Error: Cannot rename %s to %s\n", tmp_name, out_name);
        exit (1);
      }
    }
    if (debug) fprintf (stdout, "NUnique\t%llu\nNTotal\t%llu\n", header.n_words, header.total_count);
  }
  if (v) return print_error_message (v);
  for (i = 0; i < nfiles; i++) {
    az_object_shutdown (objs[i]);
  }

  delete_scouts ();
  
  return 0;
}

/* Calculates resulting freq from rule */

static unsigned int
calculate_freq (unsigned int freq1, unsigned int freq2, int rule)
{
  switch (rule) {
  case RULE_ADD:
    return freq1 + freq2;
  case RULE_SUBTRACT:
    return (freq1 > freq2) ? freq1 - freq2 : 0;
  case RULE_MIN:
    return (freq1 < freq2) ? freq1 : freq2;
  case RULE_MAX:
    return (freq1 > freq2) ? freq1 : freq2;
  case RULE_FIRST:
    return freq1;
  case RULE_SECOND:
    return freq2;
  case RULE_ONE:
    return 1;
  case RULE_TWO:
    return 2;
  default:
    break;
  }
  return 0;
}

/* Decides whether given word has to be written to union and calculates frequency */

static unsigned int
include_in_union (unsigned int freq1, unsigned int freq2, unsigned int *freq, int rule, unsigned int cutoff)
{
  if ((freq1 < cutoff) && (freq2 < cutoff)) return 0;
  if (rule == RULE_DEFAULT) rule = RULE_ADD;
  *freq = calculate_freq (freq1, freq2, rule);
  return *freq != 0;
}

static unsigned int
include_in_intersection (unsigned int freq1, unsigned int freq2, unsigned int *freq, int rule, unsigned int cutoff)
{
  if ((freq1 < cutoff) || (freq2 < cutoff)) return 0;
  if (rule == RULE_DEFAULT) rule = RULE_MIN;
  *freq = calculate_freq (freq1, freq2, rule);
  return *freq != 0;
}

static unsigned int
include_in_complement (unsigned int freq1, unsigned int freq2, unsigned int *freq, int rule, unsigned int cutoff, unsigned int subtract)
{
  if (subtract) {
     if ((freq1 != freq2) || (freq1 < cutoff)) return 0;
     *freq = freq1;
     return 1;
  }
  if ((freq1 < cutoff) || (freq2 >= cutoff)) return 0;
  if (rule == RULE_DEFAULT) rule = RULE_SUBTRACT;
  *freq = calculate_freq (freq1, freq2, rule);
  return *freq != 0;
}

static void
write_word_to_file (unsigned long long word, unsigned freq, FILE *f)
{
  fwrite (&word, sizeof (unsigned long long), 1, f);
  fwrite (&freq, sizeof (unsigned int), 1, f);
}
                        
#define TMP_BUF_SIZE (256 * 12)

static unsigned int
union_multi (AZObject *m[], unsigned int nmaps, unsigned int cutoff, int ofile, GT4ListHeader *header)
{
  GT4WordSListImplementation *impls[MAX_FILES];
  GT4WordSListInstance *insts[MAX_FILES];
  unsigned int n_sources;
  unsigned int j;
  unsigned long long word;
  unsigned char b[TMP_BUF_SIZE];
  unsigned int bp = 0;
  unsigned long long total = 0;
  double t_s, t_e;
  
  n_sources = 0;
  for (j = 0; j < nmaps; j++) {
    impls[n_sources] = (GT4WordSListImplementation *) az_object_get_interface (AZ_OBJECT(m[j]), GT4_TYPE_WORD_SLIST, (void **) &insts[n_sources]);
    if (insts[n_sources]->num_words) {
      gt4_word_slist_get_first_word (impls[n_sources], insts[n_sources]);
      total += insts[n_sources]->num_words;
      n_sources += 1;
    }
  }

  gt4_list_header_init (header, insts[0]->word_length);

  t_s = get_time ();

  if (ofile) write (ofile, header, sizeof (GT4ListHeader));

  word = 0xffffffffffffffff;
  for (j = 0; j < n_sources; j++) {
    if (insts[j]->word < word) word = insts[j]->word;
  }
    
  while (n_sources) {
    unsigned long long next = 0xffffffffffffffff;
    unsigned int freq = 0;
    j = 0;
    while (j < n_sources) {
      if (insts[j]->word == word) {
        freq += insts[j]->count;
        if (!gt4_word_slist_get_next_word (impls[j], insts[j])) {
          n_sources -= 1;
          if (n_sources > 0) {
            impls[j] = impls[n_sources];
            insts[j] = insts[n_sources];
            continue;
          } else {
            break;
          }
        }
      }
      if (insts[j]->word < next) next = insts[j]->word;
      j += 1;
    }
    
    /* Now we have word and freq */
    if (freq >= cutoff) {
      if (ofile) {
        memcpy (&b[bp], &word, 8);
        memcpy (&b[bp + 8], &freq, 4);
        bp += 12;
        if (bp >= TMP_BUF_SIZE) {
          write (ofile, b, bp);
          bp = 0;
        }
      }
      header->n_words += 1;
      header->total_count += freq;
      if (debug && !(header->n_words % 100000000)) {
        fprintf (stderr, "Words written: %uM\n", (unsigned int) (header->n_words / 1000000));
      }
    }
    word = next;
  }
  if (ofile) {
    if (bp) write (ofile, b, bp);
    pwrite (ofile, header, sizeof (GT4ListHeader), 0);
  }
  t_e = get_time ();

  if (debug > 0) {
    fprintf (stderr, "Combined %u maps: input %llu (%.3f Mwords/s) output %llu (%.3f Mwords/s)\n", nmaps, total, total / (1000000 * (t_e - t_s)), header->n_words, header->n_words / (1000000 * (t_e - t_s)));
  }
  
  return 0;
}

static unsigned int
subset (GT4WordSListImplementation *impl, GT4WordSListInstance *inst, unsigned int subset_method, unsigned long long subset_size, const char *filename)
{
  GT4ListHeader h_out;
  FILE *ofs;
  uint64_t in = 0;
  uint64_t out = subset_size;

  gt4_list_header_init (&h_out, inst->word_length);

  ofs = fopen (filename, "w");
  fwrite (&h_out, sizeof (GT4ListHeader), 1, ofs);

  gt4_word_slist_get_first_word (impl, inst);
  if (subset_method == RAND_ALL) {
    in = inst->sum_counts;
    while (out > 0) {
      unsigned int count = 0, i;
      for (i = 0; (i < inst->count) && (out > 0); i++) {
        double val = drand48 ();
        if (val <= ((double) out / in)) {
          count += 1;
          out -= 1;
        }
        in -= 1;
      }
      if (count > 0) {
        fwrite (&inst->word, sizeof (unsigned long long), 1, ofs);
        fwrite (&count, sizeof (unsigned int), 1, ofs);
        h_out.n_words += 1;
        h_out.total_count += count;
      }
      gt4_word_slist_get_next_word (impl, inst);
    }
  } else if (subset_method == RAND_UNIQUE) {
    in = inst->num_words;
    while (out > 0) {
      double val = drand48 ();
      if (val <= ((double) out / in)) {
        fwrite (&inst->word, sizeof (unsigned long long), 1, ofs);
        fwrite (&inst->count, sizeof (unsigned int), 1, ofs);
        h_out.n_words += 1;
        h_out.total_count += inst->count;
        out -= 1;
      }
      in -= 1;
      gt4_word_slist_get_next_word (impl, inst);
    }
  } else if (subset_method == RAND_WEIGHTED_UNIQUE) {
    in = inst->sum_counts;
    while (out > 0) {
      double val = drand48 ();
      if (val <= ((double) inst->count * out / in)) {
        fwrite (&inst->word, sizeof (unsigned long long), 1, ofs);
        fwrite (&inst->count, sizeof (unsigned int), 1, ofs);
        h_out.n_words += 1;
        h_out.total_count += inst->count;
        out -= 1;
      }
      in -= inst->count;
      gt4_word_slist_get_next_word (impl, inst);
    }
  }
        
  fseek (ofs, 0, SEEK_SET);
  fwrite (&h_out, sizeof (GT4ListHeader), 1, ofs);
  fclose (ofs);
  return 0;
}

static int
compare_wordmaps (AZObject *list1, AZObject *list2, int find_union, int find_intrsec, int find_diff, int find_ddiff, int subtract, int countonly, const char *out, unsigned int cutoff, int rule)
{
  GT4WordSListImplementation *impl1, *impl2;
  GT4WordSListInstance *inst1, *inst2;
  FILE *outf[4] = { 0 };
  /* the length is limited in main(..) method */
  char fname[300][4], name[256];
  GT4ListHeader h_out;

  unsigned long long word1, word2;
  unsigned int freq1, freq2;
  unsigned long long c_union = 0L, c_inters = 0L, c_diff1 = 0L, c_diff2 = 0L;
  unsigned long long freqsum_union = 0L, freqsum_inters = 0L, freqsum_diff1 = 0L, freqsum_diff2 = 0L;

  impl1 = (GT4WordSListImplementation *) az_object_get_interface (list1, GT4_TYPE_WORD_SLIST, (void **) &inst1);
  impl2 = (GT4WordSListImplementation *) az_object_get_interface (list2, GT4_TYPE_WORD_SLIST, (void **) &inst2);

  if (debug) {
     fprintf (stderr, "compare_wordmaps: List 1: %llu entries\n", inst1->num_words);
     fprintf (stderr, "compare_wordmaps; List 2: %llu entries\n", inst2->num_words);
  }

  /* creating output files */
  gt4_list_header_init (&h_out, inst1->word_length);
  if (find_union && !countonly) {
    sprintf (fname[0], "%s_%d_union.list.tmp", out, inst1->word_length);
    outf[0] = fopen (fname[0], "w");
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[0]);
  }
  if (find_intrsec && !countonly) {
    sprintf (fname[1], "%s_%d_intrsec.list.tmp", out, inst1->word_length);
    outf[1] = fopen (fname[1], "w");
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[1]);
  }
  if (find_diff && !countonly) {
    sprintf (fname[2], "%s_%d_0_diff1.list.tmp", out, inst1->word_length);
    outf[2] = fopen (fname[2], "w");
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[2]);
  }
  if (find_ddiff && !countonly) {
    sprintf (fname[3], "%s_%d_0_diff2.list.tmp", out, inst1->word_length);
    outf[3] = fopen (fname[3], "w");
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[3]);
  }

  gt4_word_slist_get_first_word (impl1, inst1);
  word1 = inst1->word;
  freq1 = inst1->count;
  gt4_word_slist_get_first_word (impl2, inst2);
  word2 = inst2->word;
  freq2 = inst2->count;

  while ((inst1->idx < inst1->num_words) || (inst2->idx < inst2->num_words)) {
    unsigned int freq = 0;
    if ((inst1->idx < inst1->num_words) && (inst2->idx < inst2->num_words) && (word1 == word2)) {
      /* Two words are equal */
      if (find_union && include_in_union (freq1, freq2, &freq, rule, cutoff)) {
        if (!countonly) write_word_to_file (word1, freq, outf[0]);
        c_union += 1;
        freqsum_union += freq;
      }
      if (find_intrsec && include_in_intersection (freq1, freq2, &freq, rule, cutoff)) {
        if (!countonly) write_word_to_file (word1, freq, outf[1]);
        c_inters += 1;
        freqsum_inters += freq;
      }
      if (find_diff && include_in_complement (freq1, freq2, &freq, rule, cutoff, subtract)) {
        if (!countonly) write_word_to_file (word1, freq, outf[2]);
        freqsum_diff1 += freq;
        c_diff1 += 1;
      }
      if (find_ddiff && include_in_complement (freq2, freq1, &freq, rule, cutoff, 0)) {
        if (!countonly) write_word_to_file (word2, freq, outf[3]);
        freqsum_diff2 += freq;
        c_diff2 += 1;
      }
      gt4_word_slist_get_next_word (impl1, inst1);
      word1 = inst1->word;
      freq1 = inst1->count;
      gt4_word_slist_get_next_word (impl2, inst2);
      word2 = inst2->word;
      freq2 = inst2->count;
    } else if ((inst1->idx < inst1->num_words) && ((inst2->idx >= inst2->num_words) || (word1 < word2))) {
      /* Second is EOF or first is smaller than second */
      if (find_union && include_in_union (freq1, 0, &freq, rule, cutoff)) {
        if (!countonly) write_word_to_file (word1, freq, outf[0]);
        c_union += 1;
        freqsum_union += freq;
      }
      if (find_diff && include_in_complement (freq1, 0, &freq, rule, cutoff, subtract)) {
        if (!countonly) write_word_to_file (word1, freq, outf[2]);
        freqsum_diff1 += freq;
        c_diff1 += 1;
        
      }
      gt4_word_slist_get_next_word (impl1, inst1);
      word1 = inst1->word;
      freq1 = inst1->count;
    } else if ((inst2->idx < inst2->num_words) && ((inst1->idx >= inst1->num_words) || (word2 < word1))) {
      /* First is EOF or second is smaller than first */
      if (find_union && include_in_union (0, freq2, &freq, rule, cutoff)) {
        if (!countonly) write_word_to_file (word2, freq, outf[0]);
        c_union += 1;
        freqsum_union += freq;
      }
      if (find_ddiff && include_in_complement (freq2, 0, &freq, rule, cutoff, 0)) {
        if (!countonly) write_word_to_file (word2, freq, outf[3]);
        freqsum_diff2 += freq;
        c_diff2 += 1;
      }
      gt4_word_slist_get_next_word (impl2, inst2);
      word2 = inst2->word;
      freq2 = inst2->count;
    }
  }

  /* add headers and close files */
  if (find_union && !countonly) {
    h_out.n_words = c_union;
    h_out.total_count = freqsum_union;
    fseek (outf[0], 0, SEEK_SET);
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[0]);
    fclose (outf[0]);
    sprintf (name, "%s_%d_union.list", out, inst1->word_length);
    rename (fname[0], name);
  } else if (find_union) {
    fprintf (stdout, "NUnique\t%llu\nNTotal\t%llu\n", c_union, freqsum_union);
  }
  if (find_intrsec && !countonly) {
    h_out.n_words = c_inters;
    h_out.total_count = freqsum_inters;
    fseek (outf[1], 0, SEEK_SET);
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[1]);
    fclose (outf[1]);
    sprintf (name, "%s_%d_intrsec.list", out, inst1->word_length);
    rename (fname[1], name);
  } else if (find_intrsec) {
    fprintf (stdout, "NUnique\t%llu\nNTotal\t%llu\n", c_inters, freqsum_inters);
  }
  if (find_diff && !countonly) {
    h_out.n_words = c_diff1;
    h_out.total_count = freqsum_diff1;
    fseek (outf[2], 0, SEEK_SET);
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[2]);
    fclose (outf[2]);
    sprintf (name, "%s_%d_0_diff1.list", out, inst1->word_length);
    rename (fname[2], name);
  } else if (find_diff) {
    fprintf (stdout, "NUnique\t%llu\nNTotal\t%llu\n", c_diff1, freqsum_diff1);
  }
  if (find_ddiff && !countonly) {
    h_out.n_words = c_diff2;
    h_out.total_count = freqsum_diff2;
    fseek (outf[3], 0, SEEK_SET);
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[3]);
    fclose (outf[3]);
    sprintf (name, "%s_%d_0_diff2.list", out, inst1->word_length);
    rename (fname[3], name);
  } else if (find_ddiff) {
    fprintf (stdout, "NUnique\t%llu\nNTotal\t%llu\n", c_diff2, freqsum_diff2);
  }
  return 0;
}

static int
compare_wordmaps_mm (AZObject *list1, AZObject *list2, int find_diff, int find_ddiff, int subtract, int countonly, const char *out, unsigned int cutoff, unsigned int nmm, int rule)
{
  GT4WordSListImplementation *impl1, *impl2;
  GT4WordSListInstance *inst1, *inst2;
  GT4WordTable tbl_diff, tbl_ddiff;
  FILE *outf[4] = { 0 };
  GT4ListHeader h_out;
  unsigned long long word1, word2;
  unsigned int freq1, freq2;
  unsigned long long c_diff1 = 0L, c_diff2 = 0L;
  unsigned long long freqsum_diff1 = 0L, freqsum_diff2 = 0L;
  /* the length is limited in main(..) method */
  char fname[2][300], name[256];

  impl1 = (GT4WordSListImplementation *) az_object_get_interface (list1, GT4_TYPE_WORD_SLIST, (void **) &inst1);
  impl2 = (GT4WordSListImplementation *) az_object_get_interface (list2, GT4_TYPE_WORD_SLIST, (void **) &inst2);

  if (debug) {
     fprintf (stderr, "compare_wordmaps: List 1: %llu entries\n", inst1->num_words);
     fprintf (stderr, "compare_wordmaps; List 2: %llu entries\n", inst2->num_words);
  }

  gt4_list_header_init (&h_out, inst1->word_length);
  if (find_diff) {
    gt4_word_table_setup (&tbl_diff, inst1->word_length, 16384, 4);
    if (!countonly) {
      sprintf (fname[0], "%s_%d_%d_diff1.list.tmp", out, inst1->word_length, nmm);
      outf[2] = fopen (fname[0], "w");
      fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[2]);
    }
  }
  if (find_ddiff) {
    gt4_word_table_setup (&tbl_ddiff, inst1->word_length, 16384, 4);
    if (!countonly) {
      sprintf (fname[1], "%s_%d_%d_diff2.list.tmp", out, inst1->word_length, nmm);
      outf[3] = fopen (fname[1], "w");
      fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[3]);
    }
  }

  gt4_word_slist_get_first_word (impl1, inst1);
  word1 = inst1->word;
  freq1 = inst1->count;
  gt4_word_slist_get_first_word (impl2, inst2);
  word2 = inst2->word;
  freq2 = inst2->count;

  if (debug) {
    fprintf (stderr, "Table 1: %llu entries\n", inst1->num_words);
    fprintf (stderr, "Table 2: %llu entries\n", inst2->num_words);
  }
  while ((inst1->idx < inst1->num_words) || (inst2->idx < inst2->num_words)) {
    unsigned int first_ge_cutoff = (freq1 >= cutoff);
    unsigned int second_ge_cutoff = (freq2 >= cutoff);
    if ((inst1->idx < inst1->num_words) && (inst2->idx < inst2->num_words) && (word1 == word2)) {
      /* Two words are equal */
      if (find_diff) {
        if (subtract) {
          if (freq1 <= freq2) {
            freq2 -= freq1;
          }
        }
        if (first_ge_cutoff && !second_ge_cutoff) {
          unsigned int freq = freq1 - freq2;
          gt4_word_table_add_word (&tbl_diff, word1, &freq);
          c_diff1 += 1;
        }
      }
      if (find_ddiff && second_ge_cutoff && !first_ge_cutoff) {
        unsigned int freq = freq2 - freq1;
        gt4_word_table_add_word (&tbl_ddiff, word2, &freq);
        c_diff2 += 1;
      }
      gt4_word_slist_get_next_word (impl1, inst1);
      word1 = inst1->word;
      freq1 = inst1->count;
      gt4_word_slist_get_next_word (impl2, inst2);
      word2 = inst2->word;
      freq2 = inst2->count;
    } else if ((inst1->idx < inst1->num_words) && ((inst2->idx >= inst2->num_words) || (word1 < word2))) {
      /* Second is EOF or first is smaller than second */
      if (find_diff && first_ge_cutoff && !subtract) {
        gt4_word_table_add_word (&tbl_diff, word1, &freq1);
        c_diff1 += 1;
      }
      gt4_word_slist_get_next_word (impl1, inst1);
      word1 = inst1->word;
      freq1 = inst1->count;
    } else if ((inst2->idx < inst2->num_words) && ((inst1->idx >= inst1->num_words) || (word2 < word1))) {
      /* First is EOF or second is smaller than first */
      if (find_ddiff && second_ge_cutoff) {
        gt4_word_table_add_word (&tbl_ddiff, word2, &freq2);
        c_diff2 += 1;
      }
      gt4_word_slist_get_next_word (impl2, inst2);
      word2 = inst2->word;
      freq2 = inst2->count;
    }
  }

  /* finding the mismatches */
  if (find_diff) {
    if (debug > 0) {
      fprintf (stderr, "Finding diff with mismatches (%llu entries)\n", tbl_diff.n_words);
    }
    c_diff1 = fetch_relevant_words (&tbl_diff, list2, list1, cutoff, nmm, outf[2], subtract, countonly, &freqsum_diff1);
  }
  if (find_ddiff) {
    c_diff2 = fetch_relevant_words (&tbl_ddiff, list1, NULL, cutoff, nmm, outf[3], subtract, countonly, &freqsum_diff2);
  }

  /* Add headers and close files */
  if (find_diff && !countonly) {
    h_out.n_words = c_diff1;
    h_out.total_count = freqsum_diff1;
    fseek (outf[2], 0, SEEK_SET);
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[2]);
    fclose (outf[2]);
    sprintf (name, "%s_%d_%d_diff1.list", out, inst1->word_length, nmm);
    rename (fname[0], name);
  } else if (find_diff) {
    fprintf (stdout, "NUnique\t%llu\nNTotal\t%llu\n", c_diff1, freqsum_diff1);
  }
  if (find_ddiff && !countonly) {
    h_out.n_words = c_diff2;
    h_out.total_count = freqsum_diff2;
    fseek (outf[3], 0, SEEK_SET);
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[3]);
    fclose (outf[3]);
    sprintf (name, "%s_%d_%d_diff2.list", out, inst1->word_length, nmm);
    rename (fname[1], name);
  } else if (find_ddiff) {
    fprintf (stdout, "NUnique\t%llu\nNTotal\t%llu\n", c_diff2, freqsum_diff2);
  }
  return 0;
}

static unsigned int 
search_query (GT4WordDictImplementation *impl_m, GT4WordDictInstance *inst_m, unsigned long long query, unsigned int n_mm, unsigned int equalmmonly,
  unsigned int dosubtraction, GT4WordDictImplementation *impl_q, GT4WordDictInstance *inst_q)
{
  static GT4WordTable mm_table = {0};
  unsigned long long i;
  unsigned int count = 0L, currentcount = 0L, querycount = 0L;

  if (!mm_table.data_size) {
    gt4_word_table_setup (&mm_table, inst_m->word_length, 16384, 0);
  }

  gt4_word_table_generate_mismatches (&mm_table, query, NULL, n_mm, 0, 0, 0, equalmmonly);
  if (debug_wordmap > 1) {
    fprintf (stderr, "MM Table size %llu\n", mm_table.n_words);
  }

  for (i = 0; i < mm_table.n_words; i++) {
    if (dosubtraction) {
      querycount = gt4_word_dict_lookup (impl_q, inst_q, mm_table.words[i], 1);
      currentcount = gt4_word_dict_lookup (impl_m, inst_m, mm_table.words[i], 1);
      if (currentcount > querycount) {
        if (debug_wordmap > 1) {
          fprintf (stderr, "%llu %llu %llu querycount %u currentcount %u\n", query, i, mm_table.words[i], querycount, currentcount);
        }
        mm_table.n_words = 0;
        return ~0L;
      }
      count += (currentcount - querycount);
    } else {
      currentcount = gt4_word_dict_lookup (impl_m, inst_m, mm_table.words[i], 1);
      count += currentcount;
    }
  }
  gt4_word_table_clear (&mm_table);
  return count;
}

static unsigned long long
fetch_relevant_words (GT4WordTable *table, AZObject *map, AZObject *querymap, unsigned int cutoff, unsigned int nmm, FILE *f, int subtract, int countonly, unsigned long long *totalfreq)
{
  GT4WordDictImplementation *impl_q = NULL, *impl_m;
  GT4WordDictInstance *inst_q = NULL, *inst_m;
  unsigned long long ri, wi, word, sumfreq = 0L, count = 0L;
  unsigned int freq, cnmm;
  unsigned int *freqs = (unsigned int *) table->data;

  if (table->n_words == 0) return 0;
  impl_m = (GT4WordDictImplementation *) az_object_get_interface (map, GT4_TYPE_WORD_DICT, (void **) &inst_m);
  if (querymap) impl_q = (GT4WordDictImplementation *) az_object_get_interface (querymap, GT4_TYPE_WORD_DICT, (void **) &inst_q);

  for (cnmm = 1; cnmm <= nmm; cnmm++) {
    wi = 0;
    for (ri = 0; ri < table->n_words; ri++) {
      if (debug > 2) {
        fprintf (stderr, "cnmm %u ri %llu wi %llu\n", cnmm, ri, wi);
      }
      word = table->words[ri];
      freq = freqs[ri];      
      sumfreq = search_query (impl_m, inst_m, word, cnmm, 1, subtract, impl_q, inst_q);
      if (cnmm == nmm && sumfreq < cutoff) {
        if (!countonly) write_word_to_file (word, freq, f);
        count += 1;
        *totalfreq += freq;
        
      } else if (sumfreq < cutoff) {
        table->words[wi] = word;
        freqs[wi] = freq;
        wi += 1;
      }
    }
    table->n_words = wi;
  }
  return count;
}

static void
print_help (int exit_value)
{
  fprintf (stdout, "glistcompare version %u.%u.%u (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_MICRO, VERSION_QUALIFIER);
  fprintf (stdout, "Usage: glistcompare INPUTLIST1 [INPUTLIST2...] METHOD [OPTIONS]\n");
  fprintf (stdout, "Options:\n");
  fprintf (stdout, "    -v, --version            - print version information and exit\n");
  fprintf (stdout, "    -h, --help               - print this usage screen and exit\n");
  fprintf (stdout, "    -u, --union              - union of input lists\n");
  fprintf (stdout, "    -i, --intersection       - intersection of input lists\n");
  fprintf (stdout, "    -d, --difference         - difference of input lists\n");
  fprintf (stdout, "    -dd, --double_difference - double difference of input lists\n");
  fprintf (stdout, "    -du, --diff_union        - subtract first list from the second and finds difference\n");
  fprintf (stdout, "    -mm, --mismatch   NUMBER - specify number of mismatches (default 0, can be used with -diff and -ddiff)\n");
  fprintf (stdout, "    -c, --cutoff NUMBER      - specify frequency cut-off (default 1)\n");
  fprintf (stdout, "    -o, --outputname STRING  - specify output name (default \"out\")\n");
  fprintf (stdout, "    -r, --rule STRING        - specify rule how final frequencies are calculated (default, add, subtract, min, max, first, second, 1, 2)\n");
  fprintf (stdout, "                               NOTE: rules min, subtract, first and second can only be used with finding the intersection.\n");
  fprintf (stdout, "    -ss, --subset METHOD SIZE - make subset with given method (rand, rand_unique, rand_weighted_unique)\n");
  fprintf (stdout, "    --seed INTEGER           - Set seed of random number generator (default uses start time)\n");
  fprintf (stdout, "    --count_only             - output count of k-mers instead of k-mers themself\n");
  fprintf (stdout, "    --disable_scouts         - disable list read-ahead in background thread\n");
  fprintf (stdout, "    --stream                 - read input as stream (do not memory map files)\n");
  fprintf (stdout, "    -D                       - increase debug level\n");
  exit (exit_value);
}
