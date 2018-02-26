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
               
#include "wordtable.h"
#include "word-list-stream.h"
#include "word-map.h"
#include "common.h"
#include "utils.h"
#include "sequence.h"
#include "listmaker-queue.h"
#include "version.h"

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
  RAND_ALL,
  RAND_UNIQUE
};

static int compare_word_map_headers (GT4ListHeader *h1, GT4ListHeader *h2);
/* Wordmap correctness is not tested */
static int compare_wordmaps (AZObject *list1, AZObject *list2, int find_union, int find_intrsec, int find_diff, int find_ddiff, int subtract, int countonly, const char *out, unsigned int cutoff, int rule);
/* No actual writing will be done if ofile is 0 */
/* Upon completion list header is filled regardless of whether actual writing was performed */
static unsigned int union_multi (AZObject *m[], unsigned int nmaps, unsigned int cutoff, int ofile, GT4ListHeader *header);
static unsigned int subset (GT4WordMap *map, unsigned int subset_method, unsigned long long subset_size, const char *filename);
static int compare_wordmaps_mm (GT4WordMap *map1, GT4WordMap *map2, int find_diff, int find_ddiff, int subtract, int countonly, const char *out, unsigned int cutoff, unsigned int nmm, int rule);
static unsigned long long fetch_relevant_words (wordtable *table, GT4WordMap *map, GT4WordMap *querymap, unsigned int cutoff, unsigned int nmm, FILE *f, int subtract, int countonly, unsigned long long *totalfreq);
static void print_help (int exitvalue);

#define MAX_FILES 1024

int debug = 0;

unsigned int use_scouts = 1;
unsigned int stream = 0;

int main (int argc, const char *argv[])
{
  int arg_idx, v;
  const char *fnames[MAX_FILES];
  unsigned int nfiles = 0;
  char *end;
  int rule = RULE_DEFAULT;

  /* default values */
  unsigned int cutoff = 1, nmm = 0;
  int find_union = 0, find_intrsec = 0, find_diff = 0, find_ddiff = 0, subtraction = 0, countonly = 0, print_operation = 0;
  int find_subset = 0;
  int subset_method = RAND_ALL;
  unsigned long long subset_size = 0;
  const char *outputname = "out";

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
      fprintf (stdout, "glistcompare version %d.%d (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_QUALIFIER);
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
      if (arg_idx >= argc) {
        fprintf (stderr, "Warning: No rule specified!");
        exit (1);
      }
      arg_idx += 1;
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
    } else if (!strcmp (argv[arg_idx], "-ss") || !strcmp (argv[arg_idx], "-subset")) {
      find_subset = 1;
      arg_idx += 1;
      if (arg_idx >= argc) {
        print_help (1);
      }
      if (!strcmp (argv[arg_idx], "rand")) {
        subset_method = RAND_ALL;
      } else if (!strcmp (argv[arg_idx], "rand_unique")) {
        subset_method = RAND_UNIQUE;
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
  
  /* check list files */
  if (debug) fprintf (stderr, "Num files: %d\n", nfiles);
  if (nfiles < 2) {
    if (find_subset) {
      GT4WordMap *map;
      char c[2048];
      map = gt4_word_map_new (fnames[0], VERSION_MAJOR, USE_SCOUTS);
      if (!map) {
        fprintf (stderr, "Error: Creating the wordmap failed!\n");
        exit (1);
      }
      if ((subset_method == RAND_UNIQUE) && (subset_size > map->header->nwords)) {
        fprintf (stderr, "Error: Unique subset size (%llu) is bigger than number of unique kmers (%llu)\n", subset_size, map->header->nwords);
        exit (1);
      }
      sprintf (c, "%s_subset", outputname);
      v = subset (map, subset_method, subset_size, c);
      exit (0);
    } else {
      fprintf(stderr, "Error: Missing one or both list files!\n");
      print_help (1);
    }
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
    GT4WordMap *map1, *map2;
    int cinf;
    map1 = gt4_word_map_new (fnames[0], VERSION_MAJOR, USE_SCOUTS);
    map2 = gt4_word_map_new (fnames[1], VERSION_MAJOR, USE_SCOUTS);
    if (!map1 || !map2) {
      fprintf (stderr, "Error: Creating the wordmap failed!\n");
      exit (1);
    }
    if (map1->header->version_major > VERSION_MAJOR || map1->header->version_minor > VERSION_MINOR) {
      fprintf (stderr, "Error: %s is created with a newer glistmaker version.\n", map1->filename);
      exit (1);
    }  
    if (map2->header->version_major > VERSION_MAJOR || map2->header->version_minor > VERSION_MINOR) {
      fprintf (stderr, "Error: %s is created with a newer glistmaker version.\n", map2->filename);
      exit (1);
    }  
    cinf = compare_word_map_headers (map1->header, map2->header);
    if (cinf == -1) {
      fprintf (stdout, "Error: %s is not a glistmaker v.4 list!\n", map1->filename);
      exit (1);
    } else if (cinf == -2) {
      fprintf (stdout, "Error: %s is not a glistmaker v.4 list!\n", map2->filename);
      exit (1);
    } else if (cinf) {
      print_error_message (cinf);
      exit (1);
    }
    v = compare_wordmaps_mm (map1, map2, find_diff, find_ddiff, subtraction, countonly, outputname, cutoff, nmm, rule);
    gt4_word_map_delete (map1);
    gt4_word_map_delete (map2);
  } else if (nfiles == 2) {
    AZObject *list1, *list2;
    if (stream) {
      list1 = (AZObject *) gt4_word_list_stream_new (fnames[0], VERSION_MAJOR);
      list2 = (AZObject *) gt4_word_list_stream_new (fnames[1], VERSION_MAJOR);
    } else {
      list1 = (AZObject *) gt4_word_map_new (fnames[0], VERSION_MAJOR, use_scouts);
      list2 = (AZObject *) gt4_word_map_new (fnames[1], VERSION_MAJOR, use_scouts);
    }
    if (!list1 || !list2) {
      fprintf (stderr, "Error: Creating the wordmap failed!\n");
      exit (1);
    }
    v = compare_wordmaps (AZ_OBJECT(list1), AZ_OBJECT(list2), find_union, find_intrsec, find_diff, find_ddiff, subtraction, countonly, outputname, cutoff, rule);
    az_object_shutdown (AZ_OBJECT (list1));
    az_object_shutdown (AZ_OBJECT (list2));
  } else {
    AZObject *maps[MAX_FILES];
    GT4WordSArrayInstance *inst;
    GT4ListHeader header;
    unsigned int i;
    char name[2048];
    int ofile = 0;
    for (i = 0; i < nfiles; i++) {
      if (stream) {
        maps[i] = (AZObject *) gt4_word_list_stream_new (fnames[i], VERSION_MAJOR);
        if (!maps[i]) exit (1);
      } else {
        maps[i] = (AZObject *) gt4_word_map_new (fnames[i], VERSION_MAJOR, use_scouts);
        if (!maps[i]) exit (1);
      }
    }
    az_object_get_interface (maps[0], GT4_TYPE_WORD_SARRAY, (void **) &inst);
    if (!countonly) {
      sprintf (name, "%s_%d_union.list", outputname, inst->word_length);
      ofile = creat (name, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
      if (ofile < 0) {
        fprintf (stderr, "Creating output file %s failed\n", name);
        exit (1);
      }
    }
    v = union_multi (maps, nfiles, cutoff, ofile, &header);
    if (ofile > 0) close (ofile);
    fprintf (stdout, "NUnique\t%llu\nNTotal\t%llu\n", header.nwords, header.totalfreq);
  }
  if (v) return print_error_message (v);

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

static int
compare_wordmaps (AZObject *list1, AZObject *list2, int find_union, int find_intrsec, int find_diff, int find_ddiff, int subtract, int countonly, const char *out, unsigned int cutoff, int rule)
{
  GT4WordSArrayImplementation *impl1, *impl2;
  GT4WordSArrayInstance *inst1, *inst2;
  FILE *outf[4] = { 0 };
  GT4ListHeader h_out;

  unsigned long long word1, word2;
  unsigned int freq1, freq2;
  unsigned long long c_union = 0L, c_inters = 0L, c_diff1 = 0L, c_diff2 = 0L;
  unsigned long long freqsum_union = 0L, freqsum_inters = 0L, freqsum_diff1 = 0L, freqsum_diff2 = 0L;
  char fname[256]; /* the length is limited in main(..) method */

  impl1 = (GT4WordSArrayImplementation *) az_object_get_interface (list1, GT4_TYPE_WORD_SARRAY, (void **) &inst1);
  impl2 = (GT4WordSArrayImplementation *) az_object_get_interface (list2, GT4_TYPE_WORD_SARRAY, (void **) &inst2);

  if (inst1->word_length != inst2->word_length) {
    fprintf (stderr, "compare_wordmaps: Word lengths of lists differ: %u != %u\n", inst1->word_length, inst2->word_length);
    return 1;
  }
  
  if (debug) {
     fprintf (stderr, "compare_wordmaps: List 1: %llu entries\n", inst1->num_words);
     fprintf (stderr, "compare_wordmaps; List 2: %llu entries\n", inst2->num_words);
  }

  /* creating output files */
  if (find_union && !countonly) {
    sprintf (fname, "%s_%d_union.list", out, inst1->word_length);
    outf[0] = fopen (fname, "w");
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[0]);
  }
  if (find_intrsec && !countonly) {
    sprintf (fname, "%s_%d_intrsec.list", out, inst1->word_length);
    outf[1] = fopen (fname, "w");
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[1]);
  }
  if (find_diff && !countonly) {
    sprintf (fname, "%s_%d_0_diff1.list", out, inst1->word_length);
    outf[2] = fopen (fname, "w");
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[2]);
  }
  if (find_ddiff && !countonly) {
    sprintf (fname, "%s_%d_0_diff2.list", out, inst1->word_length);
    outf[3] = fopen (fname, "w");
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[3]);
  }

  gt4_word_sarray_get_first_word (impl1, inst1);
  word1 = inst1->word;
  freq1 = inst1->count;
  gt4_word_sarray_get_first_word (impl2, inst2);
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
      gt4_word_sarray_get_next_word (impl1, inst1);
      word1 = inst1->word;
      freq1 = inst1->count;
      gt4_word_sarray_get_next_word (impl2, inst2);
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
      gt4_word_sarray_get_next_word (impl1, inst1);
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
      gt4_word_sarray_get_next_word (impl2, inst2);
      word2 = inst2->word;
      freq2 = inst2->count;
    }
  }

  h_out.code = GT4_LIST_CODE;
  h_out.version_major = VERSION_MAJOR;
  h_out.version_minor = VERSION_MINOR;
  h_out.wordlength = inst1->word_length;

  /* add headers and close files */
  if (find_union && !countonly) {
    h_out.nwords = c_union;
    h_out.totalfreq = freqsum_union;
    fseek (outf[0], 0, SEEK_SET);
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[0]);
    fclose (outf[0]);
  } else if (find_union) {
    fprintf (stdout, "NUnique\t%llu\nNTotal\t%llu\n", c_union, freqsum_union);
  }
  if (find_intrsec && !countonly) {
    h_out.nwords = c_inters;
    h_out.totalfreq = freqsum_inters;
    fseek (outf[1], 0, SEEK_SET);
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[1]);
    fclose (outf[1]);
  } else if (find_intrsec) {
    fprintf (stdout, "NUnique\t%llu\nNTotal\t%llu\n", c_inters, freqsum_inters);
  }
  if (find_diff && !countonly) {
    h_out.nwords = c_diff1;
    h_out.totalfreq = freqsum_diff1;
    fseek (outf[2], 0, SEEK_SET);
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[2]);
    fclose (outf[2]);
  } else if (find_diff) {
    fprintf (stdout, "NUnique\t%llu\nNTotal\t%llu\n", c_diff1, freqsum_diff1);
  }
  if (find_ddiff && !countonly) {
    h_out.nwords = c_diff2;
    h_out.totalfreq = freqsum_diff2;
    fseek (outf[3], 0, SEEK_SET);
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[3]);
    fclose (outf[3]);
  } else if (find_ddiff) {
    fprintf (stdout, "NUnique\t%llu\nNTotal\t%llu\n", c_diff2, freqsum_diff2);
  }
  return 0;
}

#define TMP_BUF_SIZE (256 * 12)

static unsigned int
union_multi (AZObject *m[], unsigned int nmaps, unsigned int cutoff, int ofile, GT4ListHeader *header)
{
  GT4WordSArrayImplementation *impls[MAX_FILES];
  GT4WordSArrayInstance *insts[MAX_FILES];
  unsigned int n_sources;
  unsigned int j;
  unsigned long long word;
  unsigned char b[TMP_BUF_SIZE];
  unsigned int bp = 0;
  unsigned long long total = 0;
  double t_s, t_e;
  
  n_sources = 0;
  for (j = 0; j < nmaps; j++) {
    impls[n_sources] = (GT4WordSArrayImplementation *) az_object_get_interface (AZ_OBJECT(m[j]), GT4_TYPE_WORD_SARRAY, (void **) &insts[n_sources]);
    if (insts[n_sources]->num_words) {
      gt4_word_sarray_get_first_word (impls[n_sources], insts[n_sources]);
      total += insts[n_sources]->num_words;
      n_sources += 1;
    }
  }

  header->code = GT4_LIST_CODE;
  header->version_major = VERSION_MAJOR;
  header->version_minor = VERSION_MINOR;
  header->wordlength = insts[0]->word_length;
  header->nwords = 0;
  header->totalfreq = 0;
  header->list_start = sizeof (GT4ListHeader);

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
        if (!gt4_word_sarray_get_next_word (impls[j], insts[j])) {
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
      header->nwords += 1;
      header->totalfreq += freq;
      if (debug && !(header->nwords % 100000000)) {
        fprintf (stderr, "Words written: %uM\n", (unsigned int) (header->nwords / 1000000));
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
    fprintf (stderr, "Combined %u maps: input %llu (%.3f Mwords/s) output %llu (%.3f Mwords/s)\n", nmaps, total, total / (1000000 * (t_e - t_s)), header->nwords, header->nwords / (1000000 * (t_e - t_s)));
  }
  
  return 0;
}

static unsigned int
subset (GT4WordMap *map, unsigned int subset_method, unsigned long long subset_size, const char *filename)
{
  wordtable *wt;
  unsigned long long i;

  wt = wordtable_new (map->header->wordlength, subset_size);
  /* fixme: We rely here on RNG cycle being at least >>32 bits */
  if (subset_method == RAND_ALL) {
    while (subset_size > 0) {
      /* Pick random KMer */
      unsigned long long lhs = (unsigned long long) (((rand () + 1.0) / RAND_MAX) * 0xffffffff);
      unsigned long long rhs = (unsigned long long) (((rand () + 1.0) / RAND_MAX) * 0xffffffff);
      unsigned long long p = ((lhs << 32) | rhs) % map->header->nwords;
      wordtable_add_word_nofreq (wt, WORDMAP_WORD(map, p), map->header->wordlength);
      subset_size -= 1;
    }
  } else if (subset_method == RAND_UNIQUE) {
    unsigned long long *q = (unsigned long long *) malloc (map->header->nwords * sizeof (unsigned long long));
    for (i = 0; i < map->header->nwords; i++) q[i] = i;
    for (i = 0; i < subset_size; i++) {
      unsigned long long lhs = (unsigned long long) (((rand () + 1.0) / RAND_MAX) * 0xffffffff);
      unsigned long long rhs = (unsigned long long) (((rand () + 1.0) / RAND_MAX) * 0xffffffff);
      unsigned long long p = ((lhs << 32) | rhs) % map->header->nwords;
      unsigned long long t = q[i];
      q[i] = q[p];
      q[p] = t;
    }
    for (i = 0; i < subset_size; i++) {
      wordtable_add_word_nofreq (wt, WORDMAP_WORD(map, q[i]), map->header->wordlength);
    }
    free (q);
  }
  wordtable_sort (wt, 0);
  wordtable_find_frequencies (wt);
  wordtable_write_to_file (wt, filename, 1);
  wordtable_delete (wt);
  return 0;
}

static int
compare_wordmaps_mm (GT4WordMap *map1, GT4WordMap *map2, int find_diff, int find_ddiff, int subtract, int countonly, const char *out, unsigned int cutoff, unsigned int nmm, int rule)
{

  FILE *outf[4] = { 0 };
  GT4ListHeader h_out;

  unsigned long long i = 0L, j = 0L, word1, word2;
  unsigned int freq1, freq2;
  unsigned long long c_diff1 = 0L, c_diff2 = 0L;
  unsigned long long freqsum_diff1 = 0L, freqsum_diff2 = 0L;
  int cinf, v;
  char fname[256]; /* the length is limited in main(..) method */

  /* wordtables for -diff and -ddiff in case we want to use mismatches */
  wordtable tables[2];
  wordtable *difftable = &tables[0];
  wordtable *ddifftable = &tables[1];

  cinf = compare_word_map_headers (map1->header, map2->header);
  if (cinf == -1) {
    fprintf (stdout, "Error: %s is not a glistmaker v.4 list!\n", map1->filename);
    return 1;
  }
  if (cinf == -2) {
    fprintf (stdout, "Error: %s is not a glistmaker v.4 list!\n", map2->filename);
    return 1;
  }
  if (cinf) v = print_error_message (cinf);
  else v = 0;
  if (v) return v;

  /* filling wordtables */
  memset (difftable, 0, sizeof (wordtable));
  memset (ddifftable, 0, sizeof (wordtable));
  if (nmm && find_diff) {
    difftable->wordlength = map1->header->wordlength;
  }
  if (nmm && find_ddiff) {
    ddifftable->wordlength = map2->header->wordlength;
  }

  /* creating output files */
  if (find_diff && !countonly) {
    sprintf (fname, "%s_%d_%d_diff1.list", out, map1->header->wordlength, nmm);
    outf[2] = fopen (fname, "w");
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[2]);
  }
  if (find_ddiff && !countonly) {
    sprintf (fname, "%s_%d_%d_diff2.list", out, map1->header->wordlength, nmm);
    outf[3] = fopen (fname, "w");
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[3]);
  }

  word1 = *((unsigned long long *) (map1->wordlist + i * WORDMAP_ELEMENT_SIZE));
  word2 = *((unsigned long long *) (map2->wordlist + j * WORDMAP_ELEMENT_SIZE));
  freq1 = *((unsigned int *) (map1->wordlist + i * WORDMAP_ELEMENT_SIZE + sizeof (unsigned long long)));
  freq2 = *((unsigned int *) (map2->wordlist + j * WORDMAP_ELEMENT_SIZE + sizeof (unsigned long long)));

  if (debug) {
    fprintf (stderr, "Table 1: %llu entries\n", map1->header->nwords);
    fprintf (stderr, "Table 2: %llu entries\n", map2->header->nwords);
        }
  while (i < map1->header->nwords || j < map2->header->nwords) {
    unsigned int first_ge_cutoff = (freq1 >= cutoff);
    unsigned int second_ge_cutoff = (freq2 >= cutoff);

    if (word1 == word2) {
            /* Two words are equal */
      if (find_diff) {
        if (subtract) {
          if (freq1 <= freq2) {
            freq2 -= freq1;
          }
        }
        if (first_ge_cutoff && !second_ge_cutoff) {
          wordtable_add_word (difftable, word1, freq1 - freq2, map1->header->wordlength);
          c_diff1 += 1;
        }
      }
      if (find_ddiff && second_ge_cutoff && !first_ge_cutoff) {
        wordtable_add_word (ddifftable, word2, freq2 - freq1, map2->header->wordlength);
        c_diff2 += 1;
      }

    /* first word is smaller */
    } else if (word1 < word2) {
      if (find_diff && first_ge_cutoff && !subtract) {
        wordtable_add_word (difftable, word1, freq1, map1->header->wordlength);
        c_diff1 += 1;
        
      }
    /* second word is smaller */
    } else {
      if (find_ddiff && second_ge_cutoff) {
        wordtable_add_word (ddifftable, word2, freq2, map2->header->wordlength);
        c_diff2 += 1;
      }
    }
    /* advance relevant indices */
    if (word1 <= word2) {
      i += 1;
      if (i < map1->header->nwords) {
        word1 = *((unsigned long long *) (map1->wordlist + i * WORDMAP_ELEMENT_SIZE));
        freq1 = *((unsigned int *) (map1->wordlist + i * WORDMAP_ELEMENT_SIZE + sizeof (unsigned long long)));
      } else {
        word1 = ~0L;
        freq1 = 0;
      }
    }
    if (word2 <= word1) {
      j += 1;
      if (j < map2->header->nwords) {
        word2 = *((unsigned long long *) (map2->wordlist + j * WORDMAP_ELEMENT_SIZE));
        freq2 = *((unsigned int *) (map2->wordlist + j * WORDMAP_ELEMENT_SIZE + sizeof (unsigned long long)));
      } else {
        word2 = ~0L;
        freq2 = 0;
      }
    }
  }

  /* finding the mismatches */
  if (find_diff) {
    if (debug > 0) fprintf (stderr, "Finding diff with mismatches (%llu entries)\n", difftable->nwords);
    c_diff1 = fetch_relevant_words (difftable, map2, map1, cutoff, nmm, outf[2], subtract, countonly, &freqsum_diff1);
  }
  if (find_ddiff) {
    c_diff2 = fetch_relevant_words (ddifftable, map1, NULL, cutoff, nmm, outf[3], subtract, countonly, &freqsum_diff2);
  }

  h_out.code = GT4_LIST_CODE;
  h_out.version_major = VERSION_MAJOR;
  h_out.version_minor = VERSION_MINOR;
  h_out.wordlength = map1->header->wordlength;

  /* add headers and close files */
  if (find_diff && !countonly) {
    h_out.nwords = c_diff1;
    h_out.totalfreq = freqsum_diff1;
    fseek (outf[2], 0, SEEK_SET);
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[2]);
    fclose (outf[2]);
  } else if (find_diff) {
    fprintf (stdout, "NUnique\t%llu\nNTotal\t%llu\n", c_diff1, freqsum_diff1);
  }
  if (find_ddiff && !countonly) {
    h_out.nwords = c_diff2;
    h_out.totalfreq = freqsum_diff2;
    fseek (outf[3], 0, SEEK_SET);
    fwrite (&h_out, sizeof (GT4ListHeader), 1, outf[3]);
    fclose (outf[3]);
  } else if (find_ddiff) {
    fprintf (stdout, "NUnique\t%llu\nNTotal\t%llu\n", c_diff2, freqsum_diff2);
  }
  return 0;
}

static unsigned long long
fetch_relevant_words (wordtable *table, GT4WordMap *map, GT4WordMap *querymap, unsigned int cutoff, unsigned int nmm, FILE *f, int subtract, int countonly, unsigned long long *totalfreq)
{
  parameters p = {0};
  unsigned long long ri, wi, word, sumfreq = 0L, count = 0L;
  unsigned int freq, cnmm;

  if (table->nwords == 0) return 0;
  p.wordlength = table->wordlength;

  for (cnmm = 1; cnmm <= nmm; cnmm++) {
    p.nmm = cnmm;
    wi = 0L;

    for (ri = 0L; ri < table->nwords; ri++) {
            if (debug > 2) {
              fprintf (stderr, "cnmm %u ri %llu wi %llu\n", cnmm, ri, wi);
            }
      word = table->words[ri];
      freq = table->frequencies[ri];      
      sumfreq = word_map_search_query (map, word, &p, 0, 1, subtract, querymap);
      
      
      if (cnmm == nmm && sumfreq < cutoff) {
        if (!countonly) write_word_to_file (word, freq, f);
        count += 1;
        *totalfreq += freq;
        
      } else if (sumfreq < cutoff) {
        table->words[wi] = word;
        table->frequencies[wi] = freq;
        wi += 1;
      }

    }
    table->nwords = wi;
  }
  return count;
}

static int
compare_word_map_headers (GT4ListHeader *h1, GT4ListHeader *h2)
{
  if (h1->code != GT4_LIST_CODE) return -1;
  if (h2->code != GT4_LIST_CODE) return -2;
  if (h1->wordlength != h2->wordlength) return GT_INCOMPATIBLE_WORDLENGTH_ERROR;
  if (h1->version_major != h2->version_major) return GT_INCOMPATIBLE_VERSION_WARNING;
  if (h1->version_minor != h2->version_minor) return GT_INCOMPATIBLE_VERSION_WARNING;
  return 0;
}

static void
print_help (int exit_value)
{
  fprintf (stdout, "glistcompare version %u.%u (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_QUALIFIER);
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
  fprintf (stdout, "    -ss, --subset METHOD SIZE - make subset with given method (rand, rand_unique)\n");
  fprintf (stdout, "    --count_only             - output count of k-mers instead of k-mers themself\n");
  fprintf (stdout, "    --disable_scouts         - disable list read-ahead in background thread\n");
  fprintf (stdout, "    --stream                 - read input as stream (do not memory map files)\n");
  fprintf (stdout, "    -D                       - increase debug level\n");
  exit (exit_value);
}
