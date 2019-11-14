/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 * 
 * Copyright (C) 2014 University of Tartu
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

#define MAX_LISTS 1024

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "common.h"
#include "file-array.h"
#include "utils.h"
#include "sequence.h"
#include "sequence-stream.h"
#include "fasta.h"
#include "version.h"
#include "index-map.h"
#include "word-list-stream.h"
#include "word-map.h"
#include "word-table.h"

struct QueryData {
  GT4WordDictImplementation *impl;
  GT4WordDictInstance *inst;
  unsigned int n_mm;
  unsigned int pm_3;
  unsigned int min_freq;
  unsigned int max_freq;
  int print_all;
};

static void search_one_query_string (GT4WordDictImplementation *impl, GT4WordDictInstance *inst, const char *querystring, unsigned int n_mm, unsigned int pm_3, unsigned int min_freq, unsigned int max_freq, int print_all);
static unsigned int search_n_query_strings (GT4WordDictImplementation *impl, GT4WordDictInstance *inst, const char *fname, unsigned int n_mm, unsigned int pm_3, unsigned int minfreq, unsigned int maxfreq, int printall);
static unsigned int search_fasta (GT4WordDictImplementation *impl, GT4WordDictInstance *inst, const char *fname, unsigned int n_mm, unsigned int pm_3, unsigned int minfreq, unsigned int maxfreq, int printall);
static unsigned int search_list (GT4WordDictImplementation *impl, GT4WordDictInstance *inst, const char *querylistfilename, unsigned int n_mm, unsigned int pm_3, unsigned int min_freq, unsigned int max_freq, int print_all);
unsigned int search_lists_multi (const char *list, const char *lists[], unsigned int n_lists);
int process_word (GT4FastaReader *reader, unsigned long long word, void *data);
int print_full_map (AZObject *obj);
void get_statistics (AZObject *obj);
void print_median (AZObject *obj);
void print_distro (AZObject *obj, unsigned int size);
void print_gc (AZObject *obj);
static int print_files (AZObject *obj);
void print_help (int exitvalue);

int debug = 0;

unsigned int use_scouts = 1;

int main (int argc, const char *argv[])
{
  int argidx, v = 0, i;
  unsigned int n_lists = 0, invalid = 0;
  const char *lists[MAX_LISTS];
  const char *querystring = NULL, *queryfilename = NULL, *seqfilename = NULL, *querylistfilename = NULL;
  parameters p = {0};
  char *end;
  int printall = 0, getstat = 0, getmed = 0;
  unsigned int minfreq = 0, maxfreq = UINT_MAX;
  unsigned int distro = 0;
  unsigned int gc = 0;
  
  for (argidx = 1; argidx < argc; argidx++) {
    if (!strcmp (argv[argidx], "-v") || !strcmp (argv[argidx], "--version")) {
      fprintf (stdout, "glistquery version %u.%u.%u (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_MICRO, VERSION_QUALIFIER);
      return 0;

    } else if (!strcmp (argv[argidx], "-h") || !strcmp (argv[argidx], "--help") || !strcmp (argv[argidx], "-?")) {
      print_help (0);

    } else if (!strcmp(argv[argidx], "-s") || !strcmp(argv[argidx], "--seqfile")) {
      if (!argv[argidx + 1] || argv[argidx + 1][0] == '-') {
        fprintf(stderr, "Warning: No sequence file name specified!\n");
        argidx += 1;
        continue;
      }
      seqfilename = argv[argidx + 1];
      argidx += 1;
    } else if (!strcmp(argv[argidx], "-l") || !strcmp(argv[argidx], "--listfile")) {
      if (!argv[argidx + 1] || argv[argidx + 1][0] == '-') {
        fprintf(stderr, "Warning: No query list file name specified!\n");
        argidx += 1;
        continue;
      }
      querylistfilename = argv[argidx + 1];
      argidx += 1;
    } else if (!strcmp(argv[argidx], "-f") || !strcmp(argv[argidx], "--queryfile")) {
      if (!argv[argidx + 1] || argv[argidx + 1][0] == '-') {
        fprintf(stderr, "Warning: No query file name specified!\n");
        argidx += 1;
        continue;
      }
      queryfilename = argv[argidx + 1];
      argidx += 1;
    } else if (!strcmp(argv[argidx], "-q") || !strcmp(argv[argidx], "--query")) {
      if (!argv[argidx + 1] || argv[argidx + 1][0] == '-') {
        fprintf(stderr, "Warning: No query specified!\n");
        argidx += 1;
        continue;
      }
      querystring = argv[argidx + 1];
      argidx += 1;
    } else if (!strcmp(argv[argidx], "-p") || !strcmp(argv[argidx], "--perfectmatch")) {
      if (!argv[argidx + 1]) {
        fprintf(stderr, "Warning: No number of 3 prime perfect matches specified! Using the default value: %d.\n", p.pm3);
        argidx += 1;
        continue;
      }
      p.pm3 = strtol (argv[argidx + 1], &end, 10);
      if (*end != 0) {
        fprintf(stderr, "Error: Invalid number of 3 prime perfect matches: %s! Must be an integer.\n", argv[argidx + 1]);
        print_help (1);
      }
      argidx += 1;
    } else if (!strcmp(argv[argidx], "-mm") || !strcmp(argv[argidx], "--mismatch")) {
      if (!argv[argidx + 1]) {
        fprintf(stderr, "Warning: No number of mismatches specified! Using the default value: %d.\n", p.nmm);
        argidx += 1;
        continue;
      }
      p.nmm = strtol (argv[argidx + 1], &end, 10);
      if (*end != 0) {
        fprintf(stderr, "Error: Invalid number of mismatches: %s! Must be an integer.\n", argv[argidx + 1]);
        print_help (1);
      }
      argidx += 1;
    } else if (!strcmp(argv[argidx], "-min") || !strcmp(argv[argidx], "--minfreq")) { 
      if (!argv[argidx + 1]) {
        fprintf(stderr, "Warning: No minimum frequency specified! Using the default value: %d.\n", minfreq);
        argidx += 1;
        continue;
      }
      minfreq = strtol (argv[argidx + 1], &end, 10);
      if (*end != 0) {
        fprintf(stderr, "Error: Invalid minimum frequency: %s! Must be a positive integer.\n", argv[argidx + 1]);
        print_help (1);
      }
      argidx += 1;
    } else if (!strcmp(argv[argidx], "-max") || !strcmp(argv[argidx], "--maxfreq")) {
      if (!argv[argidx + 1]) {
        fprintf(stderr, "Warning: No maximum frequency specified! Using the default value: %d.\n", maxfreq);
        argidx += 1;
        continue;
      }
      maxfreq = strtol (argv[argidx + 1], &end, 10);
      if (*end != 0) {
        fprintf(stderr, "Error: Invalid maximum frequency: %s! Must be a positive integer.\n", argv[argidx + 1]);
        print_help (1);
      }
      argidx += 1;
    } else if (!strcmp(argv[argidx], "-D")) {
      debug += 1;
    } else if (!strcmp(argv[argidx], "--all") || !strcmp(argv[argidx], "-all")) {
      printall = 1;
    } else if (!strcmp(argv[argidx], "--stat") || !strcmp(argv[argidx], "-stat")) {  
      getstat = 1;
    } else if (!strcmp(argv[argidx], "--median") || !strcmp(argv[argidx], "-median")) {
      getmed = 1;
    } else if (!strcmp(argv[argidx], "--distribution") || !strcmp(argv[argidx], "-distribution")) {
      if ((argidx + 1) >= argc) {
        print_help (1);
      }
      argidx += 1;
      distro = strtol (argv[argidx], &end, 10);;
    } else if (!strcmp(argv[argidx], "-gc") || !strcmp(argv[argidx], "--gc")) {
      gc = 1;
    } else if (!strcmp(argv[argidx], "--disable_scouts")) {  
      use_scouts = 0;
    } else if (argv[argidx][0] != '-') {
      lists[n_lists++] = argv[argidx];
    } else {
      fprintf(stderr, "Error: Unknown argument: %s!\n", argv[argidx]);
      print_help (1);
    }
  }
  
  if (!n_lists) {
    fprintf(stderr, "Error: Missing a list file!\n");
    print_help (1);      
  }

  /* checking the parameters */
  if (p.nmm < 0) {
    fprintf(stderr, "Error: Invalid number of mismatches: %d! Must be between 0 and the length of the query.\n", p.nmm);
    print_help (1);
  }
  if (p.pm3 < 0) {
    fprintf(stderr, "Error: Invalid number of 3 prime perfect matches: %d! Must be between 0 and the length of the query.\n", p.pm3);
    print_help (1);
  }
  if (minfreq < 0) {
    fprintf(stderr, "Error: Invalid number of minimum frequency: %d! Must be positive integer.\n", minfreq);
    print_help (1);    
  }
  if (maxfreq < 0) {
    fprintf(stderr, "Error: Invalid number of maximum frequency: %d! Must be positive integer.\n", maxfreq);
    print_help (1);    
  }

  if (querylistfilename && n_lists > 1) {
    unsigned int result = search_lists_multi (querylistfilename, lists, n_lists);
    exit (result);
  }

  AZObject *maps[MAX_LISTS + 1];
  for (i = 0; i < n_lists; i++) {
    FILE *ifs;
    unsigned int code = 0;
    ifs = fopen (lists[i], "r");
    if (!ifs) {
      fprintf (stderr, "Cannot open list %s\n", lists[i]);
      exit (1);
    }
    fread (&code, 4, 1, ifs);
    fclose (ifs);
    if (code == GT4_LIST_CODE) {
      maps[i] = (AZObject *) gt4_word_map_new (lists[i], VERSION_MAJOR, !getstat && use_scouts);
    } else if (code == GT4_INDEX_CODE) {
      maps[i] = (AZObject *) gt4_index_map_new (lists[i], VERSION_MAJOR, !getstat && use_scouts);
    } else {
      fprintf (stderr, "Error: %s is not a valid GenomeTester4 list/index file\n", lists[i]);
      invalid = 1;
    }
  }
  if (invalid) exit (1);
  /* --stat */
  if (getstat) {
    for (i = 0; i < n_lists; i++) {
      get_statistics (maps[i]);
    }
    exit (0);
  }
  /* --median */
  if (getmed) {
    for (i = 0; i < n_lists; i++) {
      print_median (maps[i]);
    }
    exit (0);
  }
  /* --distribution */
  if (distro) {
    for (i = 0; i < n_lists; i++) {
      print_distro (maps[i], distro + 1);
    }
    exit (0);
  }
  /* --gc */
  if (gc) {
    for (i = 0; i < n_lists; i++) {
      print_gc (maps[i]);
    }
    exit (0);
  }

  if (!seqfilename && !querylistfilename && !queryfilename && !querystring) {
    for (i = 0; i < n_lists; i++) {
      print_full_map (maps[i]);
    }
    exit (0);
  }

  if (!GT4_IS_WORD_MAP (maps[0])) {
    fprintf (stderr, "Error: Could not make wordmap from file %s!\n", lists[0]);
    return 1;
  }
  GT4WordMap *map = (GT4WordMap *) maps[0];

  /* glistquery options */
  if (seqfilename) {
    /* FastA input */
    GT4WordDictImplementation *impl;
    GT4WordDictInstance *inst;
    impl = (GT4WordDictImplementation *) az_object_get_interface (maps[0], GT4_TYPE_WORD_DICT, (void **) &inst);
    v = search_fasta (impl, inst, seqfilename, p.nmm, p.pm3, minfreq, maxfreq, printall);
    if (v) return v;
  } else if (querylistfilename) { /* list input */  
    GT4WordDictImplementation *impl;
    GT4WordDictInstance *inst;
    impl = (GT4WordDictImplementation *) az_object_get_interface (maps[0], GT4_TYPE_WORD_DICT, (void **) &inst);
    v = search_list (impl, inst, querylistfilename, p.nmm, p.pm3, minfreq, maxfreq, printall);
    if (v) return v;
  } else if (queryfilename) { /* list of queries */
    GT4WordDictImplementation *impl;
    GT4WordDictInstance *inst;
    impl = (GT4WordDictImplementation *) az_object_get_interface (maps[0], GT4_TYPE_WORD_DICT, (void **) &inst);
    v = search_n_query_strings (impl, inst, queryfilename, p.nmm, p.pm3, minfreq, maxfreq, printall);
    if (v) return v;
  } else if (querystring) { /* one query */
    GT4WordDictImplementation *impl;
    GT4WordDictInstance *inst;
    /* checking possible errors */
    if (map->header->word_length != strlen (querystring)) {
      fprintf (stderr, "Error: Incompatible wordlengths! Wordlength in list: %u, query length: %lu\n", map->header->word_length, strlen (querystring));
      return 1;
    }
    if (map->header->word_length - p.pm3 < p.nmm) {
      fprintf(stderr, "Error: Number or mismatches specified is too large for %s with %d nucleotides long 3 prime perfect match.\n", querystring, p.pm3);
      return 1;
    }
    impl = (GT4WordDictImplementation *) az_object_get_interface (maps[0], GT4_TYPE_WORD_DICT, (void **) &inst);
    search_one_query_string (impl, inst, querystring, p.nmm, p.pm3, 0, UINT_MAX, printall);
  }

  for (i = 0; i < n_lists; i++) {
    az_object_shutdown (maps[i]);
  }

  exit (0);
}

static int
print_files (AZObject *obj)
{
        GT4FileArrayImplementation *impl;
        GT4FileArrayInstance *inst;
        unsigned int i;
        impl = (GT4FileArrayImplementation *) az_object_get_interface (obj, GT4_TYPE_FILE_ARRAY, (void **) &inst);
        for (i = 0; i < inst->num_files; i++) {
                gt4_file_array_get_file (impl, inst, i);
                fprintf (stderr, "%u\t%s\t%llu\t%llu\n", i, inst->file_name, inst->file_size, inst->n_sequences);
        }
        return 0;
}

/* Print the whole list */
int
print_full_map (AZObject *obj)
{
        GT4WordSListImplementation *impl;
        GT4WordSListInstance *inst;
        impl = (GT4WordSListImplementation *) az_object_get_interface (obj, GT4_TYPE_WORD_SARRAY, (void **) &inst);
        gt4_word_slist_get_first_word (impl, inst);
        while (inst->idx < inst->num_words) {
                char b[64];
                word2string (b, inst->word, inst->word_length);
                fprintf (stdout, "%s\t%u\n", b, inst->count);
                gt4_word_slist_get_next_word (impl, inst);
        }
  return 0;
}

/* Single query */

static void
search_one_query_string (GT4WordDictImplementation *impl, GT4WordDictInstance *inst, const char *query, unsigned int n_mm, unsigned int pm_3, unsigned int min_freq, unsigned int max_freq, int print_all)
{
  unsigned long long word;
  word = string_to_word (query, inst->word_length);
  if (gt4_word_dict_lookup_mm (impl, inst, word, n_mm, pm_3, print_all, 0)) {
    /* Found it */
    if (!print_all && (inst->value >= min_freq) && (inst->value <= max_freq)) fprintf (stdout, "%s\t%u\n", query, inst->value);
  }
}

/* Text file, one query per line */

static unsigned int
search_n_query_strings (GT4WordDictImplementation *impl, GT4WordDictInstance *inst, const char *queryfile, unsigned int n_mm, unsigned int pm_3, unsigned int min_freq, unsigned int max_freq, int print_all)
{
  FILE *ifs;
  char querystring[256];
  int beg = 1;

  ifs = fopen (queryfile, "r");
  if (ifs == NULL) {
    fprintf (stderr, "Error: Cannot open file %s.\n", queryfile);
    return 1;
  }

  while (fscanf (ifs, "%s\n", querystring) != EOF) {
    if (beg) {
      /* Checking possible errors */
      if (inst->word_length != strlen (querystring)) {
        fprintf (stderr, "Error: Incompatible wordlengths! Wordlength in list: %u, query length: %lu\n", inst->word_length, strlen (querystring));
        return 1;
      }
      if ((inst->word_length - pm_3) < n_mm) {
        fprintf (stderr, "Error: Number or mismatches specified is too large for %s with %d nucleotides long 3 prime perfect match.\n", querystring, pm_3);
        return 1;
      }
      beg = 0;
    }
    search_one_query_string (impl, inst, querystring, n_mm, pm_3, min_freq, max_freq, print_all);
  }
  fclose (ifs);
  return 0;
}

static unsigned int
search_fasta (GT4WordDictImplementation *impl, GT4WordDictInstance *inst, const char *fname, unsigned int n_mm, unsigned int pm_3, unsigned int min_freq, unsigned int max_freq, int print_all)
{
  GT4SequenceStream *stream;
  struct QueryData qs = {0};
  int result;
  GT4FastaReader r;

  qs.impl = impl;
  qs.inst = inst;
  qs.n_mm = n_mm;
  qs.pm_3 = pm_3;
  qs.min_freq = min_freq;
  qs.max_freq = max_freq;
  qs.print_all = print_all;

  stream = gt4_sequence_stream_new (fname);
  if (!stream) {
    fprintf (stderr, "search_fasta: Cannot open %s\n", fname);
    return 1;
  }
  fasta_reader_init (&r, fname, inst->word_length, 0, GT4_SEQUENCE_STREAM_SEQUENCE_SOURCE_IMPLEMENTATION(stream), &stream->source_instance);
  result = fasta_reader_read_nwords (&r, 1000000000000ULL, NULL, NULL, NULL, NULL, process_word, (void *) &qs);
  fasta_reader_release (&r);
  az_object_shutdown (AZ_OBJECT (stream));
  return result;
}

static unsigned int
search_list (GT4WordDictImplementation *impl, GT4WordDictInstance *inst, const char *fname, unsigned int n_mm, unsigned int pm_3, unsigned int min_freq, unsigned int max_freq, int print_all)
{
  GT4WordSListImplementation *s_impl;
  GT4WordSListInstance *s_inst;
  AZObject *obj = NULL;
  FILE *ifs;
  unsigned int code = 0;

  ifs = fopen (fname, "r");
  if (!ifs) {
    fprintf (stderr, "search_list: Cannot open list %s\n", fname);
    return 1;
  }
  fread (&code, 4, 1, ifs);
  fclose (ifs);
  if (code == GT4_LIST_CODE) {
    obj = (AZObject *) gt4_word_map_new (fname, VERSION_MAJOR, 0);
  } else if (code == GT4_INDEX_CODE) {
    obj = (AZObject *) gt4_index_map_new (fname, VERSION_MAJOR, 0);
  } else {
    fprintf (stderr, "search_list: Invalid file format\n");
    return 1;
  }
  s_impl = (GT4WordSListImplementation *) az_object_get_interface (obj, GT4_TYPE_WORD_SLIST, (void **) &s_inst);

  if (inst->word_length != s_inst->word_length) {
    az_object_shutdown (obj);
    return GT_INCOMPATIBLE_WORDLENGTH_ERROR;
  }
  
  gt4_word_slist_get_first_word (s_impl, s_inst);
  while (s_inst->idx < s_inst->num_words) {
    uint64_t word = s_inst->word;
    if (gt4_word_dict_lookup_mm (impl, inst, word, n_mm, pm_3, print_all, 0)) {
      if (!print_all && (inst->value >= min_freq) && (inst->value <= max_freq)) fprintf (stdout, "%s\t%u\n", word_to_string (word, inst->word_length), inst->value);
    }
    if (!gt4_word_slist_get_next_word (s_impl, s_inst)) break;
  }
  az_object_shutdown (obj);
  return 0;
}

unsigned int
search_lists_multi (const char *list, const char *lists[], unsigned int n_lists)
{
  AZObject *maps[MAX_LISTS + 1];
  GT4WordSListImplementation *impls[MAX_LISTS + 1];
  GT4WordSListInstance *insts[MAX_LISTS + 1];
  unsigned int i;

  maps[0] = (AZObject *) gt4_word_list_stream_new (list, VERSION_MAJOR);
  if (!maps[0]) exit (1);
  impls[0] = (GT4WordSListImplementation *) az_object_get_interface (AZ_OBJECT(maps[0]), GT4_TYPE_WORD_SLIST, (void **) &insts[0]);
  if (!impls[0]) exit (1);
  gt4_word_slist_get_first_word (impls[0], insts[0]);
  for (i = 0; i < n_lists; i++) {
    maps[i + 1] = (AZObject *) gt4_word_list_stream_new (lists[i], VERSION_MAJOR);
    if (!maps[i + 1]) exit (1);
    impls[i + 1] = (GT4WordSListImplementation *) az_object_get_interface (AZ_OBJECT(maps[i + 1]), GT4_TYPE_WORD_SLIST, (void **) &insts[i + 1]);
    if (!impls[i + 1]) exit (1);
    gt4_word_slist_get_first_word (impls[i + 1], insts[i + 1]);
  }

  while (insts[0]->idx < insts[0]->num_words) {
    unsigned long long word = insts[0]->word;
    unsigned int printed = 0;
    for (i = 1; i <= n_lists; i++) {
      while ((insts[i]->idx < insts[i]->num_words) && (insts[i]->word < word)) {
        if (!gt4_word_slist_get_next_word (impls[i], insts[i])) break;
      }
      if ((insts[i]->idx < insts[i]->num_words) && (insts[i]->word == word)) {
        if (!printed) {
          char b[64];
          word2string (b, word, insts[0]->word_length);
          fprintf (stdout, "%s", b);
          printed = 1;
        }
        fprintf (stdout, "\t%u:%u", i - 1, insts[i]->count);
      }
    }
    if (printed) fprintf (stdout, "\n");
    gt4_word_slist_get_next_word (impls[0], insts[0]);
  }

  return 0;
}

int process_word (GT4FastaReader *reader, unsigned long long word, void *data)
{
  struct QueryData *qs = (struct QueryData *) data;
  if (gt4_word_dict_lookup_mm (qs->impl, qs->inst, word, qs->n_mm, qs->pm_3, qs->print_all, 0)) {
    if (!qs->print_all && (qs->inst->value >= qs->min_freq) && (qs->inst->value <= qs->max_freq)) fprintf (stdout, "%s\t%u\n", word_to_string (word, reader->wordlength), qs->inst->value);
  }
  return 0;
}

void get_statistics (AZObject *obj)
{
  GT4WordSListInstance *inst;
  az_object_get_interface (obj, GT4_TYPE_WORD_SLIST, (void **) &inst);
  if (GT4_IS_WORD_MAP (obj)) {
    GT4WordMap *map = GT4_WORD_MAP(obj);
    fprintf (stdout, "Statistics of list file %s <<Built with glistmaker version %d.%d>>\n", map->filename, map->header->version_major, map->header->version_minor);
  } else if (GT4_IS_INDEX_MAP (obj)) {
    GT4IndexMap *imap = GT4_INDEX_MAP(obj);
    fprintf (stdout, "Statistics of index file %s <<Built with glistmaker version %d.%d>>\n", imap->filename, imap->header->version_major, imap->header->version_minor);
  }
  fprintf (stdout, "Wordlength\t%u\n", inst->word_length);
  fprintf (stdout, "NUnique\t%llu\n", inst->num_words);
  fprintf (stdout, "NTotal\t%llu\n", inst->sum_counts);
  return;
}

void print_median (AZObject *obj)
{
  unsigned int min, max, med, gmin, gmax;
  unsigned long long i;
  GT4WordSArrayImplementation *impl;
  GT4WordSArrayInstance *inst;
  impl = (GT4WordSArrayImplementation *) az_object_get_interface (obj, GT4_TYPE_WORD_SARRAY, (void **) &inst);
  gmin = 0xffffffff;
  gmax = 0;
  if (debug > 0) fprintf (stderr, "Finding min/max...");
  for (i = 0; i < inst->slist_inst.num_words; i++) {
          gt4_word_sarray_get_word (impl, inst, i);
    if (inst->slist_inst.count < gmin) gmin = inst->slist_inst.count;
    if (inst->slist_inst.count > gmax) gmax = inst->slist_inst.count;
  }
  if (debug > 0) fprintf (stderr, "done (%u %u)\n", gmin, gmax);
  min = gmin;
  max = gmax;
  med = (unsigned int) (((unsigned long long) min + max) / 2);
  while (max > min) {
    unsigned long long above = 0, below = 0, equal;
    for (i = 0; i < inst->slist_inst.num_words; i++) {
            gt4_word_sarray_get_word (impl, inst, i);
      if (inst->slist_inst.count > med) above += 1;
      if (inst->slist_inst.count < med) below += 1;
    }
    equal = inst->slist_inst.num_words - above - below;
    if (debug > 0) fprintf (stderr, "Trying median %u - equal %llu, below %llu, above %llu\n", med, equal, below, above);
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
  
        if (GT4_IS_WORD_MAP (obj)) {
                GT4WordMap *map = GT4_WORD_MAP(obj);
          fprintf (stdout, "Statistics of list file %s <<Built with glistmaker version %d.%d>>\n", map->filename, map->header->version_major, map->header->version_minor);
        } else if (GT4_IS_INDEX_MAP (obj)) {
                GT4IndexMap *imap = GT4_INDEX_MAP(obj);
                fprintf (stdout, "Statistics of index file %s <<Built with glistmaker version %d.%d>>\n", imap->filename, imap->header->version_major, imap->header->version_minor);
        }
  fprintf (stdout, "Wordlength\t%u\n", inst->slist_inst.word_length);
  fprintf (stdout, "NUnique\t%llu\n", inst->slist_inst.num_words);
  fprintf (stdout, "NTotal\t%llu\n", inst->slist_inst.sum_counts);
  fprintf (stdout, "Min %u Max %u Median %u Average %.2f\n", gmin, gmax, med, (double) inst->slist_inst.sum_counts / inst->slist_inst.num_words);
}

void
print_distro (AZObject *obj, unsigned int max)
{
  unsigned long long *d;
  GT4WordSListImplementation *impl;
  GT4WordSListInstance *inst;
  unsigned int i;
  impl = (GT4WordSListImplementation *) az_object_get_interface (obj, GT4_TYPE_WORD_SLIST, (void **) &inst);
  d = (unsigned long long *) malloc (max * 8);
  memset (d, 0, max * 8);
  gt4_word_slist_get_first_word (impl, inst);
  while (inst->idx < inst->num_words) {
    if (inst->count <= max) d[inst->count - 1] += 1;
    gt4_word_slist_get_next_word (impl, inst);
  }
  for (i = 0; i < max; i++) {
    fprintf (stdout, "%u\t%llu\n", i + 1, d[i]);
  }
}

void
print_gc (AZObject *obj)
{
  unsigned long long count = 0;
  GT4WordSListImplementation *impl;
  GT4WordSListInstance *inst;
  impl = (GT4WordSListImplementation *) az_object_get_interface (obj, GT4_TYPE_WORD_SLIST, (void **) &inst);
  gt4_word_slist_get_first_word (impl, inst);
  while (inst->idx < inst->num_words) {
    unsigned long long word = inst->word;
    unsigned int j;
    for (j = 0; j < inst->word_length; j++) {
      count += inst->count * ((word ^ (word >> 1)) & 1);
      word = word >> 2;
    }
    gt4_word_slist_get_next_word (impl, inst);
  }
  printf ("GC\t%g\n", (double) count / (inst->num_words * inst->word_length));
}

void print_help (int exit_value)
{
    fprintf (stderr, "glistquery version %u.%u.%u (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_MICRO, VERSION_QUALIFIER);
  fprintf (stderr, "Usage: glistquery INPUT_LIST [OPTIONS]\n");
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "    -v, --version             - print version information and exit\n");
  fprintf (stderr, "    -h, --help                - print this usage screen and exit\n");
  fprintf (stderr, "    --stat                    - print statistics of the list file and exit\n");
  fprintf (stderr, "    --median                  - print min/max/median/average and exit\n");
  fprintf (stderr, "    --distribution MAX        - print distribution up to MAX\n");
  fprintf (stderr, "    --gc                      - print average GC content of all words\n");
  fprintf (stderr, "    -q, --query               - single query word\n");
  fprintf (stderr, "    -f, --queryfile           - list of query words in a file\n");
  fprintf (stderr, "    -s, --seqfile             - FastA/FastQ file\n");
  fprintf (stderr, "    -l, --listfile            - list file made by glistmaker\n");
  fprintf (stderr, "    -mm, --mismatch NUMBER    - specify number of mismatches (default 0)\n");
  fprintf (stderr, "    -p, --perfectmatch NUMBER - specify number of 3' perfect matches (default 0)\n");
  fprintf (stderr, "    -min, --minfreq NUMBER    - minimum frequency of the printed words (default 0)\n");
  fprintf (stderr, "    -max, --maxfreq NUMBER    - maximum frequency of the printed words (default MAX_UINT)\n");
  fprintf (stderr, "    --all                     - in case of mismatches prints all found words\n");
  fprintf (stderr, "    -D                        - increase debug level\n");
  exit (exit_value);
}


