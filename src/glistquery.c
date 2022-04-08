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
#include "set-operations.h"
#include "fasta.h"
#include "version.h"
#include "index-map.h"
#include "word-list-stream.h"
#include "word-map.h"
#include "word-table.h"

static unsigned int search_one_query_string (AZObject *obj, const char *query, unsigned int n_mm, unsigned int pm_3, unsigned int min_freq, unsigned int max_freq, int print_all_words);
static unsigned int search_n_query_strings (AZObject *obj, const char *fname, unsigned int n_mm, unsigned int pm_3, unsigned int minfreq, unsigned int maxfreq, int print_all_words);
static unsigned int search_fasta (AZObject *obj, const char *fname, unsigned int n_mm, unsigned int pm_3, unsigned int minfreq, unsigned int maxfreq, int printall);
static unsigned int search_list (AZObject *obj, const char *querylistfilename, unsigned int n_mm, unsigned int pm_3, unsigned int min_freq, unsigned int max_freq, int print_all);
unsigned int search_lists_multi (AZObject *query, AZObject *lists[], unsigned int n_lists);
static unsigned int print_full_map (AZObject *obj, unsigned int locations);
void get_statistics (AZObject *obj);
void print_median (AZObject *obj);
void print_distro (AZObject *obj, unsigned int size);
void print_gc (AZObject *obj);
static void print_files (GT4IndexMap *imap);
static void print_sequences (GT4IndexMap *imap);
void print_help (int exitvalue);

int debug = 0;

static unsigned int use_scouts = 1;
static unsigned int locations = 0;
static unsigned int use_3p = 0;
static unsigned int use_5p = 0;

enum {
  QUERY,
  STATS,
  GC,
  MEDIAN,
  DISTRO,
  FILES,
  SEQUENCES
};

typedef struct _DumpData DumpData;

struct _DumpData {
  unsigned int n_lists;
  unsigned int wlen;
};

static unsigned int
dump_callback (uint64_t word, uint32_t *counts, void *data)
{
  DumpData *dd = (DumpData *) data;
  unsigned int i;
  fprintf (stdout, "%s", word_to_string (word, dd->wlen));
  for (i = 0; i < dd->n_lists; i++) {
    fprintf (stdout, "\t%u", counts[i]);
  }
  fprintf (stdout, "\n");
  return 0;
}

static void
dump_lists (AZObject *objs[], unsigned int n_objs, unsigned int wlen, unsigned int is_union)
{
  DumpData dd;
  dd.n_lists = n_objs;
  dd.wlen = wlen;
  if (is_union) {
    gt4_is_union (objs, n_objs, dump_callback, &dd);
  } else {
    gt4_union (objs, n_objs, dump_callback, &dd);
  }
}
  
int main (int argc, const char *argv[])
{
  int argidx, v = 0, i;
  unsigned int n_lists = 0, invalid = 0;
  const char *lists[MAX_LISTS];
  const char *querystring = NULL, *queryfilename = NULL, *seqfilename = NULL, *querylistfilename = NULL;
  unsigned int nmm = 0;
  unsigned int pm3 = 0;
  char *end;
  int printall = 0, print_header = 0;
  unsigned int minfreq = 0, maxfreq = UINT_MAX;
  unsigned int distro = 0;
  unsigned int bloom = 0;
  unsigned int command = QUERY;
  unsigned int is_union = 0;

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
      argidx += 1;
      if (argidx >= argc) {
        print_help (1);
        break;
      }
      pm3 = strtol (argv[argidx], &end, 10);
      if (*end || (pm3 < 0) || (pm3 > 32)) {
        print_help (1);
        break;
      }
    } else if (!strcmp(argv[argidx], "-mm") || !strcmp(argv[argidx], "--mismatch")) {
      argidx += 1;
      if (argidx >= argc) {
        print_help (1);
        break;
      }
      nmm = strtol (argv[argidx], &end, 10);
      if (*end || (nmm < 0) || (nmm > 16)) {
        print_help (1);
        break;
      }
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
    } else if (!strcmp(argv[argidx], "--stats") || !strcmp(argv[argidx], "--stat") || !strcmp(argv[argidx], "-stat")) {  
      command = STATS;
    } else if (!strcmp(argv[argidx], "--median") || !strcmp(argv[argidx], "-median")) {
      command = MEDIAN;
    } else if (!strcmp(argv[argidx], "--distribution") || !strcmp(argv[argidx], "-distribution")) {
      if ((argidx + 1) >= argc) {
        print_help (1);
      }
      argidx += 1;
      distro = strtol (argv[argidx], &end, 10);
      command = DISTRO;
    } else if (!strcmp(argv[argidx], "-gc") || !strcmp(argv[argidx], "--gc")) {
      command = GC;
    } else if (!strcmp(argv[argidx], "--files")) {
      command = FILES;
    } else if (!strcmp(argv[argidx], "--sequences")) {
      command = SEQUENCES;
    } else if (!strcmp(argv[argidx], "--locations")) {
      locations = 1;
    } else if (!strcmp(argv[argidx], "--3p")) {
      use_3p = 1;
    } else if (!strcmp(argv[argidx], "--5p")) {
      use_5p = 1;
    } else if (!strcmp(argv[argidx], "--header")) {
      print_header = 1;
    } else if (!strcmp(argv[argidx], "--bloom")) {
      bloom = 1;
    } else if (!strcmp(argv[argidx], "--is_union")) {  
      is_union = 1;
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
    fprintf(stderr, "No list/index files specified!\n");
    print_help (1);      
  }

  /* Checking the parameters */
  if (minfreq < 0) {
    fprintf(stderr, "Error: Invalid number of minimum frequency: %d! Must be positive integer.\n", minfreq);
    print_help (1);    
  }
  if (maxfreq < 0) {
    fprintf(stderr, "Error: Invalid number of maximum frequency: %d! Must be positive integer.\n", maxfreq);
    print_help (1);    
  }

  /* Load all list/index files into maps array and test some errors */
  unsigned int use_stream = 0;
  if (querylistfilename && n_lists > 1) {
    use_stream = 1;
  }
  AZObject *maps[MAX_LISTS + 1];
  unsigned int has_lists = 0;
  unsigned int wlen = 0;
  for (i = 0; i < n_lists; i++) {
    FILE *ifs;
    unsigned int code = 0;
    maps[i] = NULL;
    ifs = fopen (lists[i], "r");
    if (!ifs) {
      fprintf (stderr, "Cannot open list %s\n", lists[i]);
      exit (1);
    }
    fread (&code, 4, 1, ifs);
    fclose (ifs);
    if (code == GT4_LIST_CODE) {
      if (use_stream) {
        maps[i] = (AZObject *) gt4_word_list_stream_new (lists[i], VERSION_MAJOR);
      } else {
        maps[i] = (AZObject *) gt4_word_map_new (lists[i], VERSION_MAJOR, (command != STATS) && use_scouts, (command != STATS) && bloom);
      }
      if (debug) fprintf (stderr, "List %s loaded\n", lists[i]);
      has_lists = 1;
    } else if (code == GT4_INDEX_CODE) {
      maps[i] = (AZObject *) gt4_index_map_new (lists[i], VERSION_MAJOR, (command != STATS) && use_scouts);
      /* has_indices = 1; */
    } else {
      fprintf (stderr, "Error: %s is not a valid GenomeTester4 list/index file\n", lists[i]);
      invalid = 1;
    }
    if (!maps[i]) {
      fprintf (stderr, "Error: %s is invalid or corrupted\n", lists[i]);
      invalid = 1;
    } else {
      GT4WordSListInstance *inst;
      az_object_get_interface (maps[i], GT4_TYPE_WORD_SLIST, (void **) &inst);
      if (!wlen) {
        wlen = inst->word_length;
      } else {
        if (inst->word_length != wlen) {
          fprintf (stderr, "Error: %s has different word length %u (first list had %u)\n", lists[i], inst->word_length, wlen);
          invalid = 1;
        }
      }
    }
  }
  /* Load query list into stream */
  AZObject *query_map = NULL;
  if (querylistfilename) {
    query_map = (AZObject *) gt4_word_list_stream_new (querylistfilename, VERSION_MAJOR);
    if (!query_map) {
      fprintf (stderr, "Error: %s is invalid or corrupted\n", querylistfilename);
      invalid = 1;
    } else {
      GT4WordSListInstance *inst;
      az_object_get_interface (query_map, GT4_TYPE_WORD_SLIST, (void **) &inst);
      if (inst->word_length != wlen) {
        fprintf (stderr, "Error: %s has different word length %u (first list had %u)\n", querylistfilename, inst->word_length, wlen);
        invalid = 1;
      }
    }
  }

  if (invalid) exit (1);

  /* Generic methods */
  if (command == STATS) {
    for (i = 0; i < n_lists; i++) {
      get_statistics (maps[i]);
    }
    exit (0);
  } else if (command == MEDIAN) {
    for (i = 0; i < n_lists; i++) {
      print_median (maps[i]);
    }
    exit (0);
  } else if (command == DISTRO) {
    for (i = 0; i < n_lists; i++) {
      print_distro (maps[i], distro + 1);
    }
    exit (0);
  } else if (command == GC) {
    for (i = 0; i < n_lists; i++) {
      print_gc (maps[i]);
    }
    exit (0);
  }

  /* Index generic queries */
  if (command == FILES) {
    if (has_lists || (n_lists > 1)) {
      fprintf (stderr, "Error: Files can only be queried from single index\n");
      exit (1);
    }
    print_files ((GT4IndexMap *) maps[0]);
    exit (0);
  } else if (command == SEQUENCES) {
    if (has_lists || (n_lists > 1)) {
      fprintf (stderr, "Error: Sequences can only be queried from single index\n");
      exit (1);
    }
    print_sequences ((GT4IndexMap *) maps[0]);
    exit (0);
  }

  /* If no options is given print all lists/indices */
  if (!seqfilename && !querylistfilename && !queryfilename && !querystring) {
    if (n_lists > 1) {
      if (print_header) {
        fprintf (stdout, "KMER");
        for (i = 0; i < n_lists; i++) {
          fprintf (stdout, "\t%s", lists[i]);
        }
        fprintf (stdout, "\n");
      }
      dump_lists (maps, n_lists, wlen, is_union);
    } else {
      for (i = 0; i < n_lists; i++) {
        print_full_map (maps[i], locations);
      }
    }
    exit (0);
  }

  /* Search one list agains multiple */
  if (query_map && n_lists > 1) {
    if (nmm || pm3) {
      fprintf (stderr, "Error: Searching multiple lists is incompatible with mismatches\n");
      exit (1);
    }
    unsigned int result = search_lists_multi (query_map, maps, n_lists);
    exit (result);
  }

  if (n_lists > 1) {
    fprintf (stderr, "Error: Query is incompatible with multiple lists/indices\n");
    exit (1);
  }
  if ((nmm + pm3) > wlen) {
    fprintf(stderr, "Error: Number of mismatches (%u) and 3' perfect match (%u) are longer than word length %u\n", nmm, pm3, wlen);
    return 1;
  }

  if (querystring) {
    /* Single query */
    v = search_one_query_string (maps[0], querystring, nmm, pm3, minfreq, maxfreq, printall);
    if (v) return v;
  } else if (queryfilename) {
    /* Multiple queries in text file */
    v = search_n_query_strings (maps[0], queryfilename, nmm, pm3, minfreq, maxfreq, printall);
    if (v) return v;
  } else if (seqfilename) {
    /* FastA input */
    v = search_fasta (maps[0], seqfilename, nmm, pm3, minfreq, maxfreq, printall);
    if (v) return v;
  } else if (querylistfilename) { /* list input */  
    v = search_list (maps[0], querylistfilename, nmm, pm3, minfreq, maxfreq, printall);
    if (v) return v;
  }

  for (i = 0; i < n_lists; i++) {
    az_object_shutdown (maps[i]);
  }

  exit (0);
}

static void
print_files (GT4IndexMap *imap)
{
  unsigned int i;
  GT4FileArrayImplementation *impl = GT4_INDEX_MAP_FILE_ARRAY_IMPL(imap);
  GT4FileArrayInstance *inst = &imap->file_array_inst;
  for (i = 0; i < inst->num_files; i++) {
    gt4_file_array_get_file (impl, inst, i);
    fprintf (stdout, "%u\t%s\t%llu\t%llu\n", i, inst->file_name, inst->file_size, inst->n_sequences);
  }
}

static void
print_sequences (GT4IndexMap *imap)
{
  unsigned int i;
  GT4FileArrayImplementation *impl = GT4_INDEX_MAP_FILE_ARRAY_IMPL(imap);
  GT4FileArrayInstance *inst = &imap->file_array_inst;
  for (i = 0; i < inst->num_files; i++) {
    unsigned int j;
    gt4_file_array_get_file (impl, inst, i);
    for (j = 0; j < inst->n_sequences; j++) {
      unsigned char b[1024];
      gt4_file_array_get_sequence (impl, inst, j);
      gt4_index_map_get_sequence_name (imap, b, 1024, i, j);
      fprintf (stdout, "%u\t%u\t%s\t%llu\t%llu\t%llu\n",i, j, b, inst->name_pos, inst->seq_pos, inst->seq_len);
    }
  }
}

static void
print_index_info (GT4WordIndexImplementation *impl, GT4WordIndexInstance *inst, unsigned int reverse)
{
  unsigned int i;
  for (i = 0; i < inst->n_locations; i++) {
    gt4_word_index_get_location (impl, inst, i);
    fprintf (stdout, "%u\t%u\t%llu\t%u\n", inst->file_idx, inst->seq_idx, (unsigned long long) inst->pos, !inst->dir != !reverse);
  }
}

/* Print the whole list */

static unsigned int
print_full_map (AZObject *obj, unsigned int locations)
{
  if (!locations || !GT4_IS_INDEX_MAP (obj)) {
    GT4WordSListImplementation *impl;
    GT4WordSListInstance *inst;
    impl = (GT4WordSListImplementation *) az_object_get_interface (obj, GT4_TYPE_WORD_SLIST, (void **) &inst);
    gt4_word_slist_get_first_word (impl, inst);
    while (inst->idx < inst->num_words) {
      char b[64];
      word2string (b, inst->word, inst->word_length);
      fprintf (stdout, "%s\t%u\n", b, inst->count);
      gt4_word_slist_get_next_word (impl, inst);
    }
  } else {
    GT4IndexMap *imap = GT4_INDEX_MAP(obj);
    GT4WordIndexImplementation *impl;
    GT4WordIndexInstance *inst;
    impl = (GT4WordIndexImplementation *) az_object_get_interface ((AZObject *) imap, GT4_TYPE_WORD_INDEX, (void **) &inst);
    gt4_word_slist_get_first_word (GT4_INDEX_MAP_SLIST_IMPLEMENTATION(imap), &imap->sarray_inst.slist_inst);
    while (imap->sarray_inst.slist_inst.idx < imap->sarray_inst.slist_inst.num_words) {
      char b[64];
      word2string (b, imap->sarray_inst.slist_inst.word, imap->sarray_inst.slist_inst.word_length);
      fprintf (stdout, "%s\t%u\n", b, imap->sarray_inst.slist_inst.count);
      print_index_info (impl, inst, 0);
      gt4_word_slist_get_next_word (GT4_INDEX_MAP_SLIST_IMPLEMENTATION(imap), &imap->sarray_inst.slist_inst);
    }
  }
  return 0;
}

typedef struct _QueryData QueryData;
struct _QueryData {
  GT4WordDictImplementation *dict_impl;
  GT4WordDictInstance *dict_inst;
  GT4WordIndexImplementation *index_impl;
  GT4WordIndexInstance *index_inst;
  unsigned int n_mm;
  unsigned int pm_3;
  unsigned int min_freq;
  unsigned int max_freq;
  int print_all_words;
  unsigned int reverse;
};

static unsigned int
cb_print (unsigned long long word, unsigned int value, void *data)
{
  QueryData *qd = (QueryData *) data;
  if (!qd->index_impl) {
    fprintf (stdout, "%s\t%u\n", word_to_string (word, qd->dict_inst->word_length), value);
  } else {
    fprintf (stdout, "%s\t%u\t%u\n", word_to_string (word, qd->dict_inst->word_length), value, qd->reverse);
    print_index_info (qd->index_impl, qd->index_inst, qd->reverse);
  }
  return 1;
}

/* Single query */

static void
search_one_word (QueryData *qd, unsigned long long word)
{
  unsigned long long rword;
  rword = get_reverse_complement (word, qd->dict_inst->word_length);
  if (rword < word) {
    word = rword;
    qd->reverse = 1;
  }
  if (qd->index_impl || qd->print_all_words) {
    if (!gt4_word_dict_lookup_mm (qd->dict_impl, qd->dict_inst, word, qd->n_mm, qd->pm_3, 1, 0, cb_print, qd) && !qd->min_freq) {
      /* If min freq is 0 always print */
      fprintf (stdout, "%s\t0\n", word_to_string (word, qd->dict_inst->word_length));
    };
  } else {
    if (gt4_word_dict_lookup_mm (qd->dict_impl, qd->dict_inst, word, qd->n_mm, qd->pm_3, 1, 0, NULL, NULL)) {
      if ((qd->dict_inst->value >= qd->min_freq) && (qd->dict_inst->value <= qd->max_freq)) {
        fprintf (stdout, "%s\t%u\n", word_to_string (word, qd->dict_inst->word_length), qd->dict_inst->value);
      }
    } else if (!qd->min_freq) {
      /* If min freq is 0 always print */
      fprintf (stdout, "%s\t0\n", word_to_string (word, qd->dict_inst->word_length));
    }
  }
}

/* One query */
static unsigned int
search_one_query_string (AZObject *obj, const char *query, unsigned int n_mm, unsigned int pm_3, unsigned int min_freq, unsigned int max_freq, int print_all_words)
{
  QueryData qd = { 0 };
  unsigned long long word;
  unsigned int len;

  qd.dict_impl = (GT4WordDictImplementation *) az_object_get_interface (obj, GT4_TYPE_WORD_DICT, (void **) &qd.dict_inst);
  len = (unsigned int) strlen (query);
  if (len != qd.dict_inst->word_length) {
    if (len < qd.dict_inst->word_length) {
      fprintf (stderr, "search_one_query_string: Word too short (%u < %u)\n", qd.dict_inst->word_length, len);
      return 1;
    } else if (use_3p) {
      word = string_to_word (query + (len - qd.dict_inst->word_length), qd.dict_inst->word_length);
    } else if (use_5p) {
      word = string_to_word (query, qd.dict_inst->word_length);
    } else {
      fprintf (stderr, "search_one_query_string: Wrong query length (%u != %u) - use --3p or --5p\n", qd.dict_inst->word_length, len);
      return 1;
    }
  } else {
    word = string_to_word (query, qd.dict_inst->word_length);
  }
  if (GT4_IS_INDEX_MAP (obj) && locations) {
    qd.index_impl = (GT4WordIndexImplementation *) az_object_get_interface (obj, GT4_TYPE_WORD_INDEX, (void **) &qd.index_inst);
  }
  qd.n_mm = n_mm;
  qd.pm_3 = pm_3;
  qd.min_freq = min_freq;
  qd.max_freq = max_freq;
  qd.print_all_words = print_all_words;

  search_one_word (&qd, word);
  return 0;
}

/* Text file, one query per line */

static unsigned int
search_n_query_strings (AZObject *obj, const char *queryfile, unsigned int n_mm, unsigned int pm_3, unsigned int min_freq, unsigned int max_freq, int print_all_words)
{
  QueryData qd = { 0 };
  FILE *ifs;

  ifs = fopen (queryfile, "r");
  if (ifs == NULL) {
    fprintf (stderr, "search_n_query_strings: Cannot open file %s.\n", queryfile);
    return 1;
  }
  qd.dict_impl = (GT4WordDictImplementation *) az_object_get_interface (obj, GT4_TYPE_WORD_DICT, (void **) &qd.dict_inst);
  if (GT4_IS_INDEX_MAP (obj) && locations) {
    qd.index_impl = (GT4WordIndexImplementation *) az_object_get_interface (obj, GT4_TYPE_WORD_INDEX, (void **) &qd.index_inst);
  }
  qd.n_mm = n_mm;
  qd.pm_3 = pm_3;
  qd.min_freq = min_freq;
  qd.max_freq = max_freq;
  qd.print_all_words = print_all_words;

  int val = fgetc (ifs);
  while (val > 0) {
    char c[256];
    unsigned int i = 0, len;
    unsigned long long word;
    while ((val > 0) && (i < 255) && (val != '\n')) {
      c[i++] = (char) val;
      val = fgetc (ifs);
    }
    c[i] = 0;
    while ((val > 0) && (val != '\n')) val = fgetc (ifs);
    while ((val > 0) && (val < 'A')) val = fgetc (ifs);
    len = (unsigned int) strlen (c);
    if (len != qd.dict_inst->word_length) {
      if (len < qd.dict_inst->word_length) {
        fprintf (stderr, "search_n_query_strings: Word too short (%u < %u)\n", qd.dict_inst->word_length, len);
        return 1;
      } else if (use_3p) {
        word = string_to_word (c + (len - qd.dict_inst->word_length), qd.dict_inst->word_length);
      } else if (use_5p) {
         word = string_to_word (c, qd.dict_inst->word_length);
      } else {
        fprintf (stderr, "search_n_query_strings: Wrong query length (%u != %u) - use --3p or --5p\n", qd.dict_inst->word_length, len);
        return 1;
      }
    } else {
      word = string_to_word (c, qd.dict_inst->word_length);
    }
    search_one_word (&qd, word);
  }
  fclose (ifs);
  return 0;
}

static int
process_word (GT4FastaReader *reader, unsigned long long word, void *data)
{
  QueryData *qd = (QueryData *) data;
  search_one_word (qd, word);
  return 0;
}

static unsigned int
search_fasta (AZObject *obj, const char *fname, unsigned int n_mm, unsigned int pm_3, unsigned int min_freq, unsigned int max_freq, int print_all_words)
{
  GT4SequenceStream *stream;
  QueryData qd = {0};
  int result;
  GT4FastaReader r;

  qd.dict_impl = (GT4WordDictImplementation *) az_object_get_interface (obj, GT4_TYPE_WORD_DICT, (void **) &qd.dict_inst);
  if (GT4_IS_INDEX_MAP (obj) && locations) {
    qd.index_impl = (GT4WordIndexImplementation *) az_object_get_interface (obj, GT4_TYPE_WORD_INDEX, (void **) &qd.index_inst);
  }
  qd.n_mm = n_mm;
  qd.pm_3 = pm_3;
  qd.min_freq = min_freq;
  qd.max_freq = max_freq;
  qd.print_all_words = print_all_words;

  stream = gt4_sequence_stream_new (fname);
  if (!stream) {
    fprintf (stderr, "search_fasta: Cannot open %s\n", fname);
    return 1;
  }
  fasta_reader_init (&r, fname, qd.dict_inst->word_length, 0, GT4_SEQUENCE_STREAM_SEQUENCE_SOURCE_IMPLEMENTATION(stream), &stream->source_instance);
  result = fasta_reader_read_nwords (&r, 1000000000000ULL, NULL, NULL, NULL, NULL, process_word, (void *) &qd);
  fasta_reader_release (&r);
  az_object_shutdown (AZ_OBJECT (stream));
  return result;
}

static unsigned int
search_list (AZObject *obj, const char *fname, unsigned int n_mm, unsigned int pm_3, unsigned int min_freq, unsigned int max_freq, int print_all_words)
{
  GT4WordSListImplementation *s_impl;
  GT4WordSListInstance *s_inst;
  AZObject *s_obj = NULL;
  QueryData qd = {0};
  FILE *ifs;
  unsigned int code = 0;

  qd.dict_impl = (GT4WordDictImplementation *) az_object_get_interface (obj, GT4_TYPE_WORD_DICT, (void **) &qd.dict_inst);
  if (GT4_IS_INDEX_MAP (obj) && locations) {
    qd.index_impl = (GT4WordIndexImplementation *) az_object_get_interface (obj, GT4_TYPE_WORD_INDEX, (void **) &qd.index_inst);
  }
  qd.n_mm = n_mm;
  qd.pm_3 = pm_3;
  qd.min_freq = min_freq;
  qd.max_freq = max_freq;
  qd.print_all_words = print_all_words;

  ifs = fopen (fname, "r");
  if (!ifs) {
    fprintf (stderr, "search_list: Cannot open list %s\n", fname);
    return 1;
  }
  fread (&code, 4, 1, ifs);
  fclose (ifs);
  if (code == GT4_LIST_CODE) {
    s_obj = (AZObject *) gt4_word_map_new (fname, VERSION_MAJOR, 0, 0);
  } else if (code == GT4_INDEX_CODE) {
    s_obj = (AZObject *) gt4_index_map_new (fname, VERSION_MAJOR, 0);
  } else {
    fprintf (stderr, "search_list: Invalid file format\n");
    return 1;
  }
  s_impl = (GT4WordSListImplementation *) az_object_get_interface (s_obj, GT4_TYPE_WORD_SLIST, (void **) &s_inst);

  if (qd.dict_inst->word_length != s_inst->word_length) {
    az_object_shutdown (s_obj);
    return GT_INCOMPATIBLE_WORDLENGTH_ERROR;
  }
  
  gt4_word_slist_get_first_word (s_impl, s_inst);
  while (s_inst->idx < s_inst->num_words) {
    uint64_t word = s_inst->word;
    search_one_word (&qd, word);
    if (!gt4_word_slist_get_next_word (s_impl, s_inst)) break;
  }
  az_object_shutdown (s_obj);
  return 0;
}

unsigned int
search_lists_multi (AZObject *query, AZObject *lists[], unsigned int n_lists)
{
  GT4WordSListImplementation *impls[MAX_LISTS + 1];
  GT4WordSListInstance *insts[MAX_LISTS + 1];
  unsigned int i;

  impls[0] = (GT4WordSListImplementation *) az_object_get_interface (query, GT4_TYPE_WORD_SLIST, (void **) &insts[0]);
  gt4_word_slist_get_first_word (impls[0], insts[0]);
  for (i = 0; i < n_lists; i++) {
    impls[i + 1] = (GT4WordSListImplementation *) az_object_get_interface (lists[i], GT4_TYPE_WORD_SLIST, (void **) &insts[i + 1]);
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

void get_statistics (AZObject *obj)
{
  GT4WordSListInstance *inst;
  az_object_get_interface (obj, GT4_TYPE_WORD_SLIST, (void **) &inst);
  if (GT4_IS_WORD_MAP (obj)) {
    GT4WordMap *map = GT4_WORD_MAP(obj);
    fprintf (stdout, "List %s: built with glistmaker version %d.%d\n", map->filename, map->header.version_major, map->header.version_minor);
  } else if (GT4_IS_INDEX_MAP (obj)) {
    GT4IndexMap *imap = GT4_INDEX_MAP(obj);
    fprintf (stdout, "Index %s: built with glistmaker version %d.%d\n", imap->filename, imap->header->version_major, imap->header->version_minor);
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
    fprintf (stdout, "List %s: built with glistmaker version %d.%d\n", map->filename, map->header.version_major, map->header.version_minor);
  } else if (GT4_IS_INDEX_MAP (obj)) {
    GT4IndexMap *imap = GT4_INDEX_MAP(obj);
    fprintf (stdout, "Index %s: built with glistmaker version %d.%d\n", imap->filename, imap->header->version_major, imap->header->version_minor);
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
  fprintf (stderr, "    -stat, --stats            - print statistics of the list file and exit\n");
  fprintf (stderr, "    --median                  - print min/max/median/average and exit\n");
  fprintf (stderr, "    --distribution MAX        - print distribution up to MAX\n");
  fprintf (stderr, "    --gc                      - print average GC content of all words\n");
  fprintf (stderr, "    -q, --query               - single query word\n");
  fprintf (stderr, "    -f, --queryfile           - list of query words in a file\n");
  fprintf (stderr, "    -s, --seqfile             - FastA/FastQ file\n");
  fprintf (stderr, "    -l, --listfile            - list file made by glistmaker\n");
  fprintf (stderr, "    -mm, --mismatch NUMBER    - specify number of mismatches (0-16; default 0)\n");
  fprintf (stderr, "    -p, --perfectmatch NUMBER - specify number of 3' perfect matches (0-32; default 0)\n");
  fprintf (stderr, "    -min, --minfreq NUMBER    - minimum frequency of the printed words (default 0)\n");
  fprintf (stderr, "    -max, --maxfreq NUMBER    - maximum frequency of the printed words (default MAX_UINT)\n");
  fprintf (stderr, "    --files                   - Print indexed files\n");
  fprintf (stderr, "    --sequences               - Print indexed subsequences\n");
  fprintf (stderr, "    --bloom                   - use bloom filter to speed up lookups\n");
  fprintf (stderr, "    --all                     - in case of mismatches prints all found words\n");
  fprintf (stderr, "    --locations               - in case of index print all word locations\n");
  fprintf (stderr, "    --3p                      - if query is longer than word use 3' end\n");
  fprintf (stderr, "    --5p                      - if query is longer than word use 5' end\n");
  fprintf (stderr, "    -D                        - increase debug level\n");
  exit (exit_value);
}
