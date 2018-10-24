#define __GT4_WORD_MAP_C__

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>

#include "sequence.h"
#include "common.h"
#include "utils.h"
#include "listmaker-queue.h"
#include "word-table.h"

#include "word-map.h"

unsigned int debug_wordmap = 0;

static void word_map_class_init (GT4WordMapClass *klass);
/* AZObject implementation */
static void word_map_shutdown (AZObject *object);

/* GT4WordSArray implementation */
unsigned int word_map_get_first_word (GT4WordSArrayImplementation *impl, GT4WordSArrayInstance *inst);
unsigned int word_map_get_next_word (GT4WordSArrayImplementation *impl, GT4WordSArrayInstance *inst);

static unsigned int word_map_type = 0;
GT4WordMapClass *gt4_word_map_class = NULL;

unsigned int
gt4_word_map_get_type (void)
{
  if (!word_map_type) {
    gt4_word_map_class = (GT4WordMapClass *) az_register_type (&word_map_type, AZ_TYPE_OBJECT, (const unsigned char *) "GT4WordMap",
      sizeof (GT4WordMapClass), sizeof (GT4WordMap),
      (void (*) (AZClass *)) word_map_class_init,
      NULL, NULL);
  }
  return word_map_type;
}

static void
word_map_class_init (GT4WordMapClass *klass)
{
  klass->object_class.shutdown = word_map_shutdown;
  az_class_set_num_interfaces ((AZClass *) klass, 1);
  az_class_declare_interface ((AZClass *) klass, 0, GT4_TYPE_WORD_SARRAY, ARIKKEI_OFFSET (GT4WordMapClass, sarray_implementation), ARIKKEI_OFFSET (GT4WordMap, sarray_instance));
  /* GT4WordSArray implementation */
  klass->sarray_implementation.get_first_word = word_map_get_first_word;
  klass->sarray_implementation.get_next_word = word_map_get_next_word;
}

static void
word_map_shutdown (AZObject *object)
{
  GT4WordMap *wmap = (GT4WordMap *) object;
  if (wmap->filename) {
    free (wmap->filename);
    wmap->filename = NULL;
  }
  if (wmap->file_map) {
    munmap ((void *) wmap->file_map, wmap->file_size);
    wmap->file_map = NULL;
    wmap->file_size = 0;
  }
  wmap->header = NULL;
  wmap->wordlist = NULL;
  wmap->user_data = NULL;
}

unsigned int
word_map_get_first_word (GT4WordSArrayImplementation *impl, GT4WordSArrayInstance *inst)
{
  GT4WordMap *wmap = GT4_WORD_MAP_FROM_SARRAY_INSTANCE(inst);
  inst->word = WORDMAP_WORD(wmap,0);
  inst->count = WORDMAP_FREQ(wmap,0);
  return 1;
}

unsigned int
word_map_get_next_word (GT4WordSArrayImplementation *impl, GT4WordSArrayInstance *inst)
{
  GT4WordMap *wmap = GT4_WORD_MAP_FROM_SARRAY_INSTANCE(inst);
  inst->word = WORDMAP_WORD(wmap,inst->idx);
  inst->count = WORDMAP_FREQ(wmap,inst->idx);
  __builtin_prefetch (&WORDMAP_WORD(wmap,inst->idx + 4), 0, 0);
  return 1;
}

GT4WordMap * 
gt4_word_map_new (const char *listfilename, unsigned int major_version, unsigned int scout)
{
  const unsigned char *cdata;
  unsigned long long csize, start;
  GT4WordMap *wmap;

  cdata = gt4_mmap (listfilename, &csize);
  if (!cdata) {
    fprintf (stderr, "gt4_word_map_new: could not mmap file %s\n", listfilename);
    return NULL;
  }

  wmap = (GT4WordMap *) az_object_new (GT4_TYPE_WORD_MAP);
  if (!wmap) {
    fprintf (stderr, "gt4_word_map_new: could not allocate map\n");
    return NULL;
  }

  wmap->filename = strdup (listfilename);
  wmap->file_map = cdata;
  wmap->file_size = csize;
  wmap->header = (GT4ListHeader *) cdata;
  if (wmap->header->code != GT4_LIST_CODE) {
    fprintf (stderr, "gt4_word_map_new: invalid file tag (%x, should be %x)\n", wmap->header->code, GT4_LIST_CODE);
    gt4_word_map_delete (wmap);
    return NULL;
  }
  if (wmap->header->version_major != major_version) {
    fprintf (stderr, "gt4_word_map_new: incompatible major version %u (required %u)\n", wmap->header->version_major, major_version);
    gt4_word_map_delete (wmap);
    return NULL;
  }
  if ((wmap->header->version_major == 4) && (wmap->header->version_minor == 0)) {
    start = sizeof (GT4ListHeader);
  } else {
    start = wmap->header->list_start;
  }
  wmap->wordlist = cdata + start;
  if (csize < start + wmap->header->nwords * 12) {
    fprintf (stderr, "gt4_word_map_new: file size too small (%llu, should be at least %llu)\n", csize, start + wmap->header->nwords * 12);
    gt4_word_map_delete (wmap);
    return NULL;
  }
  if (scout) {
    scout_mmap ((const unsigned char *) cdata, csize);
  }

  /* Set up sorted array interface */
  wmap->sarray_instance.num_words = wmap->header->nwords;
  if (wmap->sarray_instance.num_words > 0) {
    wmap->sarray_instance.word = WORDMAP_WORD(wmap,0);
    wmap->sarray_instance.count = WORDMAP_FREQ(wmap,0);
  }
  wmap->sarray_instance.word_length = wmap->header->wordlength;

  return wmap;
}

void
gt4_word_map_delete (GT4WordMap *wmap)
{
  az_object_shutdown (AZ_OBJECT (wmap));
}

unsigned int 
word_map_search_query (GT4WordMap *map, unsigned long long query, parameters *p, int printall, unsigned int equalmmonly, unsigned int dosubtraction, GT4WordMap *querymap)
{
  static GT4WordTable mm_table = {0};
  unsigned long long i, nwords = 0L;
  unsigned int count = 0L, currentcount = 0L, querycount = 0L;

  /* if no mismatches */
  if (!p->nmm) {
    return gt4_word_map_lookup (map, query);
  }

  mm_table.wordlength = p->wordlength;

  /* find and set table size */
  if (!mm_table.n_words) {
    nwords = generate_mismatches (NULL, query, p->wordlength, 0, p->nmm, p->pm3, 0, 1, 0);
    gt4_word_table_ensure_size (&mm_table, nwords);
    if (debug_wordmap > 1) {
      fprintf (stderr, "MM Table size %llu, num mismatches %llu\n", mm_table.n_words, nwords);
    }
  }
  generate_mismatches (&mm_table, query, p->wordlength, 0, p->nmm, p->pm3, 0, 0, equalmmonly);
  if (debug_wordmap > 1) {
    fprintf (stderr, "MM Table size %llu\n", mm_table.n_words);
  }

  for (i = 0; i < mm_table.n_words; i++) {
    if (dosubtraction) {
      querycount = gt4_word_map_lookup (querymap, mm_table.words[i]);
      currentcount = gt4_word_map_lookup (map, mm_table.words[i]);
      if (currentcount > querycount) {
        if (debug_wordmap > 1) {
          fprintf (stderr, "%llu %llu %llu querycount %u currentcount %u\n", query, i, mm_table.words[i], querycount, currentcount);
        }
        mm_table.n_words = 0;
        return ~0L;
      }
      count += (currentcount - querycount);
    } else {
      currentcount = gt4_word_map_lookup (map, mm_table.words[i]);
      count += currentcount;
      if (printall && currentcount > 0) {
        fprintf (stdout, "%s\t%u\n", word_to_string (mm_table.words[i], mm_table.wordlength), currentcount);
      }
    }
  }
  
  mm_table.n_words = 0;
  return count;
}

unsigned int 
gt4_word_map_lookup_canonical (GT4WordMap *map, unsigned long long query)
{
  unsigned long long word, low, high, mid;
  low = 0;
  high = map->header->nwords - 1;
  mid = (low + high) / 2;
      
  while (low <= high) {
    word = WORDMAP_WORD (map, mid);
    if (word < query) {
      low = mid + 1;
    } else if (word > query) {
      if (mid == 0) break;
      high = mid - 1;
    } else {
      return WORDMAP_FREQ (map, mid);
    }
    mid = (low + high) / 2;
  }
  return 0;
}

unsigned int 
gt4_word_map_lookup (GT4WordMap *map, unsigned long long query)
{
  unsigned long long rev = get_reverse_complement (query, map->header->wordlength);
  if (rev < query) query = rev;
  return gt4_word_map_lookup_canonical (map, query);
}




