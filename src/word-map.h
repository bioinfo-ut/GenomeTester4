#ifndef __GT4_WORD_MAP_H__
#define __GT4_WORD_MAP_H__

/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 * 
 * Copyright (C) 2014-2016 University of Tartu
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

/*
 * GT4WordMap is the most basic list container
 */

typedef struct _GT4WordMap GT4WordMap;
typedef struct _GT4WordMapClass GT4WordMapClass;

#define GT4_TYPE_WORD_MAP (gt4_word_map_get_type ())
#define GT4_WORD_MAP(o) (AZ_CHECK_INSTANCE_CAST ((o), GT4_TYPE_WORD_MAP, GT4WordMap))
#define GT4_IS_WORD_MAP(o) (AZ_CHECK_INSTANCE_TYPE ((o), GT4_TYPE_WORD_MAP))

#define GT4_WORD_MAP_FROM_SARRAY_INST(i) (GT4WordMap *) AZ_BASE_ADDRESS(GT4WordMap,sarray_inst,i)
#define GT4_WORD_MAP_SARRAY_IMPL(o) &((GT4WordMapClass *) ((AZObject *) (o))->klass)->sarray_impl
#define GT4_WORD_MAP_FROM_DICT_INST(i) (GT4WordMap *) AZ_BASE_ADDRESS(GT4WordMap,dict_inst,i)
#define GT4_WORD_MAP_DICT_IMPL(o) &((GT4WordMapClass *) ((AZObject *) (o))->klass)->dict_impl

#ifndef __GT4_WORD_MAP_C__
/* Class pointer, valid after get_type is called once */
extern GT4WordMapClass *gt4_word_map_class;
/* Defaults to 0, can be increased to print debug information */
extern unsigned int debug_wordmap;
#endif

#include <az/object.h>

#include "bloom.h"
#include "utils.h"
#include "word-array-sorted.h"
#include "word-dict.h"
#include "word-list.h"

#define WORDMAP_ELEMENT_SIZE (sizeof (unsigned long long) + sizeof (unsigned int))

struct _GT4WordMap {
  AZObject object;
  char *filename;
  const unsigned char *file_map;
  unsigned long long file_size;
  const unsigned char *wordlist;
  GT4ListHeader header;
  /* GT4WordSArray instance */
  GT4WordSArrayInstance sarray_inst;
  /* GT4WordDict instance */
  GT4WordDictInstance dict_inst;
  /* Background mapper */
  GT4Scout scout;
  /* Bloom filter */
  GT4Bloom *bloom;
  unsigned long long reject;
  unsigned long long pass;
};

struct _GT4WordMapClass {
  AZObjectClass object_class;
  /* GT4WordSArray implementation */
  GT4WordSArrayImplementation sarray_impl;
  /* GT4WordDict implementation */
  GT4WordDictImplementation dict_impl;
};

unsigned int gt4_word_map_get_type (void);

inline extern uint64_t
gt4_word_map_get_word (const GT4WordMap *wmap, uint64_t idx)
{
  return *((uint64_t *) (wmap->wordlist + 12 * idx));
}

inline extern uint32_t
gt4_word_map_get_count (const GT4WordMap *wmap, uint64_t idx)
{
  return *((uint32_t *) (wmap->wordlist + 12 * idx + 8));
}

inline extern uint64_t *
gt4_word_map_get_word_ptr (const GT4WordMap *wmap, uint64_t idx)
{
  return (uint64_t *) (wmap->wordlist + 12 * idx);
}

#define WORDMAP_WORD gt4_word_map_get_word
#define WORDMAP_FREQ gt4_word_map_get_count
/*
#define WORDMAP_WORD(w,i) (*((unsigned long long *) ((w)->wordlist + 12 * (i))))
#define WORDMAP_FREQ(w,i) (*((unsigned int *) ((w)->wordlist + 12 * (i) + 8)))
*/
/* Creates new GT4WordMap by memory-mapping file, returns NULL if error */
/* If "scout" is true, a new thread is created that sequentially prefetces the map into virtual memory */
GT4WordMap *gt4_word_map_new (const char *listfilename, unsigned int major_version, unsigned int scout, unsigned int create_bloom);
/* Releases word map and frees the structure */
void gt4_word_map_delete (GT4WordMap *map);

unsigned int word_map_search_query (GT4WordMap *wmap, unsigned long long query, unsigned int n_mm, unsigned int pm_3, int printall, unsigned int equalmmonly, unsigned int dosubtraction, GT4WordMap *querymap);

unsigned int gt4_word_map_lookup_canonical (GT4WordMap *wmap, unsigned long long query);
unsigned int gt4_word_map_lookup (GT4WordMap *wmap, unsigned long long query);

#endif /* WORDMAP_H_ */
