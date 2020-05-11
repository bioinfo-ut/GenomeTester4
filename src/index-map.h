#ifndef __GT4_INDEX_MAP_H__
#define __GT4_INDEX_MAP_H__

/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 * 
 * Copyright (C) 2014-2018 University of Tartu
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
 * GT4IndexMap is mmap-ed container for sorted index
 */

#ifndef __GT4_INDEX_MAP_C__
/* List tag "GT4C" encoded to big-endian 32-bit integer */
extern unsigned int GT4_INDEX_CODE;
#endif

typedef struct _GT4IndexMap GT4IndexMap;
typedef struct _GT4IndexMapClass GT4IndexMapClass;
typedef struct _GT4IndexHeader GT4IndexHeader;

#define GT4_TYPE_INDEX_MAP (gt4_index_map_get_type ())
#define GT4_INDEX_MAP(o) (AZ_CHECK_INSTANCE_CAST ((o), GT4_TYPE_INDEX_MAP, GT4IndexMap))
#define GT4_IS_INDEX_MAP(o) (AZ_CHECK_INSTANCE_TYPE ((o), GT4_TYPE_INDEX_MAP))

#define GT4_INDEX_MAP_FROM_SARRAY_INSTANCE(i) (GT4IndexMap *) AZ_BASE_ADDRESS(GT4IndexMap,sarray_inst,i)
#define GT4_INDEX_MAP_SLIST_IMPLEMENTATION(o) &((GT4IndexMapClass *) ((AZObject *) (o))->klass)->sarray_impl.slist_impl
#define GT4_INDEX_MAP_SARRAY_IMPLEMENTATION(o) &((GT4IndexMapClass *) ((AZObject *) (o))->klass)->sarray_impl
#define GT4_INDEX_MAP_FROM_DICT_INST(i) (GT4IndexMap *) AZ_BASE_ADDRESS(GT4IndexMap,dict_inst,i)
#define GT4_INDEX_MAP_DICT_IMPL(o) &((GT4IndexMapClass *) ((AZObject *) (o))->klass)->dict_impl
#define GT4_INDEX_MAP_FROM_INDEX_INST(i) (GT4IndexMap *) AZ_BASE_ADDRESS(GT4IndexMap,index_inst,i)
#define GT4_INDEX_MAP_INDEX_IMPL(o) &((GT4IndexMapClass *) ((AZObject *) (o))->klass)->index_impl
#define GT4_INDEX_MAP_FROM_FILE_ARRAY_INST(i) (GT4IndexMap *) AZ_BASE_ADDRESS(GT4IndexMap,file_array_inst,i)
#define GT4_INDEX_MAP_FILE_ARRAY_IMPL(o) &((GT4IndexMapClass *) ((AZObject *) (o))->klass)->file_array_impl

#ifndef __GT4_INDEX_MAP_C__
/* Class pointer, valid after get_type is called once */
extern GT4IndexMapClass *gt4_index_map_class;
/* Defaults to 0, can be increased to print debug information */
extern unsigned int debug_index_map;
#endif

#include <az/object.h>

#include "file-array.h"
#include "word-array-sorted.h"
#include "word-dict.h"
#include "word-index.h"

struct _GT4IndexHeader {
  unsigned int code;
  unsigned int version_major;
  unsigned int version_minor;
  unsigned int word_length;
  unsigned long long num_words;
  unsigned long long num_locations;
  unsigned int n_file_bits;
  unsigned int n_subseq_bits;
  unsigned int n_pos_bits;
  unsigned int filler0;
  unsigned long long files_start;
  unsigned long long kmers_start;
  unsigned long long locations_start;
};

struct _GT4IndexMap {
  AZObject object;
  char *filename;
  const unsigned char *file_map;
  unsigned long long file_size;
  const GT4IndexHeader *header;
  const unsigned char *files;
  const unsigned char *kmers;
  const unsigned char *locations;
  const unsigned char *file_ptr;
  /* Temporarily mapped source */
  unsigned int src_idx;
  const unsigned char *src_map;
  unsigned long long src_size;
  GT4WordSArrayInstance sarray_inst;
  GT4WordDictInstance dict_inst;
  GT4WordIndexInstance index_inst;
  GT4FileArrayInstance file_array_inst;
};

struct _GT4IndexMapClass {
  AZObjectClass object_class;
  GT4WordSArrayImplementation sarray_impl;
  GT4WordDictImplementation dict_impl;
  GT4WordIndexImplementation index_impl;
  GT4FileArrayImplementation file_array_impl;
};

unsigned int gt4_index_map_get_type (void);

/* Creates new GT4IndexMap by memory-mapping file, returns NULL if error */
/* If "scout" is true, a new thread is created that sequentially prefetces the map into virtual memory */
GT4IndexMap *gt4_index_map_new (const char *listfilename, unsigned int major_version, unsigned int scout);
/* Releases word map and frees the structure */
void gt4_index_map_delete (GT4IndexMap *map);

unsigned int gt4_index_map_lookup_canonical (GT4IndexMap *imap, unsigned long long query);
unsigned int gt4_index_map_lookup (GT4IndexMap *imap, unsigned long long query);

unsigned int gt4_index_map_get_location (GT4IndexMap *imap, unsigned int kmer_idx, unsigned int loc_idx, unsigned int *file_idx, unsigned int *subseq_idx, unsigned long long *pos, unsigned int *dir);

unsigned int gt4_index_map_get_sequence_name (GT4IndexMap *imap, unsigned char b[], unsigned int b_len, unsigned int file_idx, unsigned int seq_idx);

#endif
