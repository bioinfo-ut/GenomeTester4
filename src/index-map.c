#define __GT4_INDEX_MAP_C__

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
#include "utils.h"

#include "index-map.h"

unsigned int debug_index_map = 0;

unsigned int GT4_INDEX_CODE = 'G' << 24 | 'T' << 16 | '4' << 8 | 'I';

static void index_map_class_init (GT4IndexMapClass *klass);
/* AZObject implementation */
static void index_map_shutdown (AZObject *object);

/* GT4WordSList implementation */
static unsigned int index_map_get_first_word (GT4WordSListImplementation *impl, GT4WordSListInstance *inst);
static unsigned int index_map_get_next_word (GT4WordSListImplementation *impl, GT4WordSListInstance *inst);
/* GT4WordDict implementation */
static unsigned int index_map_lookup (GT4WordDictImplementation *impl, GT4WordDictInstance *inst, unsigned long long word);
/* GT4FileArrayImplementation */
static unsigned int index_map_get_file (GT4FileArrayImplementation *impl, GT4FileArrayInstance *inst, unsigned int idx);

static unsigned int index_map_type = 0;
GT4IndexMapClass *gt4_index_map_class = NULL;

unsigned int
gt4_index_map_get_type (void)
{
  if (!index_map_type) {
    gt4_index_map_class = (GT4IndexMapClass *) az_register_type (&index_map_type, AZ_TYPE_OBJECT, (const unsigned char *) "GT4IndexMap",
      sizeof (GT4IndexMapClass), sizeof (GT4IndexMap),
      (void (*) (AZClass *)) index_map_class_init,
      NULL, NULL);
  }
  return index_map_type;
}

static void
index_map_class_init (GT4IndexMapClass *klass)
{
  klass->object_class.shutdown = index_map_shutdown;
  az_class_set_num_interfaces ((AZClass *) klass, 3);
  az_class_declare_interface ((AZClass *) klass, 0, GT4_TYPE_WORD_SARRAY, ARIKKEI_OFFSET (GT4IndexMapClass, sarray_implementation), ARIKKEI_OFFSET (GT4IndexMap, sarray_instance));
  az_class_declare_interface ((AZClass *) klass, 1, GT4_TYPE_WORD_DICT, ARIKKEI_OFFSET (GT4IndexMapClass, dict_impl), ARIKKEI_OFFSET (GT4IndexMap, dict_inst));
  az_class_declare_interface ((AZClass *) klass, 2, GT4_TYPE_FILE_ARRAY, ARIKKEI_OFFSET (GT4IndexMapClass, file_array_impl), ARIKKEI_OFFSET (GT4IndexMap, file_array_inst));
  /* GT4WordSList implementation */
  klass->sarray_implementation.slist_impl.get_first_word = index_map_get_first_word;
  klass->sarray_implementation.slist_impl.get_next_word = index_map_get_next_word;
  /* GT4WordDict implementation */
  klass->dict_impl.lookup = index_map_lookup;
  /* GT4FileArrayImplementation */
  klass->file_array_impl.get_file = index_map_get_file;
}

static void
index_map_shutdown (AZObject *object)
{
  GT4IndexMap *imap = (GT4IndexMap *) object;
  if (imap->filename) {
    free (imap->filename);
    imap->filename = NULL;
  }
  if (imap->file_map) {
    munmap ((void *) imap->file_map, imap->file_size);
    imap->file_map = NULL;
    imap->file_size = 0;
  }
  imap->header = NULL;
  imap->files = NULL;
  imap->kmers = NULL;
  imap->locations = NULL;
}

static unsigned long long
imap_get_word (GT4IndexMap *imap, unsigned long long idx)
{
  return *((unsigned long long *) (imap->kmers + idx * 16));
}

static unsigned int
imap_get_count (GT4IndexMap *imap, unsigned long long idx)
{
  unsigned long long next_idx = idx + 1;
  unsigned long long loc = *((unsigned long long *) (imap->kmers + idx * 16 + 8));
  if (next_idx == imap->header->num_words) {
    return (unsigned int) (imap->header->num_locations - loc);
  } else {
    unsigned long long next_loc = *((unsigned long long *) (imap->kmers + next_idx * 16 + 8));
    return (unsigned int) (next_loc - loc);
  }
}

static unsigned int
index_map_get_first_word (GT4WordSListImplementation *impl, GT4WordSListInstance *inst)
{
  GT4IndexMap *imap = GT4_INDEX_MAP_FROM_SARRAY_INSTANCE(inst);
  inst->word = imap_get_word (imap, 0);
  inst->count = imap_get_count (imap, 0);
  return 1;
}

static unsigned int
index_map_get_next_word (GT4WordSListImplementation *impl, GT4WordSListInstance *inst)
{
  GT4IndexMap *imap = GT4_INDEX_MAP_FROM_SARRAY_INSTANCE(inst);
  inst->word = imap_get_word (imap, inst->idx);
  inst->count = imap_get_count (imap, inst->idx);
  __builtin_prefetch (imap->kmers + (inst->idx + 4) * 16);
  return 1;
}

static unsigned int
index_map_lookup (GT4WordDictImplementation *impl, GT4WordDictInstance *inst, unsigned long long word)
{
  GT4IndexMap *imap = GT4_INDEX_MAP_FROM_DICT_INST(inst);
  unsigned long long current, low, high, mid;
  low = 0;
  high = imap->header->num_words - 1;
  mid = (low + high) / 2;
      
  while (low <= high) {
    current = imap_get_word (imap, mid);
    if (current < word) {
      low = mid + 1;
    } else if (current > word) {
      if (mid == 0) break;
      high = mid - 1;
    } else {
      imap->dict_inst.value = imap_get_count (imap, mid);
      return 1;
    }
    mid = (low + high) / 2;
  }
  return 0;
}

static unsigned int
index_map_get_file (GT4FileArrayImplementation *impl, GT4FileArrayInstance *inst, unsigned int idx)
{
  GT4IndexMap *imap = GT4_INDEX_MAP_FROM_FILE_ARRAY_INST(inst);
  const unsigned char *s = imap->files + 16;
  unsigned int i;
  for (i = 0; i < idx; i++) {
    unsigned short len;
    memcpy (&len, s, 2);
    s += 2;
    s += len;
  }
  inst->file_name = s + 2;
  inst->file_size = 0;
  inst->n_sequences = 0;
  return 1;
}

GT4IndexMap * 
gt4_index_map_new (const char *listfilename, unsigned int major_version, unsigned int scout)
{
  const unsigned char *cdata;
  unsigned long long csize;
  GT4IndexMap *imap;

  cdata = gt4_mmap (listfilename, &csize);
  if (!cdata) {
    fprintf (stderr, "gt4_index_map_new: could not mmap file %s\n", listfilename);
    return NULL;
  }

  imap = (GT4IndexMap *) az_object_new (GT4_TYPE_INDEX_MAP);
  if (!imap) {
    fprintf (stderr, "gt4_index_map_new: could not allocate map\n");
    return NULL;
  }

  imap->filename = strdup (listfilename);
  imap->file_map = cdata;
  imap->file_size = csize;
  imap->header = (GT4IndexHeader *) cdata;
  if (imap->header->code != GT4_INDEX_CODE) {
    fprintf (stderr, "gt4_index_map_new: invalid file tag (%x, should be %x)\n", imap->header->code, GT4_INDEX_CODE);
    gt4_index_map_delete (imap);
    return NULL;
  }
  if (imap->header->version_major != major_version) {
    fprintf (stderr, "gt4_index_map_new: incompatible major version %u (required %u)\n", imap->header->version_major, major_version);
    gt4_index_map_delete (imap);
    return NULL;
  }
  imap->files = cdata + imap->header->files_start;
  imap->kmers = cdata + imap->header->kmers_start;
  imap->locations = cdata + imap->header->locations_start;
  if (scout) {
    scout_mmap ((const unsigned char *) cdata, csize);
  }

  /* Set up sorted array interface */
  imap->sarray_instance.slist_inst.num_words = imap->header->num_words;
  imap->sarray_instance.slist_inst.sum_counts = imap->header->num_locations;
  imap->sarray_instance.slist_inst.word_length = imap->header->word_length;
  if (imap->sarray_instance.slist_inst.num_words > 0) {
    imap->sarray_instance.slist_inst.word = imap_get_word (imap, 0);
    imap->sarray_instance.slist_inst.count = imap_get_count (imap, 0);
  }

  /* Set up file array interface */
  memcpy (&imap->file_array_inst.num_files, imap->files + 12, 4);

  return imap;
}

void
gt4_index_map_delete (GT4IndexMap *imap)
{
  az_object_shutdown (AZ_OBJECT (imap));
}

unsigned int 
gt4_index_map_lookup_canonical (GT4IndexMap *imap, unsigned long long query)
{
  if (gt4_word_dict_lookup (GT4_INDEX_MAP_DICT_IMPL(imap), &imap->dict_inst, query, 0)) {
    return imap->dict_inst.value;
  }
  return 0;
}

unsigned int 
gt4_index_map_lookup (GT4IndexMap *imap, unsigned long long query)
{
  if (gt4_word_dict_lookup (GT4_INDEX_MAP_DICT_IMPL(imap), &imap->dict_inst, query, 0)) {
    return imap->dict_inst.value;
  }
  return 0;
}
