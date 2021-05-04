#ifndef __GT4_WORD_LIST_STREAM_H__
#define __GT4_WORD_LIST_STREAM_H__

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
 * A stream-based word list container
 */

typedef struct _GT4WordListStream GT4WordListStream;
typedef struct _GT4WordListStreamClass GT4WordListStreamClass;

#define GT4_TYPE_WORD_LIST_STREAM (gt4_word_list_stream_get_type ())
#define GT4_WORD_LIST_STREAM(o) (AZ_CHECK_INSTANCE_CAST ((o), GT4_TYPE_WORD_LIST_STREAM, GT4WordListStream))
#define GT4_IS_WORD_LIST_STREAM(o) (AZ_CHECK_INSTANCE_TYPE ((o), GT4_TYPE_WORD_LIST_STREAM))

#define GT4_WORD_LIST_STREAM_FROM_SLIST_INSTANCE(i) (GT4WordListStream *) AZ_BASE_ADDRESS(GT4WordListStream,slist_instance,i)
#define GT4_WORD_LIST_STREAM_SLIST_IMPLEMENTATION(o) &((GT4WordListStreamClass *) ((AZObject *) (o))->klass)->slist_implementation

#define GT4_WORD_LIST_STREAM_BUF_SIZE (12 * 256)

#include <stdio.h>

#include <az/object.h>

#include "word-map.h"
#include "word-list-sorted.h"

struct _GT4WordListStream {
  AZObject object;
  char *filename;
  int ifile;
  unsigned char b[GT4_WORD_LIST_STREAM_BUF_SIZE];
  unsigned int bp, bsize;
  GT4ListHeader header;
  /* GT4WordSArray instance */
  GT4WordSListInstance slist_instance;
};

struct _GT4WordListStreamClass {
  AZObjectClass object_class;
  /* GT4WordSArray implementation */
  GT4WordSListImplementation slist_implementation;
};

unsigned int gt4_word_list_stream_get_type (void);

GT4WordListStream *gt4_word_list_stream_new (const char *filename, unsigned int major_version);
void gt4_word_list_stream_delete (GT4WordListStream *map);

#endif
