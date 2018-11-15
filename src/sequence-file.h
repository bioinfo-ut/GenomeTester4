#ifndef __GT4_SEQUENCE_FILE_H__
#define __GT4_SEQUENCE_FILE_H__

/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 *
 * Copyright (C) 2014-2017 University of Tartu
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
 * Subclass of GT4SequenceBlock for memory mapped whole files
 *
 */

typedef struct _GT4SequenceFile GT4SequenceFile;
typedef struct _GT4SequenceFileClass GT4SequenceFileClass;

#define GT4_TYPE_SEQUENCE_FILE (gt4_sequence_file_get_type ())
#define GT4_SEQUENCE_FILE(o) (AZ_CHECK_INSTANCE_CAST ((o), GT4_TYPE_SEQUENCE_FILE, GT4SequenceFile))
#define GT4_IS_SEQUENCE_FILE(o) (AZ_CHECK_INSTANCE_TYPE ((o), GT4_TYPE_SEQUENCE_FILE))

#define GT4_SEQUENCE_FILE_FROM_SEQUENCE_SOURCE_INSTANCE(i) (GT4SequenceFile *) AZ_BASE_ADDRESS(GT4SequenceBlock,source_inst,i)
#define GT4_SEQUENCE_FILE_SEQUENCE_SOURCE_IMPLEMENTATION(o) &((GT4SequenceBlockClass *) ((AZObject *) (o))->klass)->source_impl

#include <pthread.h>

#include "sequence-block.h"

struct _GT4SequenceFile {
  GT4SequenceBlock block;
  char *path;
  /* File length */
  unsigned long long size;
  /* Flags */
  /* Whether to use explicit locking */
  unsigned int lock : 1;
  /* For concurrent access */
  pthread_mutex_t mutex;
};

struct _GT4SequenceFileClass {
  GT4SequenceBlockClass block_class;
};

unsigned int gt4_sequence_file_get_type (void);

/* Creates new object, does not memory map or parse it */
GT4SequenceFile *gt4_sequence_file_new (const char *path, unsigned int lock);

/* These implement locking at the top of AZobject */
void gt4_sequence_file_ref (GT4SequenceFile *seqfile);
void gt4_sequence_file_unref (GT4SequenceFile *seqfile);

void gt4_sequence_file_lock (GT4SequenceFile *seqfile);
void gt4_sequence_file_unlock (GT4SequenceFile *seqfile);

/* Memory maps sequence if not already mapped */
void gt4_sequence_file_map_sequence (GT4SequenceFile *seqfile);

/* Returns new subsequence index */
unsigned int gt4_sequence_file_add_subsequence (GT4SequenceFile *seqfile, unsigned long long name_pos, unsigned int name_len);

#endif
