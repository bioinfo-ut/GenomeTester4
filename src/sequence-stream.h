#ifndef __GT4_SEQUENCE_STREAM_H__
#define __GT4_SEQUENCE_STREAM_H__

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
 * A FILE stream of FastA/FastQ data
 */

typedef struct _GT4SequenceStream GT4SequenceStream;
typedef struct _GT4SequenceStreamClass GT4SequenceStreamClass;

#define GT4_TYPE_SEQUENCE_STREAM (gt4_sequence_stream_get_type ())
#define GT4_SEQUENCE_STREAM(o) (AZ_CHECK_INSTANCE_CAST ((o), GT4_TYPE_SEQUENCE_STREAM, GT4SequenceStream))
#define GT4_IS_SEQUENCE_STREAM(o) (AZ_CHECK_INSTANCE_TYPE ((o), GT4_TYPE_SEQUENCE_STREAM))

#define GT4_SEQUENCE_STREAM_FROM_SEQUENCE_SOURCE_INSTANCE(i) (GT4SequenceStream *) AZ_BASE_ADDRESS(GT4SequenceStream,source_instance,i)
#define GT4_SEQUENCE_STREAM_SEQUENCE_SOURCE_IMPLEMENTATION(o) &((GT4SequenceStreamClass *) ((AZObject *) (o))->klass)->source_implementation

#include <stdio.h>

#include <az/object.h>

#include "sequence-source.h"

struct _GT4SequenceStream {
  AZObject object;
  char *filename;
  FILE *ifs;
  unsigned int close_stream : 1;
  /* GT4SequenceSource instance */
  GT4SequenceSourceInstance source_instance;
};

struct _GT4SequenceStreamClass {
  AZObjectClass object_class;
  /* GT4SequenceSource implementation */
  GT4SequenceSourceImplementation source_implementation;
};

unsigned int gt4_sequence_stream_get_type (void);

/* If filename is "-" stdin is used */
GT4SequenceStream *gt4_sequence_stream_new (const char *filename);
GT4SequenceStream *gt4_sequence_stream_new_from_stream (FILE *ifs, unsigned int close_stream);

#endif
