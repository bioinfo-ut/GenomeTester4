#define __GT4_SEQUENCE_STREAM_C__

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

#include <stdlib.h>
#include <string.h>

#include <libarikkei/arikkei-utils.h>

#include "sequence-stream.h"

static void sequence_stream_class_init (GT4SequenceStreamClass *klass);

/* AZObject implementation */
static void sequence_stream_shutdown (AZObject *obj);
/* GT4SequenceSource implementation */
static unsigned int sequence_stream_open (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst);
static int sequence_stream_read (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst);
static unsigned int sequence_stream_close (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst);

static unsigned int sequence_stream_type = 0;
static GT4SequenceStreamClass *sequence_stream_class = NULL;

unsigned int
gt4_sequence_stream_get_type (void)
{
  if (!sequence_stream_type) {
    sequence_stream_class = (GT4SequenceStreamClass *) az_register_type (&sequence_stream_type, AZ_TYPE_OBJECT, (const unsigned char *) "GT4SequenceStream",
      sizeof (GT4SequenceStreamClass), sizeof (GT4SequenceStream),
      (void (*) (AZClass *)) sequence_stream_class_init,
      NULL, NULL);
  }
  return sequence_stream_type;
}

static void
sequence_stream_class_init (GT4SequenceStreamClass *klass)
{
  ((AZClass *) klass)->flags |= AZ_CLASS_ZERO_MEMORY;
  az_class_set_num_interfaces ((AZClass *) klass, 1);
  az_class_declare_interface ((AZClass *) klass, 0, GT4_TYPE_SEQUENCE_SOURCE, ARIKKEI_OFFSET(GT4SequenceStreamClass,source_implementation), ARIKKEI_OFFSET(GT4SequenceStream,source_instance));
  /* AZObject implementation */
  ((AZObjectClass *) klass)->shutdown = sequence_stream_shutdown;
  /* GT4SequenceSource implementation */
  klass->source_implementation.open = sequence_stream_open;
  klass->source_implementation.read = sequence_stream_read;
  klass->source_implementation.close = sequence_stream_close;
}

static void
sequence_stream_shutdown (AZObject *obj)
{
  GT4SequenceStream *stream = (GT4SequenceStream *) obj;
  if (stream->source_instance.open && stream->close_stream) {
    gt4_sequence_source_close (GT4_SEQUENCE_STREAM_SEQUENCE_SOURCE_IMPLEMENTATION(obj), &stream->source_instance);
  }
}

static unsigned int
sequence_stream_open (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst)
{
  GT4SequenceStream *stream = GT4_SEQUENCE_STREAM_FROM_SEQUENCE_SOURCE_INSTANCE(inst);
  if (!stream->ifs) stream->ifs = fopen (stream->filename, "r");
  return stream->ifs != NULL;
}

static int
sequence_stream_read (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst)
{
  GT4SequenceStream *stream = GT4_SEQUENCE_STREAM_FROM_SEQUENCE_SOURCE_INSTANCE(inst);
  int result;
  unsigned char c;
  result = fread (&c, 1, 1, stream->ifs) ;
  if (!result) {
    if (feof (stream->ifs)) return 0;
    return -1;
  }
  return c;
}

static unsigned int
sequence_stream_close (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst)
{
  GT4SequenceStream *stream = GT4_SEQUENCE_STREAM_FROM_SEQUENCE_SOURCE_INSTANCE(inst);
  if (stream->ifs && stream->close_stream) {
    fclose (stream->ifs);
    stream->ifs = NULL;
  }
  return 1;
}

GT4SequenceStream *
gt4_sequence_stream_new (const char *filename)
{
  GT4SequenceStream *stream;
  arikkei_return_val_if_fail (filename != NULL, NULL);
  stream = (GT4SequenceStream *) az_object_new (GT4_TYPE_SEQUENCE_STREAM);
  stream->filename = strdup (filename);
  if (!strcmp (filename, "-")) {
    stream->ifs = stdin;
    stream->close_stream = 0;
  } else {
    stream->close_stream = 1;
  }
  return stream;
}

GT4SequenceStream *
gt4_sequence_stream_new_from_stream (FILE *ifs, unsigned int close_stream)
{
  GT4SequenceStream *stream;
  arikkei_return_val_if_fail (ifs != NULL, NULL);
  stream = (GT4SequenceStream *) az_object_new (GT4_TYPE_SEQUENCE_STREAM);
  stream->ifs = ifs;
  stream->close_stream = close_stream;
  return stream;
}
