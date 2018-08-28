#define __GT4_SEQUENCE_ZSTREAM_C__

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
#include <zlib.h>

#include <libarikkei/arikkei-utils.h>

#include "sequence-zstream.h"

static const unsigned int zstream_debug = 0;

#define GT4_SEQUENCE_STREAM_Z_CHUNK_SIZE 16 * 16384
#define windowBits 15
#define ENABLE_ZLIB_GZIP 32

#define GZIP_BITS (windowBits + ENABLE_ZLIB_GZIP)

static void sequence_zstream_class_init (GT4SequenceZStreamClass *klass);
static void sequence_zstream_init (GT4SequenceZStreamClass *klass, GT4SequenceZStream *stream);
static void sequence_zstream_finalize (GT4SequenceZStreamClass *klass, GT4SequenceZStream *stream);

/* AZObject implementation */
static void sequence_zstream_shutdown (AZObject *obj);
/* GT4SequenceSource implementation */
static unsigned int sequence_zstream_open (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst);
static int sequence_zstream_read (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst);
static unsigned int sequence_zstream_close (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst);

static unsigned int sequence_zstream_type = 0;
static GT4SequenceZStreamClass *sequence_zstream_class = NULL;

unsigned int
gt4_sequence_zstream_get_type (void)
{
  if (!sequence_zstream_type) {
    sequence_zstream_class = (GT4SequenceZStreamClass *) az_register_type (&sequence_zstream_type, AZ_TYPE_OBJECT, (const unsigned char *) "GT4SequenceZStream",
      sizeof (GT4SequenceZStreamClass), sizeof (GT4SequenceZStream),
      (void (*) (AZClass *)) sequence_zstream_class_init,
      (void (*) (AZImplementation *, void *)) sequence_zstream_init,
      (void (*) (AZImplementation *, void *)) sequence_zstream_finalize);
  }
  return sequence_zstream_type;
}

static void
sequence_zstream_class_init (GT4SequenceZStreamClass *klass)
{
  ((AZClass *) klass)->flags |= AZ_CLASS_ZERO_MEMORY;
  az_class_set_num_interfaces ((AZClass *) klass, 1);
  az_class_declare_interface ((AZClass *) klass, 0, GT4_TYPE_SEQUENCE_SOURCE, ARIKKEI_OFFSET(GT4SequenceZStreamClass,source_implementation), ARIKKEI_OFFSET(GT4SequenceZStream,source_instance));
  /* AZObject implementation */
  ((AZObjectClass *) klass)->shutdown = sequence_zstream_shutdown;
  /* GT4SequenceSource implementation */
  klass->source_implementation.open = sequence_zstream_open;
  klass->source_implementation.read = sequence_zstream_read;
  klass->source_implementation.close = sequence_zstream_close;
}

static void
sequence_zstream_init (GT4SequenceZStreamClass *klass, GT4SequenceZStream *stream)
{
  stream->z_in = (unsigned char *) malloc (GT4_SEQUENCE_STREAM_Z_CHUNK_SIZE);
  stream->z_out = (unsigned char *) malloc (GT4_SEQUENCE_STREAM_Z_CHUNK_SIZE);
}

static void
sequence_zstream_finalize (GT4SequenceZStreamClass *klass, GT4SequenceZStream *stream)
{
  free (stream->z_in);
  free (stream->z_out);
}

static void
sequence_zstream_shutdown (AZObject *obj)
{
  GT4SequenceZStream *stream = (GT4SequenceZStream *) obj;
  if (stream->source_instance.open && stream->close_ifs) {
    gt4_sequence_source_close (GT4_SEQUENCE_ZSTREAM_SEQUENCE_SOURCE_IMPLEMENTATION(obj), &stream->source_instance);
  }
}

static unsigned int
sequence_zstream_open (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst)
{
  GT4SequenceZStream *stream = GT4_SEQUENCE_ZSTREAM_FROM_SEQUENCE_SOURCE_INSTANCE(inst);
  int z_status;

  if (!stream->ifs) stream->ifs = fopen (stream->filename, "r");
  if (!stream->ifs) return 0;

  stream->z_strm.next_in = NULL;
  stream->z_strm.avail_in = 0;
  z_status = inflateInit2 (&stream->z_strm, GZIP_BITS);
  if (z_status < 0) return 0;
  
  return 1;
}

static int
sequence_zstream_read (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst)
{
  GT4SequenceZStream *stream = GT4_SEQUENCE_ZSTREAM_FROM_SEQUENCE_SOURCE_INSTANCE(inst);
  int result;

  while (stream->z_out_pos >= stream->z_out_len) {
    /* Outbuffer is exhausted */
    if (!stream->ifs_eof && !stream->has_input) {
      /* Read more input */
      result = fread (stream->z_in, 1, GT4_SEQUENCE_STREAM_Z_CHUNK_SIZE, stream->ifs);
      if (feof (stream->ifs)) {
        if (zstream_debug) fprintf (stderr, "sequence_zstream_read: end of stream\n");
        stream->ifs_eof = 1;
      }
      if (!result) {
        if (!stream->ifs_eof) {
          if (zstream_debug > 1) fprintf (stderr, "sequence_zstream_read: error reading stream\n");
          return -1;
        }
        return 0;
      } else {
        if (zstream_debug > 2) fprintf (stderr, "sequence_zstream_read: read %u bytes from ifs\n", result);
        stream->total_ifs += result;
        stream->z_strm.next_in = stream->z_in;
        stream->z_strm.avail_in = result;
        stream->has_input = 1;
      }
    }
    /* Try to inflate more */
    if (zstream_debug > 2) fprintf (stderr, "sequence_zstream_read: avail_in %u\n", stream->z_strm.avail_in);
    if (zstream_debug > 2) fprintf (stderr, "sequence_zstream_read: out_pos %u out_len %u\n", stream->z_out_pos, stream->z_out_len);
    stream->z_strm.avail_out = GT4_SEQUENCE_STREAM_Z_CHUNK_SIZE;
    stream->z_strm.next_out = stream->z_out;
    result = inflate (&stream->z_strm, Z_NO_FLUSH);
    if (result < 0) {
      if (result != Z_BUF_ERROR) {
        /* Error in decoding */
        if (zstream_debug) fprintf (stderr, "sequence_zstream_read: inflate error %d\n", result);
        return -1;
      } else {
        /* Probably need more data in buffer */
        if (zstream_debug > 1) fprintf (stderr, "sequence_zstream_read: buffer error (avail_out %u)\n", stream->z_strm.avail_out);
      }
    }
    if (result == Z_STREAM_END) {
      /* Given compressed stream ended */
      if (stream->has_input) {
        if (result == Z_STREAM_END) {
          /* Try to start new stream */
          if (zstream_debug > 1) fprintf (stderr, "sequence_zstream_read: not eof, trying to start new stream\n");
          result = inflateEnd (&stream->z_strm);
          if (zstream_debug > 1) fprintf (stderr, "sequence_zstream_read: inflateEnd result %d\n", result);
          if (result < 0) return -1;
          result = inflateInit2 (&stream->z_strm, GZIP_BITS);
          if (zstream_debug > 1) fprintf (stderr, "sequence_zstream_read: inflateInit result %d\n", result);
          if (result < 0) return -1;
          result = inflate (&stream->z_strm, Z_NO_FLUSH);
          if (zstream_debug > 1) fprintf (stderr, "sequence_zstream_read: inflate result %d\n", result);
          if ((result < 0) && (result != Z_BUF_ERROR)) return -1;
        }
      }
    }
    if (zstream_debug > 1) fprintf (stderr, "sequence_zstream_read: uncompressed %u bytes\n", GT4_SEQUENCE_STREAM_Z_CHUNK_SIZE - stream->z_strm.avail_out);
    stream->z_out_len = GT4_SEQUENCE_STREAM_Z_CHUNK_SIZE - stream->z_strm.avail_out;
    stream->z_out_pos = 0;
    if (stream->z_strm.avail_out < GT4_SEQUENCE_STREAM_Z_CHUNK_SIZE) {
      /* Wrote some data */
      break;
    } else {
      if (zstream_debug > 1) fprintf (stderr, "sequence_zstream_read: need more input\n");
      stream->has_input = 0;
      if (stream->ifs_eof) return 0;
    }
  }
  return stream->z_out[stream->z_out_pos++];
}

static unsigned int
sequence_zstream_close (GT4SequenceSourceImplementation *impl, GT4SequenceSourceInstance *inst)
{
  GT4SequenceZStream *stream = GT4_SEQUENCE_ZSTREAM_FROM_SEQUENCE_SOURCE_INSTANCE(inst);

  inflateEnd (&stream->z_strm);

  if (stream->ifs && stream->close_ifs) {
    fclose (stream->ifs);
    stream->ifs = NULL;
  }
  return 1;
}

GT4SequenceZStream *
gt4_sequence_zstream_new (const char *filename)
{
  GT4SequenceZStream *stream;
  arikkei_return_val_if_fail (filename != NULL, NULL);
  stream = (GT4SequenceZStream *) az_object_new (GT4_TYPE_SEQUENCE_ZSTREAM);
  stream->filename = strdup (filename);
  if (!strcmp (filename, "-")) {
    stream->ifs = stdin;
    stream->close_ifs = 0;
  } else {
    stream->close_ifs = 1;
  }
  return stream;
}

GT4SequenceZStream *
gt4_sequence_zstream_new_from_stream (FILE *ifs, unsigned int close_stream)
{
  GT4SequenceZStream *stream;
  arikkei_return_val_if_fail (ifs != NULL, NULL);
  stream = (GT4SequenceZStream *) az_object_new (GT4_TYPE_SEQUENCE_ZSTREAM);
  stream->ifs = ifs;
  stream->close_ifs = close_stream;
  return stream;
}
