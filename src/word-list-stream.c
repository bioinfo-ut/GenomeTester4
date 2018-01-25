#define __GT4_WORD_LIST_STREAM_C__

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

#include <malloc.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#include "word-list-stream.h"

static void word_list_stream_class_init (GT4WordListStreamClass *klass);
/* AZObject implementation */
static void word_list_stream_shutdown (AZObject *object);

/* GT4WordSArray implementation */
unsigned int word_list_stream_get_first_word (GT4WordSArrayImplementation *impl, GT4WordSArrayInstance *inst);
unsigned int word_list_stream_get_next_word (GT4WordSArrayImplementation *impl, GT4WordSArrayInstance *inst);

static unsigned int word_list_stream_type = 0;
GT4WordListStreamClass *gt4_word_list_stream_class = NULL;

unsigned int
gt4_word_list_stream_get_type (void)
{
  if (!word_list_stream_type) {
    gt4_word_list_stream_class = (GT4WordListStreamClass *) az_register_type (&word_list_stream_type, AZ_TYPE_OBJECT, (const unsigned char *) "GT4WordListStream",
      sizeof (GT4WordListStreamClass), sizeof (GT4WordListStream),
      (void (*) (AZClass *)) word_list_stream_class_init,
      NULL, NULL);
  }
  return word_list_stream_type;
}

static void
word_list_stream_class_init (GT4WordListStreamClass *klass)
{
  klass->object_class.shutdown = word_list_stream_shutdown;
  az_class_set_num_interfaces ((AZClass *) klass, 1);
  az_class_declare_interface ((AZClass *) klass, 0, GT4_TYPE_WORD_SARRAY, ARIKKEI_OFFSET (GT4WordListStreamClass, sarray_implementation), ARIKKEI_OFFSET (GT4WordListStream, sarray_instance));
  /* GT4WordSArray implementation */
  klass->sarray_implementation.get_first_word = word_list_stream_get_first_word;
  klass->sarray_implementation.get_next_word = word_list_stream_get_next_word;
}

static void
word_list_stream_shutdown (AZObject *object)
{
  GT4WordListStream *stream = (GT4WordListStream *) object;
  if (stream->filename) {
    free (stream->filename);
    stream->filename = NULL;
  }
  if (stream->ifs) {
    fclose (stream->ifs);
    stream->ifs = NULL;
  }
}

unsigned int
word_list_stream_get_first_word (GT4WordSArrayImplementation *impl, GT4WordSArrayInstance *inst)
{
  GT4WordListStream *stream = GT4_WORD_LIST_STREAM_FROM_SARRAY_INSTANCE(inst);
  fseeko (stream->ifs, sizeof (GT4ListHeader), SEEK_SET);
  fread (&stream->sarray_instance.word, 8, 1, stream->ifs);
  fread (&stream->sarray_instance.count, 4, 1, stream->ifs);
  return 1;
}

unsigned int
word_list_stream_get_next_word (GT4WordSArrayImplementation *impl, GT4WordSArrayInstance *inst)
{
  GT4WordListStream *stream = GT4_WORD_LIST_STREAM_FROM_SARRAY_INSTANCE(inst);
  fread (&stream->sarray_instance.word, 8, 1, stream->ifs);
  fread (&stream->sarray_instance.count, 4, 1, stream->ifs);
  return 1;
}

GT4WordListStream * 
gt4_word_list_stream_new (const char *filename, unsigned int major_version)
{
  GT4WordListStream *stream;
  FILE *ifs;
  
  arikkei_return_val_if_fail (filename != NULL, NULL);

  ifs = fopen (filename, "r");
  if (!ifs) {
    fprintf (stderr, "gt4_word_list_stream_new: could not open file %s\n", filename);
    return NULL;
  }

  stream = (GT4WordListStream *) az_object_new (GT4_TYPE_WORD_LIST_STREAM);
  if (!stream) {
    fprintf (stderr, "gt4_word_list_stream_new: could not allocate object\n");
    return NULL;
  }

  stream->filename = strdup (filename);
  stream->ifs = ifs;

  if (fread (&stream->header, 1, sizeof (GT4ListHeader), ifs) != sizeof (GT4ListHeader)) {
    fprintf (stderr, "gt4_word_list_stream_new: could not read list header\n");
    gt4_word_list_stream_delete (stream);
    return NULL;
  }
  if (stream->header.code != GT4_LIST_CODE) {
    fprintf (stderr, "gt4_word_list_stream_new: invalid file tag (%x, should be %x)\n", stream->header.code, GT4_LIST_CODE);
    gt4_word_list_stream_delete (stream);
    return NULL;
  }
  if (stream->header.version_major > major_version) {
    fprintf (stderr, "gt4_word_list_stream_new: incompatible major version %u (required %u)\n", stream->header.version_major, major_version);
    gt4_word_list_stream_delete (stream);
    return NULL;
  }

  /* Set up sorted array interface */
  stream->sarray_instance.num_words = stream->header.nwords;
  if (stream->sarray_instance.num_words > 0) {
    fread (&stream->sarray_instance.word, 8, 1, stream->ifs);
    fread (&stream->sarray_instance.count, 4, 1, stream->ifs);
  }
  stream->sarray_instance.word_length = stream->header.wordlength;

  return stream;
}

void
gt4_word_list_stream_delete (GT4WordListStream *stream)
{
  az_object_shutdown (AZ_OBJECT (stream));
}

