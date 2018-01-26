#define __GT4_LISTMAKER_QUEUE_C__

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

unsigned int listmaker_queue_debug = 0;

#include <string.h>
#include <sys/mman.h>

#include "sequence-zstream.h"
#include "utils.h"
#include "listmaker-queue.h"

static void listmaker_queue_class_init (GT4ListMakerQueueClass *klass);
static void listmaker_queue_init (GT4ListMakerQueueClass *klass, GT4ListMakerQueue *queue);
static void listmaker_queue_finalize (GT4ListMakerQueueClass *klass, GT4ListMakerQueue *queue);

static unsigned int listmaker_queue_type = 0;

unsigned int
gt4_listmaker_queue_get_type (void)
{
  if (!listmaker_queue_type) {
    az_register_type (&listmaker_queue_type, GT4_TYPE_QUEUE, (const unsigned char *) "GT4ListMakerQueue",
      sizeof (GT4ListMakerQueueClass), sizeof (GT4ListMakerQueue),
      (void (*) (AZClass *)) listmaker_queue_class_init,
      (void (*) (AZImplementation *, void *)) listmaker_queue_init,
      (void (*) (AZImplementation *, void *)) listmaker_queue_finalize);
  }
  return listmaker_queue_type;
}

static void
listmaker_queue_class_init (GT4ListMakerQueueClass *klass)
{
}

static void
listmaker_queue_init (GT4ListMakerQueueClass *klass, GT4ListMakerQueue *queue)
{
}

static void
listmaker_queue_finalize (GT4ListMakerQueueClass *klass, GT4ListMakerQueue *queue)
{
}

void
maker_queue_setup (GT4ListMakerQueue *mq, unsigned int nthreads)
{
  az_instance_init (mq, GT4_TYPE_LISTMAKER_QUEUE);
  gt4_queue_setup (&mq->queue, nthreads);
}

void
maker_queue_release (GT4ListMakerQueue *mq)
{
  az_instance_finalize (mq, GT4_TYPE_LISTMAKER_QUEUE);
}

void
maker_queue_add_file (GT4ListMakerQueue *mq, const char *filename)
{
        TaskFile *tf;
#if 0
        FILE *ifs;
        ifs = fopen (filename, "r");
        tf = task_file_new_from_stream (ifs, filename, 1);
#else
  unsigned int len = strlen (filename);
  if ((len > 3) && !strcmp (filename + len - 3, ".gz")) {
    if (listmaker_queue_debug) fprintf (stderr, "Opening compressed stream %s\n", filename);
    GT4SequenceZStream *zstream = gt4_sequence_zstream_new (filename);
    TaskFile *tf = task_file_new_from_source (AZ_OBJECT (zstream), filename, 0);
    tf->next = mq->files;
    //tf->idx = mq->nfiles++;
    mq->files = tf;
  } else if (!strcmp (filename, "-")) {
    tf = task_file_new_from_stream (stdin, filename, 0);
    tf->next = mq->files;
    //tf->idx = mq->nfiles++;
    mq->files = tf;
  } else {
    unsigned int nseqs, i;
    GT4SequenceBlock *seqs[32];
    GT4SequenceFile *seqf = gt4_sequence_file_new (filename, 1);
    gt4_sequence_file_map_sequence (seqf);
    nseqs = (unsigned int) (seqf->block.csize / 10000000000ULL) + 1;
    if (nseqs > 32) nseqs = 32;
    nseqs = gt4_sequence_block_split (&seqf->block, seqs, nseqs);
    for (i = 0; i < nseqs; i++) {
      TaskFile *tf = task_file_new_from_source (AZ_OBJECT (seqs[i]), "block", 0);
      if (listmaker_queue_debug) fprintf (stderr, "%s:%u from %llu to %llu\n", filename, i, (unsigned long long) seqs[i]->cdata - (unsigned long long) seqf->block.cdata, (unsigned long long) seqs[i]->cdata - (unsigned long long) seqf->block.cdata + seqs[i]->csize);
      tf->next = mq->files;
      //tf->idx = mq->nfiles++;
      mq->files = tf;
    }
    az_object_unref (AZ_OBJECT (seqf));
  }


        //tf = task_file_new (filename, 0);
#endif
        //tf->next = mq->files;
        //mq->files = tf;
}

wordtable *
queue_get_smallest_table (GT4ListMakerQueue *queue)
{
	unsigned int i, min;
	unsigned long long minslots;
	wordtable *t;
	min = 0;
	minslots = queue->available[0]->nwordslots;
	for (i = 1; i < queue->navailable; i++) {
		if (queue->available[i]->nwordslots < minslots) {
			min = i;
			minslots = queue->available[i]->nwordslots;
		}
	}
	t = queue->available[min];
	queue->available[min] = queue->available[queue->navailable - 1];
	queue->navailable -= 1;
	return t;
}

wordtable *
queue_get_largest_table (GT4ListMakerQueue *queue)
{
	unsigned int i, max;
	unsigned long long maxslots;
	wordtable *t;
	max = 0;
	maxslots = queue->available[0]->nwordslots;
	for (i = 1; i < queue->navailable; i++) {
		if (queue->available[i]->nwordslots > maxslots) {
			max = i;
			maxslots = queue->available[i]->nwordslots;
		}
	}
	t = queue->available[max];
	queue->available[max] = queue->available[queue->navailable - 1];
	queue->navailable -= 1;
	return t;
}

wordtable *
queue_get_sorted (GT4ListMakerQueue *queue)
{
	queue->nsorted -= 1;
	return queue->sorted[queue->nsorted];
}

wordtable *
queue_get_smallest_sorted (GT4ListMakerQueue *queue)
{
	unsigned int i, min;
	unsigned long long minwords;
	wordtable *t;
	min = 0;
	minwords = queue->sorted[0]->nwords;
	for (i = 1; i < queue->nsorted; i++) {
		if (queue->sorted[i]->nwords < minwords) {
			min = i;
			minwords = queue->sorted[i]->nwords;
		}
	}
	t = queue->sorted[min];
	queue->sorted[min] = queue->sorted[queue->nsorted - 1];
	queue->nsorted -= 1;
	return t;
}

wordtable *
queue_get_mostavailable_sorted (GT4ListMakerQueue *queue)
{
	unsigned int i, max;
	unsigned long long maxavail;
	wordtable *t;
	max = 0;
	maxavail = queue->sorted[0]->nwordslots - queue->sorted[0]->nwords;
	for (i = 1; i < queue->nsorted; i++) {
		if ((queue->sorted[i]->nwordslots - queue->sorted[0]->nwords) > maxavail) {
			max = i;
			maxavail = queue->sorted[i]->nwordslots - queue->sorted[i]->nwords;
		}
	}
	t = queue->sorted[max];
	queue->sorted[max] = queue->sorted[queue->nsorted - 1];
	queue->nsorted -= 1;
	return t;
}

static unsigned int finish = 0;

static void *
scout_map (void *arg)
{
  MapData *map = (MapData *) arg;
  unsigned long long i, val = 0;
  for (i = 0; i < map->csize; i += 1000) {
    val += map->cdata[i];
    if (finish) break;
  }
  free (map);
  return (void *) val;
}

void
scout_mmap (const unsigned char *cdata, unsigned long long csize)
{
  pthread_t thread;
  MapData *map = (MapData *) malloc (sizeof (MapData));
  map->cdata = cdata;
  map->csize = csize;
  pthread_create (&thread, NULL, scout_map, map);
}

void
delete_scouts () {
  finish = 1;
}

/* File parsing tasks */

TaskFile *
task_file_new (const char *filename, unsigned int scout)
{
  TaskFile *tf = (TaskFile *) malloc (sizeof (TaskFile));
  memset (tf, 0, sizeof (TaskFile));
  tf->seqfile = gt4_sequence_file_new (filename, 1);
  tf->source = (AZObject *) tf->seqfile;
  az_object_ref (AZ_OBJECT(tf->source));
  tf->scout = scout;
  return tf;
}

TaskFile *
task_file_new_from_stream (FILE *ifs, const char *filename, unsigned int close_on_delete)
{
  TaskFile *tf = (TaskFile *) malloc (sizeof (TaskFile));
  memset (tf, 0, sizeof (TaskFile));
  tf->seqfile = gt4_sequence_file_new (filename, 1);
  tf->stream = gt4_sequence_stream_new_from_stream (ifs, close_on_delete);
  tf->source = (AZObject *) tf->stream;
  az_object_ref (AZ_OBJECT(tf->source));
  tf->close_on_delete = close_on_delete;
  return tf;
}

TaskFile *
task_file_new_from_source (AZObject *source, const char *filename, unsigned int close_on_delete)
{
  TaskFile *tf = (TaskFile *) malloc (sizeof (TaskFile));
  memset (tf, 0, sizeof (TaskFile));
  tf->seqfile = gt4_sequence_file_new (filename, 1);
  tf->stream = NULL;
  tf->source = source;
  az_object_ref (AZ_OBJECT(tf->source));
  tf->close_on_delete = close_on_delete;
  return tf;
}

void
task_file_delete (TaskFile *tf)
{
  if (tf->has_reader) {
    fasta_reader_release (&tf->reader);
  }
  gt4_sequence_file_unref (tf->seqfile);
  if (tf->stream) {
    az_object_shutdown (AZ_OBJECT (tf->stream));
  }
  az_object_unref (AZ_OBJECT (tf->source));
  free (tf);
}

/* Frontend to mmap and FastaReader */

unsigned int
task_file_read_nwords (TaskFile *tf, unsigned long long maxwords, unsigned int wordsize,
  /* Called as soon as the full sequence name is known */
  int (*start_sequence) (FastaReader *, void *),
  /* Called when the full sequence has been parsed */
  int (*end_sequence) (FastaReader *, void *),
  int (*read_character) (FastaReader *, unsigned int character, void *),
  int (*read_nucleotide) (FastaReader *, unsigned int nucleotide, void *),
  int (*read_word) (FastaReader *, unsigned long long word, void *),
  void *data)
{
  if (!tf->has_reader) {
    GT4SequenceSourceImplementation *impl;
    GT4SequenceSourceInstance *inst;
    impl = (GT4SequenceSourceImplementation *) az_object_get_interface (AZ_OBJECT(tf->source), GT4_TYPE_SEQUENCE_SOURCE, (void **) &inst);
    if (!inst->open) {
      if (!gt4_sequence_source_open (impl, inst)) {
        fprintf (stderr, "Cannot open sequence source of %s\n", tf->seqfile->path);
        return 1;
      }
      /* fixme: move scouting inside SequenceFile */
      if (tf->scout && (tf->source == (AZObject *) tf->seqfile)) scout_mmap (tf->seqfile->block.cdata, tf->seqfile->block.csize);
    }
    fasta_reader_init (&tf->reader, wordsize, 1, impl, inst);
    tf->has_reader = 1;
  }
  return fasta_reader_read_nwords (&tf->reader, maxwords, start_sequence, end_sequence, read_character, read_nucleotide, read_word, data);
}
