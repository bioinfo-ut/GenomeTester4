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

unsigned int listmaker_queue_debug = 1;

#define BLOCK_SIZE 10000000000ULL

#include <string.h>
#include <sys/mman.h>

#include "sequence-file.h"
#include "sequence-stream.h"
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
listmaker_queue_init (GT4ListMakerQueueClass *klass, GT4ListMakerQueue *mq)
{
}

static void
listmaker_queue_finalize (GT4ListMakerQueueClass *klass, GT4ListMakerQueue *mq)
{
  free (mq->free_s_tables);
  free (mq->used_s_tables);
}

void
maker_queue_setup (GT4ListMakerQueue *mq, unsigned int n_threads, unsigned int wlen, unsigned int n_tmp_tables, unsigned int tmp_table_size)
{
  unsigned int i;
  az_instance_init (mq, GT4_TYPE_LISTMAKER_QUEUE);
  gt4_queue_setup (&mq->queue, n_threads);
  mq->wordlen = wlen;
  mq->free_s_tables = (wordtable **) malloc (n_tmp_tables * sizeof (wordtable *));
  mq->used_s_tables = (wordtable **) malloc (n_tmp_tables * sizeof (wordtable *));
  for (i = 0; i < n_tmp_tables; i++) {
    mq->free_s_tables[mq->n_free_s_tables++] = wordtable_new (wlen, tmp_table_size);
  }
}

void
maker_queue_release (GT4ListMakerQueue *mq)
{
  az_instance_finalize (mq, GT4_TYPE_LISTMAKER_QUEUE);
}

static void
maker_queue_add_source (GT4ListMakerQueue *mq, AZObject *src, const char *name)
{
  TaskFile *tf = task_file_new_from_source (AZ_OBJECT (src), name, 0);
  az_object_unref (AZ_OBJECT (src));
  tf->next = mq->files;
  mq->files = tf;
  mq->n_files_waiting += 1;

  TaskRead *tr = task_read_new (mq, src);
  gt4_queue_add_task (&mq->queue, &tr->task, 0);
}

void
maker_queue_add_file (GT4ListMakerQueue *mq, const char *filename, unsigned int stream)
{
  unsigned int len = (unsigned int) strlen (filename);
  if ((len > 3) && !strcmp (filename + len - 3, ".gz")) {
    /* gzipped file, add as ZStream */
    GT4SequenceZStream *zstream;
    if (listmaker_queue_debug) fprintf (stderr, "maker_queue_add_file: Compressed stream %s\n", filename);
    zstream = gt4_sequence_zstream_new (filename);
    if (!zstream) {
      fprintf (stderr, "Cannot open compressed stream %s\n", filename);
      return;
    }
    maker_queue_add_source (mq, AZ_OBJECT (zstream), filename);
  } else if (!strcmp (filename, "-")) {
    /* stdin, add as Stream */
    GT4SequenceStream *stream;
    if (listmaker_queue_debug) fprintf (stderr, "maker_queue_add_file: stdin\n");
    stream = gt4_sequence_stream_new_from_stream (stdin, 0);
    maker_queue_add_source (mq, AZ_OBJECT (stream), "stdin");
  } else if (stream) {
    /* standard file, add as Stream */
    GT4SequenceStream *stream;
    if (listmaker_queue_debug) fprintf (stderr, "maker_queue_add_file: Stream %s\n", filename);
    stream = gt4_sequence_stream_new (filename);
    maker_queue_add_source (mq, AZ_OBJECT (stream), filename);
  } else {
    unsigned int nblocks, i;
    GT4SequenceFile *seqf;
    GT4SequenceBlock *blocks[32];
    /* standard file, break into SequenceBlocks */
    seqf = gt4_sequence_file_new (filename, 1);
    gt4_sequence_file_map_sequence (seqf);
    nblocks = (unsigned int) (seqf->block.csize / BLOCK_SIZE) + 1;
    if (nblocks > 32) nblocks = 32;
    nblocks = gt4_sequence_block_split (&seqf->block, blocks, nblocks);
    if (listmaker_queue_debug) fprintf (stderr, "maker_queue_add_file: File %s as %u blocks\n", filename, nblocks);
    for (i = 0; i < nblocks; i++) {
      char c[32];
      snprintf (c, 32, "Block %u", i);
      if (listmaker_queue_debug) fprintf (stderr, "  Block %u from %llu to %llu\n", i, (unsigned long long) (blocks[i]->cdata - seqf->block.cdata), (unsigned long long) ((blocks[i]->cdata + blocks[i]->csize) - seqf->block.cdata));
      maker_queue_add_source (mq, AZ_OBJECT (blocks[i]), c);
    }
    az_object_unref (AZ_OBJECT (seqf));
  }
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

/* Tasks */

TaskRead *
task_read_new (GT4ListMakerQueue *mq, AZObject *source)
{
  TaskRead *tr;
  GT4SequenceSourceImplementation *impl;
  GT4SequenceSourceInstance *inst;
  tr = (TaskRead *) malloc (sizeof (TaskRead));
  memset (tr, 0, sizeof (TaskRead));
  tr->task.queue = &mq->queue;
  tr->task.type = TASK_READ;
  tr->task.priority = 10;
  tr->source = source;
  az_object_ref (AZ_OBJECT(tr->source));
  impl = (GT4SequenceSourceImplementation *) az_object_get_interface (AZ_OBJECT(tr->source), GT4_TYPE_SEQUENCE_SOURCE, (void **) &inst);
  gt4_sequence_source_open (impl, inst);
  fasta_reader_init (&tr->reader, mq->wordlen, 1, impl, inst);
  return tr;
}

void
task_read_delete (TaskRead *tr)
{
  fasta_reader_release (&tr->reader);
  az_object_unref (AZ_OBJECT (tr->source));
  free (tr);
}

TaskCollate *
task_collate_new (GT4ListMakerQueue *mq, unsigned int max_tables)
{
  TaskCollate *tc;
  tc = (TaskCollate *) malloc (sizeof (TaskCollate) + (max_tables - 1) * sizeof (wordtable *));
  memset (tc, 0, sizeof (TaskCollate));
  tc->task.queue = &mq->queue;
  tc->task.type = TASK_COLLATE;
  return tc;
}

void
task_collate_delete (TaskCollate *tc)
{
  free (tc);
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
task_file_new_from_source (AZObject *source, const char *name, unsigned int close)
{
  TaskFile *tf = (TaskFile *) malloc (sizeof (TaskFile));
  memset (tf, 0, sizeof (TaskFile));
  tf->seqfile = gt4_sequence_file_new (name, 1);
  tf->stream = NULL;
  tf->source = source;
  az_object_ref (AZ_OBJECT(tf->source));
  tf->close_source = close;
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
  if (tf->close_source) {
    GT4SequenceSourceImplementation *impl;
    GT4SequenceSourceInstance *inst;
    impl = (GT4SequenceSourceImplementation *) az_object_get_interface (AZ_OBJECT(tf->source), GT4_TYPE_SEQUENCE_SOURCE, (void **) &inst);
    gt4_sequence_source_close (impl, inst);
  }
  az_object_unref (AZ_OBJECT (tf->source));
  free (tf);
}

/* Frontend to mmap and GT4FastaReader */

unsigned int
task_file_read_nwords (TaskFile *tf, unsigned long long maxwords, unsigned int wordsize,
  /* Called as soon as the full sequence name is known */
  int (*start_sequence) (GT4FastaReader *, void *),
  /* Called when the full sequence has been parsed */
  int (*end_sequence) (GT4FastaReader *, void *),
  int (*read_character) (GT4FastaReader *, unsigned int character, void *),
  int (*read_nucleotide) (GT4FastaReader *, unsigned int nucleotide, void *),
  int (*read_word) (GT4FastaReader *, unsigned long long word, void *),
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
