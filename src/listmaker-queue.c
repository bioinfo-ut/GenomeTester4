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
  unsigned int i;
  free (mq->free_s_tables);
  free (mq->used_s_tables);
  for (i = 0; i < mq->n_tmp_files; i++) free (mq->tmp_files[i]);
  for (i = 0; i < mq->n_final_files; i++) free (mq->final_files[i]);
}

void
maker_queue_setup (GT4ListMakerQueue *mq, unsigned int n_threads, unsigned int wlen, unsigned int n_tmp_tables, unsigned int tmp_table_size)
{
  unsigned int i;
  az_instance_init (mq, GT4_TYPE_LISTMAKER_QUEUE);
  gt4_queue_setup (&mq->queue, n_threads);
  mq->wordlen = wlen;
  mq->free_s_tables = (GT4WordTable **) malloc (n_tmp_tables * sizeof (GT4WordTable *));
  mq->used_s_tables = (GT4WordTable **) malloc (n_tmp_tables * sizeof (GT4WordTable *));
  for (i = 0; i < n_tmp_tables; i++) {
    mq->free_s_tables[mq->n_free_s_tables++] = gt4_word_table_new (wlen, tmp_table_size);
  }
}

void
maker_queue_release (GT4ListMakerQueue *mq)
{
  unsigned int i;
  for (i = 0; i < mq->n_free_s_tables; i++) gt4_word_table_delete (mq->free_s_tables[i]);
  for (i = 0; i < mq->n_used_s_tables; i++) gt4_word_table_delete (mq->used_s_tables[i]);
  az_instance_finalize (mq, GT4_TYPE_LISTMAKER_QUEUE);
}

static void
maker_queue_add_source (GT4ListMakerQueue *mq, AZObject *source, const char *name)
{
  az_object_ref (AZ_OBJECT(source));
  TaskRead *tr = task_read_new (&mq->queue, source, mq->n_blocks++, mq->wordlen);
  gt4_queue_add_task (&mq->queue, &tr->task, 0);
  mq->n_files_waiting += 1;
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

/* Tasks */

void
gt4_task_read_setup (TaskRead *tr, GT4Queue *queue, AZObject *source, unsigned int idx, unsigned int wordlen)
{
  GT4SequenceSourceImplementation *impl;
  GT4SequenceSourceInstance *inst;
  memset (tr, 0, sizeof (TaskRead));
  tr->task.queue = queue;
  tr->task.type = TASK_READ;
  tr->task.priority = 10;
  tr->source = source;
  az_object_ref (AZ_OBJECT(source));
  tr->idx = idx;
  impl = (GT4SequenceSourceImplementation *) az_object_get_interface (AZ_OBJECT(source), GT4_TYPE_SEQUENCE_SOURCE, (void **) &inst);
  gt4_sequence_source_open (impl, inst);
  fasta_reader_init (&tr->reader, wordlen, 1, impl, inst);
}

void
gt4_task_read_release (TaskRead *tr)
{
  fasta_reader_release (&tr->reader);
  az_object_unref (AZ_OBJECT (tr->source));
}

TaskRead *
task_read_new (GT4Queue *queue, AZObject *source, unsigned int idx, unsigned int wordlen)
{
  TaskRead *tr;
  tr = (TaskRead *) malloc (sizeof (TaskRead));
  gt4_task_read_setup (tr, queue, source, idx, wordlen);
  return tr;
}

void
task_read_delete (TaskRead *tr)
{
  gt4_task_read_release (tr);
  free (tr);
}

TaskCollateTables *
task_collate_tables_new (GT4ListMakerQueue *mq, unsigned int max_tables)
{
  TaskCollateTables *tc;
  tc = (TaskCollateTables *) malloc (sizeof (TaskCollateTables) + (max_tables - 1) * sizeof (GT4WordTable *));
  memset (tc, 0, sizeof (TaskCollateTables));
  tc->task.queue = &mq->queue;
  tc->task.type = TASK_COLLATE_TABLES;
  tc->task.priority = 9;
  return tc;
}

void
task_collate_tables_delete (TaskCollateTables *tc)
{
  free (tc);
}

TaskCollateFiles *
task_collate_files_new (GT4ListMakerQueue *mq, unsigned int max_files)
{
  TaskCollateFiles *tc;
  tc = (TaskCollateFiles *) malloc (sizeof (TaskCollateFiles) + (max_files - 1) * sizeof (char *));
  memset (tc, 0, sizeof (TaskCollateFiles));
  tc->task.queue = &mq->queue;
  tc->task.type = TASK_COLLATE_FILES;
  tc->task.priority = 8;
  return tc;
}

void
task_collate_files_delete (TaskCollateFiles *tc)
{
  unsigned int i;
  for (i = 0; i < tc->n_files; i++) {
    free (tc->files[i]);
  }
  free (tc);
}

