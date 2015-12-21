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

#include <string.h>
#include <sys/mman.h>

#include "utils.h"
#include "queue.h"

/* Initialize/destroy mutexes and conds */
/* Return 0 on success */

unsigned int
queue_init (Queue *queue)
{
  pthread_mutex_init (&queue->mutex, NULL);
  pthread_cond_init (&queue->cond, NULL);
  return 0;
}

unsigned int
queue_finalize (Queue *queue)
{
  pthread_cond_destroy (&queue->cond);
  pthread_mutex_destroy (&queue->mutex);
  return 0;
}

unsigned int
queue_lock (Queue *queue)
{
  pthread_mutex_lock (&queue->mutex);
  return 0;
}

unsigned int
queue_unlock (Queue *queue)
{
  pthread_mutex_unlock (&queue->mutex);
  return 0;
}

unsigned int
queue_wait (Queue *queue)
{
  pthread_cond_wait (&queue->cond, &queue->mutex);
  return 0;
}

unsigned int
queue_broadcast (Queue *queue)
{
  pthread_cond_broadcast (&queue->cond);
  return 0;
}

void
maker_queue_setup (MakerQueue *mq)
{
	memset (mq, 0, sizeof (MakerQueue));
	pthread_mutex_init (&mq->queue.mutex, NULL);
	pthread_cond_init (&mq->queue.cond, NULL);
}

void
maker_queue_release (MakerQueue *mq)
{
	pthread_cond_destroy (&mq->queue.cond);
	pthread_mutex_destroy (&mq->queue.mutex);
}

void
maker_queue_add_file (MakerQueue *mq, const char *filename)
{
        TaskFile *task;
        task = (TaskFile *) malloc (sizeof (TaskFile));
        memset (task, 0, sizeof (TaskFile));
        task->filename = filename;
        task->next = mq->files;
        mq->files = task;
}

wordtable *
queue_get_smallest_table (MakerQueue *queue)
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
queue_get_largest_table (MakerQueue *queue)
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
queue_get_sorted (MakerQueue *queue)
{
	queue->nsorted -= 1;
	return queue->sorted[queue->nsorted];
}

wordtable *
queue_get_smallest_sorted (MakerQueue *queue)
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
queue_get_mostavailable_sorted (MakerQueue *queue)
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
task_file_new (const char *filename)
{
  TaskFile *tf = (TaskFile *) malloc (sizeof (TaskFile));
  memset (tf, 0, sizeof (TaskFile));
  tf->filename = filename;
  return tf;
}

void
task_file_delete (TaskFile *tf)
{
  if (tf->has_reader) {
    fasta_reader_release (&tf->reader);
  }
  if (tf->cdata) {
    munmap ((void *) tf->cdata, tf->csize);
  }
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
  if (!tf->cdata) {
    size_t csize;
    tf->cdata = (const unsigned char *) mmap_by_filename ((const char *) tf->filename, &csize);
    if (!tf->cdata) {
      fprintf (stderr, "Cannot mmap %s\n", tf->filename);
      return 1;
    }
    tf->csize = csize;
    scout_mmap (tf->cdata, tf->csize);
    fasta_reader_init_from_data (&tf->reader, wordsize, 1, tf->cdata, tf->csize);
    tf->has_reader = 1;
  }
  fasta_reader_read_nwords (&tf->reader, maxwords, start_sequence, end_sequence, read_character, read_nucleotide, read_word, data);
  return 0;
}

