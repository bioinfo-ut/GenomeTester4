#ifndef QUEUE_H_
#define QUEUE_H_

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

#include <stdio.h>
#include <pthread.h>

#include "wordtable.h"
#include "fasta.h"

/*
 * Queue stub
 *
 */

typedef union _QValue QValue;
union _QValue {
  unsigned long long lval;
  double dval;
};

#define NUM_ID_TOKENS 32

typedef struct _Queue Queue;
struct _Queue {
  /* Threads */
  unsigned int nthreads_total;
  unsigned int nthreads_running;
  pthread_t *threads;
  /* Single mutex and cond for queue management */
  pthread_mutex_t mutex;
  pthread_cond_t cond;
  /* Debug */
  QValue tokens[NUM_ID_TOKENS];
};

/* Initialize/destroy mutexes and conds */
/* Return 0 on success */
/* Thread 0 is main, queue will be created with nthreads - 1 entries */
unsigned int queue_init (Queue *queue, unsigned int nthreads);
unsigned int queue_create_threads (Queue *queue, void (*process) (Queue *, unsigned int, void *), void *data);
unsigned int queue_finalize (Queue *queue);
unsigned int queue_lock (Queue *queue);
unsigned int queue_unlock (Queue *queue);
/* Should be called with mutex locked, returns mutex locked */
unsigned int queue_wait (Queue *queue);
unsigned int queue_broadcast (Queue *queue);

#define MAX_TABLES 256

#define TASK_READ 0
#define TASK_SORT 1
#define TASK_MERGE 2
#define NUM_TASK_TYPES 3

/*
 * Possible tasks are:
 *   - read segment from FastA file and build table
 *   - sort table
 *   - merge two tables
 */

/* Task for parsing FastA file */
/* These form linked list */

typedef struct _TaskFile TaskFile;
typedef struct _MakerQueue MakerQueue;

struct _MakerQueue {
        Queue queue;

        /* Parameters */
        unsigned int wordlen;
        unsigned long long tablesize;
        unsigned int cutoff;
        /* Number of worker tasks */
        unsigned int ntasks[NUM_TASK_TYPES];
        /* Input files unread or partially read  */
        TaskFile *files;
        /* Total number of tables created */
        unsigned int ntablescreated;
        /* Unsorted tables */
        unsigned int nunsorted;
        wordtable *unsorted[MAX_TABLES];
        /* Sorted tables */
        unsigned int nsorted;
        wordtable *sorted[MAX_TABLES];
        /* Available tables */
        unsigned int navailable;
        wordtable *available[MAX_TABLES];
};

void maker_queue_setup (MakerQueue *mq, unsigned int nthreads);
void maker_queue_release (MakerQueue *mq);

/* File parsing tasks */

struct _TaskFile {
        TaskFile *next;
        const char *filename;
        /* File index */
        unsigned int idx;
        FILE *ifs;
        unsigned int close_on_delete;
        const unsigned char *cdata;
        unsigned long long csize;
        unsigned int scout;
        unsigned int has_reader;
        FastaReader reader;
};

TaskFile *task_file_new (const char *filename, unsigned int scout);
TaskFile *task_file_new_from_stream (FILE *ifs, const char *filename, unsigned int close_on_delete);
void task_file_delete (TaskFile *tf);
/* Frontend to mmap and FastaReader */
unsigned int task_file_read_nwords (TaskFile *tf, unsigned long long maxwords, unsigned int wordsize,
	/* Called as soon as the full sequence name is known */
	int (*start_sequence) (FastaReader *, void *),
	/* Called when the full sequence has been parsed */
	int (*end_sequence) (FastaReader *, void *),
	int (*read_character) (FastaReader *, unsigned int character, void *),
	int (*read_nucleotide) (FastaReader *, unsigned int nucleotide, void *),
	int (*read_word) (FastaReader *, unsigned long long word, void *),
	void *data);

/* Add new file task to queue (not thread-safe) */
void maker_queue_add_file (MakerQueue *mq, const char *filename);
/* Get smallest table */
wordtable *queue_get_smallest_table (MakerQueue *queue);
/* Get largest table */
wordtable *queue_get_largest_table (MakerQueue *queue);

wordtable *queue_get_sorted (MakerQueue *queue);
wordtable *queue_get_smallest_sorted (MakerQueue *queue);
wordtable *queue_get_mostavailable_sorted (MakerQueue *queue);

/* MMap scouting */

typedef struct _MapData MapData;
struct _MapData {
  const unsigned char *cdata;
  unsigned long long csize;
};

/* Create new thread that reads mmap into resident memory */
void scout_mmap (const unsigned char *cdata, unsigned long long csize);
void delete_scouts ();

#endif /* SEQUENCE_H_ */
