#ifndef __GT4_LISTMAKER_QUEUE_H__
#define __GT4_LISTMAKER_QUEUE_H__

/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 * 
 * Copyright (C) 2014-18 University of Tartu
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

typedef struct _GT4ListMakerQueue GT4ListMakerQueue;
typedef struct _GT4ListMakerQueueClass GT4ListMakerQueueClass;

#define GT4_TYPE_LISTMAKER_QUEUE gt4_listmaker_queue_get_type()

#include "fasta.h"
#include "sequence-file.h"
#include "sequence-stream.h"
#include "wordtable.h"

#include "queue.h"

#define MAX_TABLES 256

#define TASK_READ 0
#define TASK_SORT 1
#define TASK_COLLATE_TABLES 2
#define TASK_COLLATE_FILES 3
#define TASK_MERGE 4
#define NUM_TASK_TYPES 5

/*
 * Queue stub
 *
 */

#define NUM_ID_TOKENS 32

typedef union _QValue QValue;

union _QValue {
  unsigned long long lval;
  double dval;
};

/*
 * Possible tasks are:
 *   - read segment from FastA file and build table
 *   - sort table
 *   - merge two tables
 */

/* Task for parsing FastA file */
/* These form linked list */

typedef struct _TaskFile TaskFile;
typedef struct _TaskRead TaskRead;
typedef struct _TaskCollateTables TaskCollateTables;
typedef struct _TaskCollateFiles TaskCollateFiles;

struct _GT4ListMakerQueue {
        GT4Queue queue;

        /* Parameters */
        unsigned int wordlen;
        unsigned long long tablesize;
        unsigned int cutoff;
        /* Waiting files */
        TaskFile *files;

        unsigned int n_files_waiting;
        unsigned int n_files_reading;
        unsigned int n_tables_collating;

        unsigned int n_running;
        unsigned int n_free_s_tables;
        wordtable **free_s_tables;
        unsigned int n_used_s_tables;
        wordtable **used_s_tables;

        unsigned int n_tmp_files;
        char *tmp_files[4096];
        unsigned int n_final_files;
        char *final_files[256];

        /* Debug */
        QValue tokens[NUM_ID_TOKENS];
};

struct _GT4ListMakerQueueClass {
  GT4QueueClass queue_class;
};

unsigned int gt4_listmaker_queue_get_type (void);

void maker_queue_setup (GT4ListMakerQueue *mq, unsigned int n_threads, unsigned int w_len, unsigned int n_tmp_tables, unsigned int tmp_table_size);
void maker_queue_release (GT4ListMakerQueue *mq);

/* Add new file task to queue (not thread-safe) */
void maker_queue_add_file (GT4ListMakerQueue *mq, const char *filename, unsigned int stream);

/* Tasks */

struct _TaskRead {
  GT4Task task;
  AZObject *source;
  GT4FastaReader reader;
};

TaskRead *task_read_new (GT4ListMakerQueue *mq, AZObject *source);
void task_read_delete (TaskRead *tr);

struct _TaskCollateTables {
  GT4Task task;
  unsigned int n_tables;
  wordtable *tables[1];
};

TaskCollateTables *task_collate_tables_new (GT4ListMakerQueue *mq, unsigned int max_tables);
void task_collate_tables_delete (TaskCollateTables *tc);

struct _TaskCollateFiles {
  GT4Task task;
  unsigned int n_files;
  char *files[1];
};

TaskCollateFiles *task_collate_files_new (GT4ListMakerQueue *mq, unsigned int max_files);
void task_collate_files_delete (TaskCollateFiles *tc);

/* File parsing tasks */

struct _TaskFile {
        TaskFile *next;
        GT4SequenceFile *seqfile;
        GT4SequenceStream *stream;
        /* Sequence source */
        AZObject *source;
        unsigned int close_source : 1;
        /* File index */
        unsigned int idx;
        unsigned int scout;
        unsigned int has_reader;
        GT4FastaReader reader;
};

TaskFile *task_file_new (const char *filename, unsigned int scout);
TaskFile *task_file_new_from_source (AZObject *source, const char *filename, unsigned int close);
void task_file_delete (TaskFile *tf);
/* Frontend to mmap and GT4FastaReader */
unsigned int task_file_read_nwords (TaskFile *tf, unsigned long long maxwords, unsigned int wordsize,
	/* Called as soon as the full sequence name is known */
	int (*start_sequence) (GT4FastaReader *, void *),
	/* Called when the full sequence has been parsed */
	int (*end_sequence) (GT4FastaReader *, void *),
	int (*read_character) (GT4FastaReader *, unsigned int character, void *),
	int (*read_nucleotide) (GT4FastaReader *, unsigned int nucleotide, void *),
	int (*read_word) (GT4FastaReader *, unsigned long long word, void *),
	void *data);

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
