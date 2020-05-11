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
#include "word-table.h"

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

typedef struct _TaskRead TaskRead;
typedef struct _TaskCollateTables TaskCollateTables;
typedef struct _TaskCollateFiles TaskCollateFiles;
typedef struct _GT4LMQSource GT4LMQSource;

struct _GT4LMQSource {
  unsigned int file_idx;
  unsigned int block_idx;
  unsigned long long start;
  unsigned long long length;
  /* Subsequence data */
  unsigned int n_subseqs;
  GT4SubSequence *subseqs;
  /* For collation */
  unsigned int first_subseq;
  /* Bookkeeping */
  unsigned int size_subseqs;
};

struct _GT4ListMakerQueue {
  GT4Queue queue;

  /* Parameters */
  unsigned int wordlen;
  /* Bookkeeping */
  unsigned int n_files_waiting;
  unsigned int n_files_reading;
  unsigned int n_tables_collating;
  unsigned int n_running;
  /* KMer tables */        
  unsigned int n_free_s_tables;
  GT4WordTable **free_s_tables;
  unsigned int n_used_s_tables;
  GT4WordTable **used_s_tables;
  /* Temporary files */
  unsigned int n_tmp_files;
  char *tmp_files[4096];
  unsigned int n_final_files;
  char *final_files[256];
  
  /* Source data for indexing */
  unsigned int n_sources;
  GT4LMQSource sources[1024];
  unsigned int subseq_block_size;

  /* Debug */
  QValue tokens[NUM_ID_TOKENS];
};

struct _GT4ListMakerQueueClass {
  GT4QueueClass queue_class;
};

unsigned int gt4_listmaker_queue_get_type (void);

void maker_queue_setup (GT4ListMakerQueue *mq, unsigned int n_threads, unsigned int w_len, unsigned int n_tmp_tables, unsigned int tmp_table_size, unsigned int data_size);
void maker_queue_release (GT4ListMakerQueue *mq);

/* Add new file task to queue (not thread-safe) */
void maker_queue_add_file (GT4ListMakerQueue *mq, const char *filename, unsigned int stream, unsigned int file_idx);

/* Not thread-safe (Each source should have their own thread) */
unsigned int maker_queue_add_subsequence (GT4ListMakerQueue *mq, unsigned int source_idx, unsigned long long name_pos, unsigned int name_len);

/* Tasks */

struct _TaskRead {
  GT4Task task;
  AZObject *source;
  GT4FastaReader reader;
  unsigned int idx;
  void *data;
};

void gt4_task_read_setup (TaskRead *tr, const char *id, GT4Queue *queue, AZObject *source, unsigned int idx, unsigned int wordlen);
void gt4_task_read_release (TaskRead *tr);
TaskRead *task_read_new (GT4Queue *queue, const char *id, AZObject *source, unsigned int idx, unsigned int wordlen);
void task_read_delete (TaskRead *tr);

struct _TaskCollateTables {
  GT4Task task;
  unsigned int n_tables;
  GT4WordTable *tables[1];
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

#endif
