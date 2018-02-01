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
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <pthread.h>
#include <unistd.h>

#include "common.h"
#include "listmaker-queue.h"
#include "utils.h"
#include "fasta.h"
#include "wordtable.h"
#include "sequence.h"
#include "sequence-stream.h"
#include "sequence-source.h"
#include "version.h"

#define MAX_FILES 200

#define DEFAULT_NUM_THREADS 8
#define DEFAULT_TABLE_SIZE 500000000
#define DEFAULT_MAX_TABLES 32

#define BSIZE 10000000

#define MAX_MERGED_TABLES 1024

#define NUM_READ 0
#define TIME_READ 1
#define NUM_SORT 2
#define TIME_SORT 3
#define NUM_WRITE_TMP 4
#define TIME_WRITE_TMP 5
#define TIME_MERGE 6
#define TIME_FF 7

#define TMP_TABLE_SIZE (1024 * 1024)
#define NUM_TMP_TABLES (24 * 256)
#define TMP_MERGE_SIZE 128

/* Main thread loop */
static void process (GT4Queue *queue, unsigned int idx, void *arg);
static void process2 (GT4Queue *queue, unsigned int idx, void *arg);

/* Merge tables directly to disk */
static unsigned long long merge_write_multi_nofreq (wordtable *t[], unsigned int ntables, int ofile);
static unsigned int merge_write_multi (wordtable **t, unsigned int ntables, const char *filename, unsigned int cutoff);

/* */
int process_word (GT4FastaReader *reader, unsigned long long word, void *data);

/* Print usage and help menu */
void print_help (int exitvalue);

int debug = 0;
int ntables = 0;
const char *outputname = "out";
const char *tmpdir = ".";

int 
main (int argc, const char *argv[])
{
  const char *files[1024];
  unsigned int n_files = 0;
  unsigned int argidx, i;
  char *end;

  /* default values */
  unsigned int wordlength = 16;
  unsigned int cutoff = 1;
  unsigned int nthreads = 0;
  unsigned long long tablesize = 0;
  unsigned int stream = 0;

  /* parsing commandline arguments */
  for (argidx = 1; argidx < argc; argidx++) {

    if (!strcmp (argv[argidx], "-v") || !strcmp (argv[argidx], "--version")) {
      fprintf (stdout, "glistmaker version %d.%d (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_QUALIFIER);
      return 0;

    } else if (!strcmp (argv[argidx], "-h") || !strcmp (argv[argidx], "--help") || !strcmp (argv[argidx], "-?")) {
      print_help (0);
    } else if (!strcmp (argv[argidx], "-o") || !strcmp (argv[argidx], "--outputname")) {
      if (!argv[argidx + 1] || argv[argidx + 1][0] == '-') {
        fprintf (stderr, "Warning: No output name specified!\n");
        argidx += 1;
        continue;
      }
      outputname = argv[argidx + 1];
      argidx += 1;
    } else if (!strcmp (argv[argidx], "-w") || !strcmp (argv[argidx], "--wordlength")) {
      if (!argv[argidx + 1]  || argv[argidx + 1][0] == '-') {
        fprintf (stderr, "Warning: No word-length specified! Using the default value: %d.\n", wordlength);
        argidx += 1;
        continue;
      }
      wordlength = strtol (argv[argidx + 1], &end, 10);
      if (*end != 0) {
        fprintf (stderr, "Error: Invalid word-length: %s! Must be an integer.\n", argv[argidx + 1]);
        print_help (1);
      }
      argidx += 1;
    } else if (!strcmp (argv[argidx], "-c") || !strcmp (argv[argidx], "--cutoff")) {
      if (!argv[argidx + 1] || argv[argidx + 1][0] == '-') {
        fprintf (stderr, "Warning: No frequency cut-off specified! Using the default value: %d.\n", cutoff);
        argidx += 1;
        continue;
      }
      cutoff = strtol (argv[argidx + 1], &end, 10);
      if (*end != 0) {
        fprintf (stderr, "Error: Invalid frequency cut-off: %s! Must be an integer.\n", argv[argidx + 1]);
        print_help (1);
      }
      argidx += 1;
    } else if (!strcmp (argv[argidx], "--num_threads")) {
      if (!argv[argidx + 1]  || argv[argidx + 1][0] == '-') {
        fprintf (stderr, "Warning: No num-threads specified! Using the default value: %d.\n", DEFAULT_NUM_THREADS);
        argidx += 1;
        continue;
      }
      nthreads = strtol (argv[argidx + 1], &end, 10);
      if (*end != 0) {
        fprintf (stderr, "Error: Invalid num-threads: %s! Must be an integer.\n", argv[argidx + 1]);
        print_help (1);
      }
      argidx += 1;
    } else if (!strcmp (argv[argidx], "--max_tables")) {
      if (!argv[argidx + 1]  || argv[argidx + 1][0] == '-') {
        fprintf (stderr, "Warning: No max_tables specified! Using the default value: %d.\n", DEFAULT_MAX_TABLES);
        argidx += 1;
        continue;
      }
      ntables = strtol (argv[argidx + 1], &end, 10);
      if (*end != 0) {
        fprintf (stderr, "Error: Invalid max_tables: %s! Must be an integer.\n", argv[argidx + 1]);
        print_help (1);
      }
      argidx += 1;
    } else if (!strcmp (argv[argidx], "--table_size")) {
      if (!argv[argidx + 1]  || argv[argidx + 1][0] == '-') {
        fprintf (stderr, "Warning: No table-size specified! Using the default value: %llu.\n", (unsigned long long) DEFAULT_TABLE_SIZE);
        argidx += 1;
        continue;
      }
      tablesize = strtoll (argv[argidx + 1], &end, 10);
      if (*end != 0) {
        fprintf (stderr, "Error: Invalid table-size: %s! Must be an integer.\n", argv[argidx + 1]);
        print_help (1);
      }
      argidx += 1;
    } else if (!strcmp (argv[argidx], "--tmpdir")) {
      argidx += 1;
      if (argidx >= argc) print_help (1);
      tmpdir = argv[argidx];
    } else if (!strcmp (argv[argidx], "--stream")) {
      stream = 1;
    } else if (!strcmp (argv[argidx], "-D")) {
      debug += 1;
    } else {
      /* Input file */
      if (n_files >= 1024) continue;
      files[n_files++] = argv[argidx];
    }
  }

  /* debug_tables = debug; */
  
  if (!nthreads) {
    nthreads = DEFAULT_NUM_THREADS;
    if (nthreads > ((3 * n_files + 1) >> 1)) nthreads = (3 * n_files + 1) >> 1;
  }
  if (!tablesize) {
    tablesize = DEFAULT_TABLE_SIZE;
  }
  if (!ntables) {
    ntables = DEFAULT_MAX_TABLES;
    if (ntables > ((3 * n_files + 1) >> 1) + 1) ntables = ((3 * n_files + 1) >> 1) + 1;
  }
  if (ntables > MAX_TABLES) ntables = MAX_TABLES;

  /* checking parameter values */
  if (!n_files) {
    fprintf (stderr, "Error: No FastA/FastQ file specified!\n");
    print_help (1);
  }
  if (wordlength < 1 || wordlength > 32) {
    fprintf (stderr, "Error: Invalid word-length: %d!\n", wordlength);
    print_help (1);
  }
  if (cutoff < 1) {
    fprintf (stderr, "Error: Invalid frequency cut-off: %d! Must be positive.\n", cutoff);
    print_help (1);
  }
  if (strlen (outputname) > 200) {
    fprintf (stderr, "Error: Output name exceeds the 200 character limit.");
    return 1;
  }
  if (nthreads < 1) nthreads = 1;
  if (nthreads > 256) nthreads = 256;
  for (i = 0; i < n_files; i++) {
    if (!strcmp (files[i], "-")) continue;
    struct stat s;
    if (stat (files[i], &s)) {
      fprintf (stderr, "Error: Cannot stat %s\n", files[i]);
      exit (1);
    }
  }
  
  if (nthreads > 1) {
    /* CASE: SEVERAL THREADS */
    GT4ListMakerQueue mq;
    int rc;
    unsigned int finished = 0;

    maker_queue_setup (&mq, nthreads, wordlength, NUM_TMP_TABLES, TMP_TABLE_SIZE);

    for (i = 0; i < n_files; i++) {
      maker_queue_add_file (&mq, files[i], stream);
    }

    mq.tablesize = tablesize;
    mq.cutoff = cutoff;

    if (debug) {
      fprintf (stderr, "Num threads is %d\n", nthreads);
      fprintf (stderr, "Num tables is %d\n", ntables);
      fprintf (stderr, "Table size is %lld\n", mq.tablesize);
    }

    rc = gt4_queue_create_threads (&mq.queue, process2, &mq);
    if (rc) {
         fprintf (stderr, "ERROR; return code from pthread_create() is %d\n", rc);
         exit (-1);
    }

    process2 (&mq.queue, 0, &mq);

    while (!finished) {
      gt4_queue_lock (&mq.queue);
      gt4_queue_broadcast (&mq.queue);
      if (mq.queue.nthreads_running < 2) finished = 1;
      gt4_queue_unlock (&mq.queue);
      /* process (&mq.queue, 0, &mq); */
      sleep (1);

    }
    if (mq.nsorted > 0) {
      /* write the final list into a file */
      if (debug > 0) fprintf (stderr, "Writing list %s\n", outputname);
      wordtable_write_to_file (mq.sorted[0], outputname, cutoff);
    }

    if (debug) {
      fprintf (stderr, "Read %llu words at %.2f (%u words/s)\n", (unsigned long long) mq.tokens[NUM_READ].dval, mq.tokens[TIME_READ].dval, (unsigned int) (mq.tokens[NUM_READ].dval /  mq.tokens[TIME_READ].dval));
      fprintf (stderr, "Sort %llu words at %.2f (%u words/s)\n", (unsigned long long) mq.tokens[NUM_SORT].dval, mq.tokens[TIME_SORT].dval, (unsigned int) (mq.tokens[NUM_SORT].dval /  mq.tokens[TIME_SORT].dval));
      fprintf (stderr, "Write tmp %llu words at %.2f (%u words/s)\n", (unsigned long long) mq.tokens[NUM_WRITE_TMP].dval, mq.tokens[TIME_WRITE_TMP].dval, (unsigned int) (mq.tokens[NUM_WRITE_TMP].dval /  mq.tokens[TIME_WRITE_TMP].dval));
      fprintf (stderr, "Collate %.2f\n", mq.tokens[TIME_FF].dval);
      fprintf (stderr, "Merge %.2f\n", mq.tokens[TIME_MERGE].dval);
    }

    maker_queue_release (&mq);
  } else {
    /* CASE: ONE THREAD */
    wordtable *table, *temptable;
    int v;
    
    /* creating initial tables */
    table = wordtable_new (wordlength, 20000);
    temptable = wordtable_new (wordlength, 20000);
    
    for (i = 0; i < n_files; i++) {
      GT4FastaReader reader;
      AZObject *obj;

      temptable->wordlength = wordlength;

      if (!strcmp (files[i], "-")) {
  /* stdin */
  GT4SequenceStream *stream = gt4_sequence_stream_new_from_stream (stdin, 0);
  fasta_reader_init (&reader, wordlength, 1, GT4_SEQUENCE_STREAM_SEQUENCE_SOURCE_IMPLEMENTATION(stream), &stream->source_instance);
  obj = AZ_OBJECT (stream);
      } else {
  GT4SequenceFile *seqf = gt4_sequence_file_new (files[i], 0);
  gt4_sequence_file_map_sequence (seqf);
  fasta_reader_init (&reader, wordlength, 1, GT4_SEQUENCE_FILE_SEQUENCE_SOURCE_IMPLEMENTATION(seqf), &seqf->block.source_instance);
  obj = AZ_OBJECT (seqf);
      }

      /* reading words from FastA/FastQ */
      v = fasta_reader_read_nwords (&reader, 0xffffffffffffffffULL, NULL, NULL, NULL, NULL, process_word, (void *) temptable);
      if (v) return print_error_message (v);
      fasta_reader_release (&reader);

      /* radix sorting */
      wordtable_sort (temptable, 0);
      v = wordtable_find_frequencies (temptable);
      if (v) return print_error_message (v);


      /* merging two tables */
      if (i > 0) {
  v = wordtable_merge (table, temptable);
  if (v) return print_error_message (v);
      } else {
  wordtable *t = table;
  table = temptable;
  temptable = t;
      }

      /* empty the temporary table */
      wordtable_empty (temptable);

      az_object_unref (obj);
    }

    /* write the final list into a file */
    if (debug > 0) fprintf (stderr, "Writing list %s\n", outputname);
    if (wordtable_write_to_file (table, outputname, cutoff)) {
      fprintf (stderr, "Cannot write list to file\n");
    }
  }

  /*wordtable_delete (temptable);
  wordtable_delete (table);*/

  pthread_exit (NULL);
}

/* Starts and ends with queue locked */

static void
collate_tables (GT4ListMakerQueue *mq, TaskCollate *tc)
{
  double t0, t1;
  unsigned int i;
  unsigned long long nwritten;
  char *tmpname;
  int ofile;
  gt4_queue_unlock (&mq->queue);
  tmpname = (char *) malloc (strlen (tmpdir) + 256);
  sprintf (tmpname, "%s/GLM4_XXXXXX.list", tmpdir);
  ofile = mkstemps (tmpname, 5);
  if (debug) fprintf (stderr, "Collating %u tables to %s\n", tc->n_tables, tmpname);
  t0 = get_time ();
  nwritten = merge_write_multi_nofreq (tc->tables, tc->n_tables, ofile);
  gt4_queue_lock (&mq->queue);
  for (i = 0; i < tc->n_tables; i++) {
    wordtable_empty (tc->tables[i]);
    mq->free_s_tables[mq->n_free_s_tables++] = tc->tables[i];
  }
  close (ofile);
  t1 = get_time ();
  task_collate_delete (tc);
  mq->tokens[NUM_WRITE_TMP].dval += nwritten;
  mq->tokens[TIME_WRITE_TMP].dval += (t1 - t0);
}

/* Returning 1 means that we have to wait for other threads */

static unsigned int
read_table (GT4ListMakerQueue *mq, TaskRead *tr)
{
  double t0, t1, t2;
  wordtable *tbl;
  if (!mq->n_free_s_tables) {
    gt4_queue_add_task (&mq->queue, &tr->task, 0);
    return 1;
  }
  tbl = mq->free_s_tables[--mq->n_free_s_tables];
  mq->n_files_waiting -= 1;
  mq->n_files_reading += 1;
  gt4_queue_unlock (&mq->queue);
  tbl->wordlength = mq->wordlen;
  t0 = get_time ();
  fasta_reader_read_nwords (&tr->reader, TMP_TABLE_SIZE, NULL, NULL, NULL, NULL, process_word, tbl);
  t1 = get_time ();
  if (tbl->nwords) wordtable_sort (tbl, 0);
  t2 = get_time ();
  gt4_queue_lock (&mq->queue);
  mq->n_files_reading -= 1;
  mq->tokens[NUM_READ].dval += tbl->nwords;
  mq->tokens[TIME_READ].dval += (t1 - t0);
  mq->tokens[NUM_SORT].dval += tbl->nwords;
  mq->tokens[TIME_SORT].dval += (t2 - t1);
  if (tr->reader.in_eof) {
    task_read_delete (tr);
  } else {
    gt4_queue_add_task (&mq->queue, &tr->task, 0);
    mq->n_files_waiting += 1;
  }
  if (tbl->nwords) {
    mq->used_s_tables[mq->n_used_s_tables++] = tbl;
  } else {
    mq->free_s_tables[mq->n_free_s_tables++] = tbl;
  }
  if ((mq->n_used_s_tables >= TMP_MERGE_SIZE) || (!mq->n_files_reading && !mq->n_files_waiting)) {
    unsigned int n_tables, i;
    /* Create collation task */
    n_tables = mq->n_used_s_tables;
    if (n_tables > TMP_MERGE_SIZE) n_tables = TMP_MERGE_SIZE;
    TaskCollate *tc = task_collate_new (mq, n_tables);
    for (i = 0; i < n_tables; i++) {
      tc->tables[tc->n_tables++] = mq->used_s_tables[--mq->n_used_s_tables];
    }
    gt4_queue_add_task (&mq->queue, &tc->task, 0);
  }
  return 0;
}

static void
process2 (GT4Queue *queue, unsigned int idx, void *arg)
{
  GT4ListMakerQueue *mq = (GT4ListMakerQueue *) arg;
  unsigned int finished = 0;

  if (debug > 1) {
    gt4_queue_lock (queue);
    fprintf (stderr, "Thread %d started (total %d)\n", idx, queue->nthreads_running);
    gt4_queue_unlock (queue);
  }

  /* Do work */
  while (!finished) {
    GT4Task *task;
    unsigned int wait = 0;
    gt4_queue_lock (queue);
    if (!queue->tasks && !mq->n_running) {
      /* Finish */
      finished = 1;
    } else if (!queue->tasks) {
      /* Wait */
      wait = 1;
    } else {
      task = queue->tasks;
      queue->tasks = task->next;
      if (debug > 1) fprintf (stderr, "Thread %u task %u\n", idx, task->type);
      if (task->type == TASK_READ) {
        TaskRead *tr = (TaskRead *) task;
        mq->n_running += 1;
        wait = read_table (mq, tr);
        mq->n_running -= 1;
      } else if (task->type == TASK_COLLATE) {
        TaskCollate *tc = (TaskCollate *) task;
        mq->n_running += 1;
        collate_tables (mq, tc);
        mq->n_running -= 1;
      }
    }
    if (wait) {
      gt4_queue_wait (queue);
    } else {
      gt4_queue_broadcast (queue);
    }
    gt4_queue_unlock (queue);
  }
}

static void
process (GT4Queue *queue, unsigned int idx, void *arg)
{
  GT4ListMakerQueue *mq = (GT4ListMakerQueue *) arg;
  unsigned int finished = 0;
  double s_t, e_t, d_t;

  if (debug > 1) {
    gt4_queue_lock (queue);
    fprintf (stderr, "Thread %d started (total %d)\n", idx, queue->nthreads_running);
    gt4_queue_unlock (queue);
  }

  /* Do work */
  while (!finished) {
    unsigned int has_files, has_unsorted, has_unmerged, sorted_tables;
    /* Get exclusive lock on queue */
    pthread_mutex_lock (&mq->queue.mutex);
    if (debug > 1) fprintf (stderr, "Thread %d: FileTasks %u Unsorted %u Sorted %u\n", idx, mq->ntasks[TASK_READ], mq->nunsorted, mq->nsorted);
    
    has_files = mq->files || mq->ntasks[TASK_READ];
    has_unsorted = mq->nunsorted || mq->ntasks[TASK_SORT];
    has_unmerged = mq->ntasks[TASK_MERGE];
    sorted_tables = mq->nsorted + mq->ntasks[TASK_MERGE];

    /* If all files have been read and sorted and there is small enough number of sorted files */
    if (!has_files && !has_unsorted && sorted_tables && (sorted_tables <= MAX_MERGED_TABLES)) {
      if (!has_unmerged) {
        /* Merge to disk */
        wordtable *t[MAX_MERGED_TABLES];
        unsigned int ntables;
        char c[1024];
        ntables = 0;
        while (mq->nsorted) {
          t[ntables++] = queue_get_sorted (mq);
        }
        wordtable_build_filename (t[0], c, 1024, outputname);
        mq->ntasks[TASK_MERGE] += 1;
        if (debug) {
          unsigned int i;
          fprintf (stderr, "Merging %u tables: %s", ntables, t[0]->id);
          for (i = 1; i < ntables; i++) {
	    fprintf (stderr, ",%s", t[i]->id);
	  }
	  fprintf (stderr, " to %s\n", c);
	}
	/* Now we can release mutex */
	pthread_mutex_unlock (&mq->queue.mutex);
	/* merge_write (table, other, c, queue->cutoff); */
	if (merge_write_multi (t, ntables, c, mq->cutoff)) {
	fprintf (stderr, "Cannot write list to file\n");
      }
  pthread_mutex_lock (&mq->queue.mutex);
  mq->ntasks[TASK_MERGE] -= 1;
  pthread_cond_broadcast (&mq->queue.cond);
  pthread_mutex_unlock (&mq->queue.mutex);
  finished = 1;
      } else {
  /* Waiting merging to finish */
        if (debug > 1) fprintf (stderr, "Thread %d: Waiting merging to finish\n", idx);
        pthread_cond_wait (&mq->queue.cond, &mq->queue.mutex);
        pthread_mutex_unlock (&mq->queue.mutex);
      }
      continue;
    }
    if (mq->nsorted > 1) {
      /* Task 1 - merge sorted tables */
      wordtable *table, *other;
      int result;

      other = queue_get_smallest_sorted (mq);
      table = queue_get_mostavailable_sorted (mq);
      if (table->nwordslots < other->nwordslots) {
        wordtable *t = table;
        table = other;
        other = t;
      }
      mq->ntasks[TASK_MERGE] += 1;
      /* Now we can release mutex */
      pthread_mutex_unlock (&mq->queue.mutex);
      if (debug > 0) fprintf (stderr, "Thread %d: Merging tables %s (%llu/%llu) + %s (%llu/%llu) -> %s\n", idx, table->id, table->nwords, table->nwordslots, other->id, other->nwords, other->nwordslots, table->id);
      s_t = get_time ();
      result = wordtable_merge (table, other);
      e_t = get_time ();
      d_t = e_t - s_t;
      /* fixme: Error processing */
      if (result) {
        print_error_message (result);
      }
      /* Lock mutex */
      pthread_mutex_lock (&mq->queue.mutex);
      /* Add merged table to sorted list */
      mq->sorted[mq->nsorted++] = table;
      wordtable_empty (other);
      other->wordlength = mq->wordlen;
      mq->available[mq->navailable++] = other;
      mq->tokens[TIME_MERGE].dval += d_t;
      /* Release mutex */
      mq->ntasks[TASK_MERGE] -= 1;
      pthread_cond_broadcast (&mq->queue.cond);
      pthread_mutex_unlock (&mq->queue.mutex);
      if (debug > 0) fprintf (stderr, "Thread %d: Finished merging %s (%llu/%llu)\n", idx, table->id, table->nwords, table->nwordslots);
    } else if (mq->nunsorted > 0) {
      /* Task 2 - sort table */
      wordtable *table;
      int result;
      
      table = mq->unsorted[--mq->nunsorted];
      /* Now we can release mutex */
      mq->ntasks[TASK_SORT] += 1;
      gt4_queue_unlock (&mq->queue);
      if (debug > 0) fprintf (stderr, "Thread %d: Sorting table %s (%llu/%llu)\n", idx, table->id, table->nwords, table->nwordslots);
      s_t = get_time ();
      wordtable_sort (table, 0);
      e_t = get_time ();
      d_t = e_t - s_t;
      s_t = get_time ();
      result = wordtable_find_frequencies (table);
      e_t = get_time ();
      /* fixme: Error processing */
      if (result) {
        print_error_message (result);
      }
      /* Lock mutex */
      gt4_queue_lock (&mq->queue);
      /* Add sorted table to sorted list */
      mq->sorted[mq->nsorted++] = table;
      mq->tokens[TIME_SORT].dval += d_t;
      d_t = e_t - s_t;
      mq->tokens[TIME_FF].dval += d_t;
      /* Release mutex */
      mq->ntasks[TASK_SORT] -= 1;
      gt4_queue_broadcast (&mq->queue);
      gt4_queue_unlock (&mq->queue);
      if (debug > 0) fprintf (stderr, "Thread %d: Finished sorting %s (%llu/%llu)\n", idx, table->id, table->nwords, table->nwordslots);
    } else if (mq->files && (mq->ntasks[TASK_READ] < MAX_FILES) && (mq->navailable || (mq->ntablescreated < ntables))) {
      /* Task 3 - read input file */
      TaskFile *task;
      wordtable *table;
      int result;
      unsigned long long readsize;

      task = mq->files;
      mq->files = task->next;
      mq->ntasks[TASK_READ] += 1;
      if (mq->navailable > 0) {
        /* Has to create new word table */
        table = queue_get_largest_table (mq);
      } else {
        table = wordtable_new (mq->wordlen, 10000000);
        table->wordlength = mq->wordlen;
        mq->ntablescreated += 1;
        if (debug > 0) fprintf (stderr, "Thread %d: Created table %s\n", idx, table->id);
      }
      /* Now we can release mutex */
      pthread_mutex_unlock (&mq->queue.mutex);
      
      /* Process file reader task */
      if (debug > 1) fprintf (stderr, "Thread %d: Processign file %s (%llu) -> %s (%llu)\n", idx, task->seqfile->path, (unsigned long long) task->reader.cpos, table->id, table->nwordslots);
      readsize = (mq->tablesize < table->nwordslots) ? table->nwordslots : mq->tablesize;
      if (debug > 0) fprintf (stderr, "Thread %d: Reading %lld bytes from %s, position %llu/%llu\n", idx, readsize, task->seqfile->path, (unsigned long long) task->reader.cpos, (unsigned long long) task->seqfile->block.csize);
      s_t = get_time ();
      result = task_file_read_nwords (task, readsize, mq->wordlen, NULL, NULL, NULL, NULL, process_word, table);
      e_t = get_time ();
      d_t = e_t - s_t;
      if (result) {
        /* fixme: Error processing */
        print_error_message (result);
      }
      /* Lock mutex */
      pthread_mutex_lock (&mq->queue.mutex);
      /* Add generated table to unsorted list */
      mq->unsorted[mq->nunsorted++] = table;
      /* fixme: We ignore error here - maybe should quit */
      if (result || task->reader.in_eof) {
        /* Finished this task */
        if (debug > 0) fprintf (stderr, "Thread %d: FastaReader for %s finished\n", idx, task->seqfile->path);
        task_file_delete (task);
      } else {
        /* Reshedule task */
        task->next = mq->files;
        mq->files = task;
      }
      mq->ntasks[TASK_READ] -= 1;
      mq->tokens[TIME_READ].dval += d_t;
      /* Release mutex */
      pthread_cond_broadcast (&mq->queue.cond);
      pthread_mutex_unlock (&mq->queue.mutex);
      if (debug > 0) fprintf (stderr, "Thread %d: Finished reading %s (%llu/%llu)\n", idx, table->id, table->nwords, table->nwordslots);
    } else if (!has_files && !mq->nunsorted && (mq->nsorted < 2)) {
      /* Nothing to do */
      /* Release mutex */
      pthread_cond_broadcast (&mq->queue.cond);
      pthread_mutex_unlock (&mq->queue.mutex);
      finished = 1;
    } else {
      if (debug > 1) fprintf (stderr, "Thread %d: Waiting\n", idx);
      /* Release mutex */
      /* pthread_mutex_unlock (&mq->queue.mutex); */
      /* fixme: Semaphore */
      /* sleep (1); */
      pthread_cond_wait (&mq->queue.cond, &mq->queue.mutex);
      pthread_mutex_unlock (&mq->queue.mutex);
    }
  }

  /* Exit if everything is done */
  if (debug) {
    gt4_queue_lock (queue);
    if (debug > 1) fprintf (stderr, "Thread %u exiting (remaining %d)\n", idx, queue->nthreads_running);
    gt4_queue_unlock (queue);
  }
}

int 
process_word (GT4FastaReader *reader, unsigned long long word, void *data)
{
  wordtable *table = (wordtable *) data;
#if 1
  wordtable_add_word_nofreq (table, word, reader->wordlength);
#endif
  return 0;
}

#define TMP_BUF_SIZE (256 * 12)

static unsigned long long
merge_write_multi_nofreq (wordtable *tables[], unsigned int ntables_in, int ofile)
{
  wordtable *t[MAX_MERGED_TABLES];
  unsigned long long i[MAX_MERGED_TABLES];
  unsigned int ntables, j;
  unsigned long long word;

  unsigned char b[TMP_BUF_SIZE];
  unsigned int bp = 0;

  GT4ListHeader h;

  h.code = GT4_LIST_CODE;
  h.version_major = VERSION_MAJOR;
  h.version_minor = VERSION_MINOR;
  h.wordlength = tables[0]->wordlength;
  h.nwords = 0;
  h.totalfreq = 0;
  h.padding = sizeof (GT4ListHeader);

  write (ofile, &h, sizeof (GT4ListHeader));

  ntables = 0;
  for (j = 0; j < ntables_in; j++) {
    if (!tables[j]->nwords) continue;
    t[ntables] = tables[j];
    i[ntables] = 0;
    ntables += 1;
  }

  word = 0xffffffffffffffff;
  /* Find smallest word */
  for (j = 0; j < ntables; j++) {
    if (t[j]->words[0] < word) word = t[j]->words[0];
  }
    
  while (ntables) {
    unsigned long long next = 0xffffffffffffffff;
    unsigned int freq = 0;
    /* Count freq */
    j = 0;
    while (j < ntables) {
      __builtin_prefetch (&t[j]->words[i[j] + 4], 0, 0);
      if (t[j]->words[i[j]] == word) {
        freq += 1;
        i[j] += 1;
        if (i[j] >= t[j]->nwords) {
          ntables -= 1;
          if (ntables > 0) {
            t[j] = t[ntables];
          } else {
            break;
          }
        }
      } else {
        if (t[j]->words[i[j]] < next) next = t[j]->words[i[j]];
        j += 1;
      }
    }
    memcpy (&b[bp], &word, 8);
    memcpy (&b[bp + 8], &freq, 4);
    bp += 12;
    if (bp >= TMP_BUF_SIZE) {
      write (ofile, b, bp);
      bp = 0;
    }
    //write (ofile, &word, 8);
    //write (ofile, &freq, 4);
    h.nwords += 1;
    h.totalfreq += freq;
    word = next;
  }
  if (bp > 0) {
    write (ofile, b, bp);
  }
  pwrite (ofile, &h, sizeof (GT4ListHeader), 0);
  if (debug > 1) {
    for (j = 0; j < ntables; j++) fprintf (stderr, " %llu", i[j]);
    fprintf (stderr, "\n");
  }

  return h.nwords;
}

static unsigned int
merge_write_multi (wordtable *t[], unsigned int ntables, const char *filename, unsigned int cutoff)
{
  unsigned long long nwords[MAX_MERGED_TABLES];
  unsigned long long i[MAX_MERGED_TABLES];
  unsigned int nfinished;

  GT4ListHeader h;

  FILE *ofs;
  char *b;
  unsigned int bp, j;

  unsigned long long word;
  unsigned int freq;
  double t_s, t_e;

  h.code = GT4_LIST_CODE;
  h.version_major = VERSION_MAJOR;
  h.version_minor = VERSION_MINOR;
  h.wordlength = t[0]->wordlength;
  h.nwords = 0;
  h.totalfreq = 0;
  h.padding = sizeof (GT4ListHeader);

  b = (char *) malloc (BSIZE + 12);
  
  ofs = fopen (filename, "w");
  if (!ofs) {
    fprintf (stderr, "Cannot open output file %s\n", filename);
    return 1;
  }

  t_s = get_time ();
  fwrite (&h, sizeof (GT4ListHeader), 1, ofs);

  bp = 0;
  nfinished = 0;
  for (j = 0; j < ntables; j++) {
    i[j] = 0;
    nwords[j] = t[j]->nwords;
    if (!nwords[j]) nfinished += 1;
  }
    
  memset (i, 0, sizeof (i));
  
  while (nfinished < ntables) {
    word = 0xffffffffffffffff;
    freq = 0;
    /* Find smalles word and total freq */
    for (j = 0; j < ntables; j++) {
      if (i[j] < nwords[j]) {
  /* This table is not finished */
  if (t[j]->words[i[j]] < word) {
    /* This table has smaller word */
    word = t[j]->words[i[j]];
    freq = t[j]->frequencies[i[j]];
  } else if (t[j]->words[i[j]] == word) {
    /* This table has equal word */
    freq += t[j]->frequencies[i[j]];
  }
  __builtin_prefetch (&t[j]->words[i[j]] + 16);
  __builtin_prefetch (&t[j]->frequencies[i[j]] + 16);
      }
    }
    /* Now we have word and freq */
    if (freq >= cutoff) {
      memcpy (b + bp, &word, 8);
      bp += 8;
      memcpy (b + bp, &freq, 4);
      bp += 4;
      if (bp >= BSIZE) {
  fwrite (b, 1, bp, ofs);
  bp = 0;
      }
      h.nwords += 1;
      h.totalfreq += freq;
    }
    /* Update pointers */
    for (j = 0; j < ntables; j++) {
      if (i[j] < nwords[j]) {
  /* This table is not finished */
  if (t[j]->words[i[j]] == word) {
    i[j] += 1;
    if (i[j] >= nwords[j]) {
      nfinished += 1;
    }
  }
      }
    }
  }
  if (bp) {
    fwrite (b, 1, bp, ofs);
  }
  fseek (ofs, 0, SEEK_SET);
  fwrite (&h, sizeof (GT4ListHeader), 1, ofs);
  fclose (ofs);
  t_e = get_time ();
  if (debug > 0) fprintf (stderr, "Writing %d tables with merging %.2f\n", ntables, t_e - t_s);
  
  free (b);

  return 0;
}

void 
print_help (int exitvalue)
{
  fprintf (stderr, "glistmaker version %u.%u (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_QUALIFIER);
  fprintf (stderr, "Usage: glistmaker <INPUTFILES> [OPTIONS]\n");
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "    -v, --version           - print version information and exit\n");
  fprintf (stderr, "    -h, --help              - print this usage screen and exit\n");
  fprintf (stderr, "    -w, --wordlength NUMBER - specify index wordsize (1-32) (default 16)\n");
  fprintf (stderr, "    -c, --cutoff NUMBER     - specify frequency cut-off (default 1)\n");
  fprintf (stderr, "    -o, --outputname STRING - specify output name (default \"out\")\n");
  fprintf (stderr, "    --num_threads           - number of threads the program is run on (default MIN(8, num_input_files))\n");
  fprintf (stderr, "    --max_tables            - maximum number of temporary tables (default MAX(num_threads, 2))\n");
  fprintf (stderr, "    --table_size            - maximum size of the temporary table (default 500000000)\n");
  fprintf (stderr, "    --stream                - read files as streams instead of memory-mapping (slower but uses less virtual memory)\n");
  fprintf (stderr, "    -D                      - increase debug level\n");
  exit (exitvalue);
}

