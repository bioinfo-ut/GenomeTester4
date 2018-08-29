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

#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <pthread.h>
#include <unistd.h>

#include "common.h"
#include "listmaker-queue.h"
#include "utils.h"
#include "fasta.h"
#include "sequence.h"
#include "sequence-stream.h"
#include "sequence-source.h"
#include "set-operations.h"
#include "version.h"
#include "word-list-stream.h"
#include "word-table.h"

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
#define NUM_TMP_TABLES (32 * 128)
#define TMP_MERGE_SIZE 64
#define FILE_MERGE_SIZE 32

/* Main thread loop */
static void process (GT4Queue *queue, unsigned int idx, void *arg);

/* Merge tables directly to disk */
static unsigned long long merge_write_multi_nofreq (GT4WordTable *t[], unsigned int ntables, int ofile);

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
      fprintf (stdout, "glistmaker version %u.%u.%u (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_MICRO, VERSION_QUALIFIER);
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
      if ((argv[argidx][0] == '-') && argv[argidx][1]) {
        print_help (1);
      }
      if (n_files >= 1024) continue;
      files[n_files++] = argv[argidx];
    }
  }

  /* debug_tables = debug; */
  
  if (!nthreads) {
    nthreads = DEFAULT_NUM_THREADS;
  }
  if (!tablesize) {
    tablesize = DEFAULT_TABLE_SIZE;
  }
  if (!ntables) {
    ntables = DEFAULT_MAX_TABLES;
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
  
  if (nthreads > 0) {
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

    rc = gt4_queue_create_threads (&mq.queue, process, &mq);
    if (rc) {
         fprintf (stderr, "ERROR; return code from pthread_create() is %d\n", rc);
         exit (-1);
    }

    process (&mq.queue, 0, &mq);

    while (!finished) {
      gt4_queue_lock (&mq.queue);
      gt4_queue_broadcast (&mq.queue);
      if (mq.queue.nthreads_running < 2) finished = 1;
      gt4_queue_unlock (&mq.queue);
      /* process (&mq.queue, 0, &mq); */
      sleep (1);
    }

    if (mq.n_final_files > 0) {
      AZObject *objs[4096];
      char c[1024];
      unsigned int i;
      GT4ListHeader header;
      int ofile;
      for (i = 0; i < mq.n_final_files; i++) {
        objs[i] = (AZObject *) gt4_word_list_stream_new (mq.final_files[i], VERSION_MAJOR);
      }
      snprintf (c, 1024, "%s_%u.list", outputname, wordlength);
      ofile = creat (c, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
      gt4_write_union (objs, mq.n_final_files, 1, ofile, &header);
      close (ofile);
      for (i = 0; i < mq.n_final_files; i++) {
        gt4_word_list_stream_delete (GT4_WORD_LIST_STREAM (objs[i]));
        unlink (mq.final_files[i]);
      }
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
#if 0
    /* CASE: ONE THREAD */
    GT4WordTable *table, *temptable;
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
  GT4WordTable *t = table;
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
#endif
  }

  /*wordtable_delete (temptable);
  wordtable_delete (table);*/

  pthread_exit (NULL);
}

/* Starts and ends with queue locked */

static void
collate_files (GT4ListMakerQueue *mq, TaskCollateFiles *tc)
{
  AZObject *objs[256];
  unsigned int i;
  char *tmpname;
  int ofile;
  GT4ListHeader header;

  gt4_queue_unlock (&mq->queue);
  tmpname = (char *) malloc (strlen (tmpdir) + 256);
  sprintf (tmpname, "%s/GLM4_F_XXXXXX.list", tmpdir);
  for (i = 0; i < tc->n_files; i++) {
    objs[i] = (AZObject *) gt4_word_list_stream_new (tc->files[i], VERSION_MAJOR);
  }
  ofile = mkstemps (tmpname, 5);
  if (debug) fprintf (stderr, "Collating %u files to %s\n", tc->n_files, tmpname);
  gt4_write_union (objs, tc->n_files, 1, ofile, &header);
  close (ofile);
  for (i = 0; i < tc->n_files; i++) {
    gt4_word_list_stream_delete (GT4_WORD_LIST_STREAM (objs[i]));
    unlink (tc->files[i]);
  }
  gt4_queue_lock (&mq->queue);

  mq->final_files[mq->n_final_files++] = tmpname;

  if (mq->n_final_files >= 16) {
    /* Create new collation task */
    TaskCollateFiles *tc = task_collate_files_new (mq, 16);
    for (i = 0; i < 16; i++) {
      tc->files[tc->n_files++] = mq->final_files[--mq->n_final_files];
    }
    gt4_queue_add_task (&mq->queue, &tc->task, 0);
  }

  task_collate_files_delete (tc);
}

static void
collate_tables (GT4ListMakerQueue *mq, TaskCollateTables *tc)
{
  double t0, t1;
  unsigned int i;
  unsigned long long nwritten;
  char *tmpname;
  int ofile;

  mq->n_tables_collating += 1;

  gt4_queue_unlock (&mq->queue);
  tmpname = (char *) malloc (strlen (tmpdir) + 256);
  sprintf (tmpname, "%s/GLM4_T_XXXXXX.list", tmpdir);
  t0 = get_time ();
  ofile = mkstemps (tmpname, 5);
  if (debug) fprintf (stderr, "Collating %u tables to %s\n", tc->n_tables, tmpname);
  nwritten = merge_write_multi_nofreq (tc->tables, tc->n_tables, ofile);
  close (ofile);
  t1 = get_time ();
  gt4_queue_lock (&mq->queue);

  mq->n_tables_collating -= 1;

  for (i = 0; i < tc->n_tables; i++) {
    wordtable_empty (tc->tables[i]);
    mq->free_s_tables[mq->n_free_s_tables++] = tc->tables[i];
  }
  mq->tmp_files[mq->n_tmp_files++] = tmpname;
  if ((mq->n_tmp_files >= FILE_MERGE_SIZE) || (!mq->n_files_reading && !mq->n_files_waiting && !mq->n_tables_collating)) {
    unsigned int n_files, i;
    /* Create collation task */
    n_files = mq->n_tmp_files;
    if (n_files > FILE_MERGE_SIZE) n_files = FILE_MERGE_SIZE;
    TaskCollateFiles *tc = task_collate_files_new (mq, n_files);
    for (i = 0; i < n_files; i++) {
      tc->files[tc->n_files++] = mq->tmp_files[--mq->n_tmp_files];
    }
    gt4_queue_add_task (&mq->queue, &tc->task, 0);
  }
  mq->tokens[NUM_WRITE_TMP].dval += nwritten;
  mq->tokens[TIME_WRITE_TMP].dval += (t1 - t0);
  task_collate_tables_delete (tc);
}

/* Returning 1 means that we have to wait for other threads */

static unsigned int
read_table (GT4ListMakerQueue *mq, TaskRead *tr)
{
  double t0, t1, t2;
  GT4WordTable *tbl;
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
    TaskCollateTables *tc = task_collate_tables_new (mq, n_tables);
    for (i = 0; i < n_tables; i++) {
      tc->tables[tc->n_tables++] = mq->used_s_tables[--mq->n_used_s_tables];
    }
    gt4_queue_add_task (&mq->queue, &tc->task, 0);
  }
  mq->tokens[NUM_READ].dval += tbl->nwords;
  mq->tokens[TIME_READ].dval += (t1 - t0);
  mq->tokens[NUM_SORT].dval += tbl->nwords;
  mq->tokens[TIME_SORT].dval += (t2 - t1);
  return 0;
}

static void
process (GT4Queue *queue, unsigned int idx, void *arg)
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
    unsigned int wait = 0;
    gt4_queue_lock (queue);
    if (!queue->tasks && !mq->n_running) {
      /* Finish */
      finished = 1;
    } else if (!queue->tasks) {
      /* Wait */
      wait = 1;
    } else {
      GT4Task *task;
      for (task = queue->tasks; task; task = task->next) {
        if (task->type == TASK_READ) {
          if (mq->n_free_s_tables) {
            gt4_queue_remove_task (queue, task, 0);
            mq->n_running += 1;
            wait = read_table (mq, (TaskRead *) task);
            mq->n_running -= 1;
            break;
          }
        } else if (task->type == TASK_COLLATE_TABLES) {
          gt4_queue_remove_task (queue, task, 0);
          mq->n_running += 1;
          collate_tables (mq, (TaskCollateTables *) task);
          mq->n_running -= 1;
          break;
        } else if (task->type == TASK_COLLATE_FILES) {
          gt4_queue_remove_task (queue, task, 0);
          mq->n_running += 1;
          collate_files (mq, (TaskCollateFiles *) task);
          mq->n_running -= 1;
          break;
        }
      }
      if (!task) wait = 1;
    }
    if (wait) {
      gt4_queue_wait (queue);
    } else {
      gt4_queue_broadcast (queue);
    }
    gt4_queue_unlock (queue);
  }
}

int 
process_word (GT4FastaReader *reader, unsigned long long word, void *data)
{
  GT4WordTable *table = (GT4WordTable *) data;
  wordtable_add_word_nofreq (table, word, reader->wordlength);
  return 0;
}

#define TMP_BUF_SIZE (256 * 12)

static unsigned long long
merge_write_multi_nofreq (GT4WordTable *tables[], unsigned int ntables_in, int ofile)
{
  GT4WordTable *t[MAX_MERGED_TABLES];
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
  h.list_start = sizeof (GT4ListHeader);

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
            i[j] = i[ntables];
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
    h.nwords += 1;
    h.totalfreq += freq;
    word = next;
  }
  if (bp > 0) {
    write (ofile, b, bp);
  }
  pwrite (ofile, &h, sizeof (GT4ListHeader), 0);
  if (debug) {
    fprintf (stderr, "Words %llu, unique %llu\n", h.totalfreq, h.nwords);
  }

  return h.nwords;
}

void 
print_help (int exitvalue)
{
  fprintf (stderr, "glistmaker version %u.%u.%u (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_MICRO, VERSION_QUALIFIER);
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
  fprintf (stderr, "    --tmpdir                - directory for temporary files (may need an order of magnitude more space than the size of the final list)\n");
  fprintf (stderr, "    --stream                - read files as streams instead of memory-mapping (slower but uses less virtual memory)\n");
  fprintf (stderr, "    -D                      - increase debug level\n");
  exit (exitvalue);
}

