#define __GMERCOUNTER_C__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <sys/mman.h>

#include "utils.h"
#include "trie.h"
#include "sequence.h"
#include "fasta.h"
#include "queue.h"
#include "wordmap.h"
#include "database.h"

#define MAX_LINES 10000000000
#define MAX_FILESIZE 10000000000

#define BLOCK_SIZE (1024 * 1024 * 10)
#define DEFAULT_NUM_THREADS 24
#define DEFAULT_NUM_TABLES 24

#define MAX_KMERS 25

unsigned int debug = 0;

typedef struct _TaskTable TaskTable;
struct _TaskTable {
  TaskTable *next;
  unsigned int nwords;
  unsigned long long *words;
  unsigned int *alleles;
};

TaskTable *task_table_new ();
void task_table_free (TaskTable *tt);

typedef struct _SNPQueue SNPQueue;
struct _SNPQueue {
  Queue queue;
  /* Files to process */
  unsigned int nfiles;
  TaskFile *files;
  /* Free tables */
  TaskTable *free_tables;
  /* Filled tables */
  TaskTable *full_tables;
  /* Data */
  KMerDB *db;
};

/* Main thread loop */
static void *process (void *arg);
static int start_sequence (FastaReader *reader, void *data);
static int read_word_2 (FastaReader *reader, unsigned long long word, void *data);
static int compare_counts (const void *lhs, const void *rhs);
static unsigned int get_pair_median (KMerDB *db);

static void
print_usage (FILE *ofs) {
  fprintf (ofs, "Usage:\n");
  fprintf (ofs, "  gmer_counter ARGUMENTS SEQUENCES...\n");
  fprintf (ofs, "Arguments:\n");
  fprintf (ofs, "    -db DATABASE     - SNP/KMER database file\n");
  fprintf (ofs, "    -dbb DBBINARY    - binary database file\n");
  fprintf (ofs, "    -w FILENAME      - write binary database to file\n");
  fprintf (ofs, "    -32              - use 32-bit integeres for counts (default 16-bit)\n");
  fprintf (ofs, "    --max_kmers NUM  - maximum number of kmers per node\n");
  fprintf (ofs, "    --header         - print header row\n");
  fprintf (ofs, "    --total          - print the total number of kmers per node\n");
  fprintf (ofs, "    --unique         - print the number of nonzero kmers per node\n");
  fprintf (ofs, "    --kmers          - print individual kmer counts (default if no other output)\n");
  fprintf (ofs, "    --distribution NUM  - print kmer distribution (up to given number)\n");
  fprintf (ofs, "    --num_threads    - number of worker threads (default %u)\n", DEFAULT_NUM_THREADS);
  fprintf (ofs, "    --low_memory     - optimize for low memory usage\n");
  fprintf (ofs, "    -D               - increase debug level\n");
}

int
main (int argc, const char *argv[])
{
  const char *db_name = NULL;
  const char *dbb = NULL;
  const char *wdb = NULL;
  unsigned int max_kmers_per_node = 1000000000;
  unsigned int header = 0, total = 0, unique = 0, kmers = 0, distro = 0, big = 0, dm = 0;
  unsigned int lowmem = 0;
  unsigned int nseqs = 0;
  const char *seqnames[1024];
  unsigned long long i;

  unsigned int nthreads = DEFAULT_NUM_THREADS;
  SNPQueue snpq;
  pthread_t threads[256];

  KMerDB db;

  
  for (i = 1; i < argc; i++) {
    if (!strcmp (argv[i], "-db")) {
      /* Database */
      i += 1;
      if (i >= argc) {
        print_usage (stderr);
        exit (1);
      }
      db_name = argv[i];
    } else if (!strcmp (argv[i], "-dbb")) {
      /* Binary database */
      i += 1;
      if (i >= argc) {
        print_usage (stderr);
        exit (1);
      }
      dbb = argv[i];
    } else if (!strcmp (argv[i], "-w")) {
      /* Write database */
      i += 1;
      if (i >= argc) {
        print_usage (stderr);
        exit (1);
      }
      wdb = argv[i];
    } else if (!strcmp (argv[i], "--max_kmers")) {
      i += 1;
      if (i >= argc) {
        print_usage (stderr);
        exit (1);
      }
      max_kmers_per_node = strtol (argv[i], NULL, 10);
    } else if (!strcmp (argv[i], "--header")) {
      header = 1;
    } else if (!strcmp (argv[i], "--total")) {
      total = 1;
    } else if (!strcmp (argv[i], "--unique")) {
      unique = 1;
    } else if (!strcmp (argv[i], "--kmers")) {
      kmers = 1;
    } else if (!strcmp (argv[i], "-32")) {
      big = 1;
    } else if (!strcmp (argv[i], "--double_median")) {
      dm = 1;
    } else if (!strcmp (argv[i], "--distribution")) {
      i += 1;
      if (i >= argc) {
        print_usage (stderr);
        exit (1);
      }
      distro = strtol (argv[i], NULL, 10);
    } else if (!strcmp (argv[i], "--num_threads")) {
      i += 1;
      if (i >= argc) {
        print_usage (stderr);
        exit (1);
      }
      nthreads = strtol (argv[i], NULL, 10);
    } else if (!strcmp (argv[i], "--low_memory")) {
      lowmem = 1;
    } else if (!strcmp (argv[i], "-D")) {
      /* Debug */
      debug += 1;
    } else {
      if (nseqs >= 1024) {
        fprintf (stderr, "Maximum number of input sequence files is 1024\n");
        exit (1);
      }
      seqnames[nseqs++] = argv[i];
    }
  }

  if (!nseqs && !wdb) {
    fprintf (stderr, "Nothing to do!\n");
    print_usage (stderr);
    exit (1);
  }
  if (db_name && dbb) {
    fprintf (stderr, "Both text and binary database specifed\n");
    print_usage (stderr);
    exit (1);
  }
  if (dbb && wdb) {
    fprintf (stderr, "Binary database read and written\n");
    print_usage (stderr);
    exit (1);
  }
  if (!total && !unique && !distro) {
    kmers = 1;
  }
  if (distro > 65536) {
    distro = 65536;
  }

  memset (&db, 0, sizeof db);

  if (db_name) {
    /* Read text database */
    const unsigned char *cdata;
    unsigned long long csize;
    cdata = gt4_mmap (db_name, &csize);
    if (!cdata) {
      fprintf (stderr, "Cannot mmap database file %s\n", db_name);
      exit (1);
    }
    if (!lowmem) scout_mmap (cdata, csize);
    if (!read_db_from_text (&db, cdata, csize, max_kmers_per_node, (big) ? 32 : 16)) {
      fprintf (stderr, "Cannot read text database %s\n", dbb);
      exit (1);
    }
  }

  if (dbb) {
    /* Read binary database */
    const unsigned char *cdata;
    unsigned long long csize;

    if (debug) fprintf (stderr, "Loading binary database %s\n", dbb);
    cdata = gt4_mmap (dbb, &csize);
    if (!cdata) {
      fprintf (stderr, "Cannot mmap %s\n", dbb);
      exit (1);
    }
    if (!lowmem) scout_mmap (cdata, csize);
    if (!read_database_from_binary (&db, cdata, csize)) {
      fprintf (stderr, "Cannot read binary database %s\n", dbb);
      exit (1);
    }
    if (debug) fprintf (stderr, "Finished loading binary database\n");
  }

  if (wdb) {
    /* Write binary database */
    FILE *ofs;
    if (debug) {
      fprintf (stderr, "Writing binary database to %s\n", wdb);
    }
    ofs = (fopen (wdb, "w+"));
    if (!ofs) {
      fprintf (stderr, "Cannot open %s for writing\n", wdb);
      exit (1);
    }
    write_db_to_file (&db, ofs, 0);
    fclose (ofs);
    if (debug) {
      fprintf (stderr, "Done\n");
    }
  }

  if (nseqs > 0) {
    /* Read files */
    /* Initialize queue */
    memset (&snpq, 0, sizeof (snpq));
    snpq.db = &db;
    for (i = 0; i < nseqs; i++) {
      TaskFile *tf;
      if (!strcmp (seqnames[i], "-")) {
        tf = task_file_new_from_stream (stdin, seqnames[i], 0);
      } else {
        tf = task_file_new (seqnames[i], !lowmem);
      }
      tf->next = snpq.files;
      snpq.files = tf;
      snpq.nfiles += 1;
    }
    for (i = 0; i < DEFAULT_NUM_TABLES; i++) {
      TaskTable *tt = task_table_new ();
      tt->next = snpq.free_tables;
      snpq.free_tables = tt;
    }
    /* Initialize main mutex and cond */
    queue_init (&snpq.queue);
    /* Lock the mutex */
    queue_lock (&snpq.queue);
    snpq.queue.nthreads = 0;
    for (i = 1; i < (int) nthreads; i++){
      int rc;
      if (debug > 1) fprintf (stderr, "Creating thread %llu\n", i);
      rc = pthread_create (&threads[i], NULL, process, &snpq);
      if (rc) {
        fprintf (stderr, "ERROR; return code from pthread_create() is %d\n", rc);
        exit (-1);
      }
    }
    queue_unlock (&snpq.queue);
    process (&snpq);
    queue_lock (&snpq.queue);
    while (snpq.queue.nthreads > 0) {
      queue_wait (&snpq.queue);
    }
    queue_unlock (&snpq.queue);
    queue_finalize (&snpq.queue);

    if (debug) {
      fprintf (stderr, "Finished reading files\n");
    }
    
    if (dm) {
      unsigned int med = get_pair_median (&db);
      fprintf (stdout, "#PairMedian\t%u\n", med);
    }
    
    if (header) {
      fprintf (stdout, "NODE\tN_KMERS");
      if (total) fprintf (stdout, "\tTOTAL");
      if (unique) fprintf (stdout, "\tUNIQUE");
      if (kmers) fprintf (stdout, "\tKMERS");
      if (distro) fprintf (stdout, "\tDISTRIBUTION");
      fprintf (stdout, "\n");
    }
  
    for (i = 0; i < db.n_nodes; i++) {
      unsigned int j;
      fprintf (stdout, "%s\t%u", db.names + db.nodes[i].name, db.nodes[i].nkmers);
      if (total) {
        unsigned long long total = 0;
        for (j = 0; j < db.nodes[i].nkmers; j++) {
          if (db.count_bits == 16) {
            total += db.kmers_16[db.nodes[i].kmers + j];
          } else {
            total += db.kmers_32[db.nodes[i].kmers + j];
          }
        }
        fprintf (stdout, "\t%llu", total);
      }
      if (unique) {
        unsigned int uniq = 0;
        for (j = 0; j < db.nodes[i].nkmers; j++) {
          if (db.count_bits == 16) {
            if (db.kmers_16[db.nodes[i].kmers + j]) uniq += 1;
          } else {
            if (db.kmers_16[db.nodes[i].kmers + j]) uniq += 1;
          }
        }
        fprintf (stdout, "\t%u", uniq);
      }
      if (kmers) {
        for (j = 0; j < db.nodes[i].nkmers; j++) {
          if (db.count_bits == 16) {
            fprintf (stdout, "\t%u", db.kmers_16[db.nodes[i].kmers + j]);
          } else {
            fprintf (stdout, "\t%u", db.kmers_32[db.nodes[i].kmers + j]);
          }
        }
      }
      if (distro) {
        static unsigned int c_len = 0;
        static unsigned int *c = NULL;
        if (c_len < db.nodes[i].nkmers) {
          c_len = c_len << 1;
          if (c_len < db.nodes[i].nkmers) c_len = db.nodes[i].nkmers;
          c = (unsigned int *) realloc (c, c_len * 4);
        }
        unsigned int current, count;
        if (db.count_bits == 16) {
          for (j = 0; j < db.nodes[i].nkmers; j++) c[j] = db.kmers_16[db.nodes[i].kmers + j];
        } else {
          memcpy (c, db.kmers_32 + db.nodes[i].kmers, db.nodes[i].nkmers * 4);
        }
        qsort (c, db.nodes[i].nkmers, 4, compare_counts);
        current = 0;
        j = 0;
        while (current <= distro) {
          count = 0;
          while ((j < db.nodes[i].nkmers) && (c[j] == current)) {
            count += 1;
            j += 1;
          }
          fprintf (stdout, "\t%u", count);
          current += 1;
        }
      }
      fprintf (stdout, "\n");
    }
  }
  
  return 0;
}

static void *
process (void *arg)
{
  SNPQueue *snpq;
  KMerDB *db;
  int idx = -1;
  unsigned int finished;

  snpq = (SNPQueue *) arg;
  db = snpq->db;

  finished = 0;
        
  /* Do work */
  while (!finished) {
    /* Get exclusive lock on queue */
    pthread_mutex_lock (&snpq->queue.mutex);
    if (idx < 0) {
      /* Thread is started, increase counter and get idx */
      idx = snpq->queue.nthreads;
      snpq->queue.nthreads += 1;
      if (debug > 1) fprintf (stderr, "Thread %d started (total %d)\n", idx, snpq->queue.nthreads);
    }
    if ((snpq->files) && (snpq->free_tables)) {
      TaskFile *tf;
      TaskTable *tt;
      tf = snpq->files;
      snpq->files = tf->next;
      tt = snpq->free_tables;
      snpq->free_tables = tt->next;
      pthread_mutex_unlock (&snpq->queue.mutex);
      /* Read words from file */
      tt->nwords = 0;
      if (debug > 0) fprintf (stderr, "Thread %d: reading file %s from %llu\n", idx, tf->filename, tf->reader.cpos);
      if (task_file_read_nwords (tf, BLOCK_SIZE, snpq->db->wordsize, start_sequence, NULL, NULL, NULL, read_word_2, tt)) {
        fprintf (stderr, "Cannot create FastaReader fro %s\n", tf->filename);
        exit (1);
      }
      if (debug > 1) fprintf (stderr, "Thread %d: finished reading %s at %llu\n", idx, tf->filename, tf->reader.cpos);
      pthread_mutex_lock (&snpq->queue.mutex);
      if (tf->reader.in_eof) {
        task_file_delete (tf);
        snpq->nfiles -= 1;
      } else {
        tf->next = snpq->files;
        snpq->files = tf;
      }
      tt->next = snpq->full_tables;
      snpq->full_tables = tt;
      pthread_cond_broadcast (&snpq->queue.cond);
      pthread_mutex_unlock (&snpq->queue.mutex);
    } else if (snpq->full_tables) {
      TaskTable *tt;
      unsigned int i;
      tt = (TaskTable *) snpq->full_tables;
      snpq->full_tables = tt->next;

      pthread_mutex_unlock (&snpq->queue.mutex);
      if (debug > 1) fprintf (stderr, "Thread %d: table lookup\n", idx);
      for (i = 0; i < tt->nwords; i++) {
        if (debug > 1) {
          if ((i % 10000) == 0) fprintf (stderr, ".");
        }
        tt->alleles[i] = trie_lookup (&snpq->db->trie, tt->words[i]);
      }
      if (debug > 1) fprintf (stderr, "Thread %d: finished lookup\n", idx);
      pthread_mutex_lock (&snpq->queue.mutex);

      /* fixme: Create separate task / mutex */
      for (i = 0; i <  tt->nwords; i++) {
        unsigned int code, node, kmer;
        code = tt->alleles[i];
        if (!code) continue;
        node = (code >> db->kmer_bits) - 1;
        if (node >= db->n_nodes) {
          fprintf (stderr, "DB inconsistency: Node index %u is bigger than the number of nodes %llu\n", node, db->n_nodes);
          break;
        }
        kmer = code & ((1 << db->kmer_bits) - 1);
        if (kmer >= db->nodes[node].nkmers) {
          fprintf (stderr, "DB inconsistency: KMer index %u is bigger than the number of kmers %u\n", kmer, db->nodes[node].nkmers);
          break;
        }
        /* Increase kmer count */
        if (db->count_bits == 16) {
          if (db->kmers_16[db->nodes[node].kmers + kmer] < 65535) db->kmers_16[db->nodes[node].kmers + kmer] += 1;
        } else {
          if (db->kmers_32[db->nodes[node].kmers + kmer] < 0xffffffff) db->kmers_32[db->nodes[node].kmers + kmer] += 1;
        }
      }
      tt->next = snpq->free_tables;
      snpq->free_tables = tt;
      pthread_cond_broadcast (&snpq->queue.cond);
      pthread_mutex_unlock (&snpq->queue.mutex);
    } else if (snpq->nfiles) {
      /* All tables are in processing, wait */
      if (debug > 1) fprintf (stderr, "Thread %d: Waiting\n", idx);
      pthread_cond_wait (&snpq->queue.cond, &snpq->queue.mutex);
      if (debug > 1) fprintf (stderr, "Thread %d: Woke up\n", idx);
      pthread_mutex_unlock (&snpq->queue.mutex);
    } else {
      /* Nothing to do */
      if (debug > 1) fprintf (stderr, "Thread %d: Exiting\n", idx);
      pthread_cond_broadcast (&snpq->queue.cond);
      pthread_mutex_unlock (&snpq->queue.mutex);
      finished = 1;
    }
  }
  /* Exit if everything is done */
  pthread_mutex_lock (&snpq->queue.mutex);
  snpq->queue.nthreads -= 1;
  if (debug > 1) fprintf (stderr, "Thread %u exiting (remaining %d)\n", idx, snpq->queue.nthreads);
  pthread_cond_broadcast (&snpq->queue.cond);
  pthread_mutex_unlock (&snpq->queue.mutex);

  return 0;
}

static int
start_sequence (FastaReader *reader, void *data)
{
  //TaskTable *tt = (TaskTable *) data;
  if (debug > 2) fprintf (stderr, "%s\n", reader->name);
  return 0;
}

static int
read_word_2 (FastaReader *reader, unsigned long long word, void *data)
{
  TaskTable *tt = (TaskTable *) data;
  tt->words[tt->nwords++] = word;
  return 0;
}

static int
compare_counts (const void *lhs, const void *rhs) {
  if (*((unsigned int *) lhs) < *((unsigned int *) rhs)) return -1;
  if (*((unsigned int *) lhs) == *((unsigned int *) rhs)) return 0;
  return 1;
}

static unsigned int
get_pair_median (KMerDB *db)
{
  unsigned long long i;
  unsigned int total, min, max, med, j;

  max = 0;
  min = 0xffffffff;
  total = 0;
  for (i = 0; i < db->n_nodes; i++) {
    total += db->nodes[i].nkmers / 2;
    for (j = 0; j < db->nodes[i].nkmers; j += 2) {
      unsigned int sum;
      if (db->count_bits == 16) {
        sum = db->kmers_16[db->nodes[i].kmers + j] + db->kmers_16[db->nodes[i].kmers + j + 1];
      } else {
        sum = db->kmers_32[db->nodes[i].kmers + j] + db->kmers_32[db->nodes[i].kmers + j + 1];
      }
      if (sum > max) max = sum;
      if (sum < min) min = sum;
    }
  }
  med = (min + max) / 2;
  while (max > min) {
    unsigned int above = 0, below = 0, equal;
    for (i = 0; i < db->n_nodes; i++) {
      for (j = 0; j < db->nodes[i].nkmers; j += 2) {
        unsigned int sum;
        if (db->count_bits == 16) {
          sum = db->kmers_16[db->nodes[i].kmers + j] + db->kmers_16[db->nodes[i].kmers + j + 1];
        } else {
          sum = db->kmers_32[db->nodes[i].kmers + j] + db->kmers_32[db->nodes[i].kmers + j + 1];
        }
        if (sum > med) above += 1;
        if (sum < med) below += 1;
      }
    }

    equal = total - above - below;
    if (debug > 1) fprintf (stderr, "Trying median %u (%u) - equal %u, below %u, above %u\n", med / 6, med, equal, below, above);
    /* Special case: min == med, max == med + 1 */
    if (max == (min + 1)) {
      if (above > (below + equal)) {
        /* Max is true median */
        med = max;
      }
      break;
    }
    if (above > below) {
      if ((above - below) < equal) break;
      min = med;
    } else if (below > above) {
      if ((below - above) < equal) break;
      max = med;
    } else {
      break;
    }
    med = (min + max) / 2;
  }
  return med;
}

TaskTable *
task_table_new ()
{
  TaskTable *tt = (TaskTable *) malloc (sizeof (TaskTable));
  memset (tt, 0, sizeof (TaskTable));
  tt->words = (unsigned long long *) malloc (BLOCK_SIZE * 8);
  tt->alleles = (unsigned int *) malloc (BLOCK_SIZE * 4);
  return tt;
}

void
task_table_free (TaskTable *tt)
{
  free (tt->words);
  free (tt->alleles);
  free (tt);
}
