#define __GMER_COUNTER_C__

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <sys/mman.h>

#include "listmaker-queue.h"
#include "utils.h"
#include "trie.h"
#include "sequence.h"
#include "fasta.h"
#include "listmaker-queue.h"
#include "word-map.h"
#include "database.h"
#include "sequence-zstream.h"
#include "version.h"

#define MAX_LINES 10000000000
#define MAX_FILESIZE 10000000000
#define MAX_NUM_THREADS 256

#define BLOCK_SIZE (1024 * 1024 * 10)
#define DEFAULT_NUM_THREADS 24
#define DEFAULT_NUM_TABLES 24

unsigned int debug = 0;

/* File parsing tasks */

typedef struct _TaskRead TaskFile;
typedef struct _SNPQueue SNPQueue;
typedef struct _SNPTable SNPTable;

struct _SNPTable {
  unsigned int nwords;
  unsigned long long *words;
  unsigned int *alleles;
  /* Stats */
  unsigned long long n_nucl;
  unsigned long long n_gc;
  /* Index */
  Read *reads;
};

#define TASK_TABLE (TASK_READ + 1)

typedef struct _TaskTable TaskTable;

struct _TaskTable {
  GT4Task task;
  TaskTable *next;
  unsigned int idx;
  SNPTable *tbl;
};

TaskTable *task_table_new (GT4Queue *queue, unsigned int compile_index);
void task_table_delete (TaskTable *tt);
static SNPTable *snp_table_new (unsigned int compile_index);
static void snp_table_free (SNPTable *tbl);

struct _SNPQueue {
  GT4ListMakerQueue lmq;
  
  unsigned int n_free_tables;
  SNPTable **free_tables;
  unsigned int n_full_tables;
  SNPTable **full_tables;
  /* Stats */
  unsigned long long n_nucl;
  unsigned long long n_gc;
  unsigned long long n_kmers_total;
  unsigned long long n_kmers;
  unsigned long long n_kmer_gc;
  /* Data */
  KMerDB *db;
  /* Read lists */
  ReadList **reads;
};

/* Main thread loop */
static void write_index (SNPQueue *snpq, KMerDB *db, const char *files[], unsigned int n_files);
static void print_counts (SNPQueue *snpq, KMerDB *db);
static void process (GT4Queue *queue, unsigned int idx, void *arg);
static int start_sequence (GT4FastaReader *reader, void *data);
static int end_sequence (GT4FastaReader *reader, void *data);
static int read_nucleotide (GT4FastaReader *reader, unsigned int nucleotide, void *data);
static int read_word (GT4FastaReader *reader, unsigned long long word, void *data);
static int compare_counts (const void *lhs, const void *rhs);
static unsigned int get_pair_median (KMerDB *db);

static void
print_usage (FILE *ofs) {
  fprintf (ofs, "gmer_counter version %u.%u.%u (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_MICRO, VERSION_QUALIFIER);
  fprintf (ofs, "Usage:\n");
  fprintf (ofs, "  gmer_counter ARGUMENTS SEQUENCES...\n");
  fprintf (ofs, "Arguments:\n");
  fprintf (ofs, "    -v | --version   - Print version information and exit\n");
  fprintf (ofs, "    -db DATABASE     - SNP/KMER database file\n");
  fprintf (ofs, "    -dbb DBBINARY    - binary database file\n");
  fprintf (ofs, "    -w FILENAME      - write binary database to file\n");
  fprintf (ofs, "    -32              - use 32-bit integeres for counts (default 16-bit)\n");
  fprintf (ofs, "    --max_kmers NUM  - maximum number of kmers per node\n");
  fprintf (ofs, "    --silent         - do not output kmer counts (useful if only compiling db or index is needed\n");
  fprintf (ofs, "    --header         - print header row\n");
  fprintf (ofs, "    --total          - print the total number of kmers per node\n");
  fprintf (ofs, "    --unique         - print the number of nonzero kmers per node\n");
  fprintf (ofs, "    --kmers          - print individual kmer counts (default if no other output)\n");
  fprintf (ofs, "    --compile_index FILENAME - Add read index to database and write it to file\n");
  fprintf (ofs, "    --distribution NUM  - print kmer distribution (up to given number)\n");
  fprintf (ofs, "    --num_threads    - number of worker threads (default %u)\n", DEFAULT_NUM_THREADS);
  fprintf (ofs, "    --prefetch       - prefetch memory mapped files (faster on high-memory systems)\n");
  fprintf (ofs, "    --recover        - recover from FastA/FastQ errors (useful for corrupted streams)\n");
  fprintf (ofs, "    --stats          - print some statistics about sequence and kmers\n");
  fprintf (ofs, "    -D               - increase debug level\n");
  fprintf (ofs, "    -DDB             - increase database debug level\n");
}

static unsigned int recover = 0;
static unsigned int stats = 0;
unsigned int header = 0, total = 0, unique = 0, kmers = 0, distro = 0;
const char *index_name = NULL;
unsigned int dump_index = 0;

int
main (int argc, const char *argv[])
{
  const char *db_name = NULL;
  const char *dbb = NULL;
  const char *wdb = NULL;
  unsigned int max_kmers_per_node = 1000000000;
  unsigned int silent = 0, big = 0, dm = 0;
  unsigned int lowmem = 1;
  unsigned int nseqs = 0;
  const char *seqnames[1024];
  unsigned long long i;
  double start_time, last_time;

  unsigned int nthreads = DEFAULT_NUM_THREADS;
  SNPQueue snpq;

  KMerDB db;

  
  for (i = 1; i < argc; i++) {
    if (!strcmp (argv[i], "-v") || !strcmp (argv[i], "--version")) {
      fprintf (stdout, "gmer_counter version %u.%u.%u (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_MICRO, VERSION_QUALIFIER);
      exit (0);
    } else if (!strcmp (argv[i], "-h") || !strcmp (argv[i], "--help")) {
      print_usage (stdout);
      exit (0);
    } else if (!strcmp (argv[i], "-db")) {
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
    } else if (!strcmp (argv[i], "--silent")) {
      silent = 1;
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
    } else if (!strcmp (argv[i], "--compile_index")) {
      /* Write index */
      i += 1;
      if (i >= argc) {
        print_usage (stderr);
        exit (1);
      }
      index_name = argv[i];
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
    } else if (!strcmp (argv[i], "--prefetch")) {
      lowmem = 0;
    } else if (!strcmp (argv[i], "--recover")) {
      recover = 1;
    } else if (!strcmp (argv[i], "--stats")) {
      stats = 1;
    } else if (!strcmp (argv[i], "--count_trie_allocations")) {
      gt4_trie_debug |= GT4_TRIE_COUNT_ALLOCATIONS;
    } else if (!strcmp (argv[i], "-D")) {
      /* Debug */
      debug += 1;
    } else if (!strcmp (argv[i], "-DDB")) {
      /* Debug */
      db_debug += 1;
    } else if (!strcmp (argv[i], "--dump_index")) {
      /* Debug */
      dump_index = 1;
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

  start_time = last_time = get_time();

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
    if (debug) fprintf (stderr, "Loading text database %s\n", db_name);
    if (!read_db_from_text (&db, cdata, csize, max_kmers_per_node, (big) ? 32 : 16)) {
      fprintf (stderr, "Cannot read text database %s\n", dbb);
      exit (1);
    }
    if (debug) {
      fprintf (stderr, "Loading time (text): %.1fs\n", get_time() - last_time);
    }
    last_time = get_time();
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
    if (debug) fprintf (stderr, "Finished loading binary database (index = %u)\n", db.index.read_blocks != NULL);
    if (dump_index) {
      gt4_db_dump (&db, stdout);
      exit (0);
    }
    if (debug) {
      fprintf (stderr, "Loading time (binary): %.1fs\n", get_time() - last_time);
    }
    last_time = get_time();
  }

  if (gt4_trie_debug & GT4_TRIE_COUNT_ALLOCATIONS) {
    fprintf (stderr, "Trie: %u allocations, total memory %llu MiB\n", db.trie.num_allocations, db.trie.total_memory / (1024 * 1024));
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
    if (debug) {
      fprintf (stderr, "Writing time (database): %.1fs\n", get_time() - last_time);
    }
    last_time = get_time();
  }

  if (nseqs > 0) {
    /* Set up queue */
    memset (&snpq, 0, sizeof (SNPQueue));
    az_instance_init (&snpq.lmq, GT4_TYPE_LISTMAKER_QUEUE);
    maker_queue_setup (&snpq.lmq, nthreads, db.wordsize, 0, 0, 0);
    snpq.n_free_tables = DEFAULT_NUM_TABLES;
    snpq.free_tables = (SNPTable **) malloc (DEFAULT_NUM_TABLES * sizeof (SNPTable *));
    for (i = 0; i < DEFAULT_NUM_TABLES; i++) {
      snpq.free_tables[i] = snp_table_new (index_name != NULL);
    }
    snpq.full_tables = (SNPTable **) malloc (DEFAULT_NUM_TABLES * sizeof (SNPTable *));
    /* Read files */
    snpq.db = &db;
    for (i = 0; i < nseqs; i++) {
      maker_queue_add_file (&snpq.lmq, seqnames[i], lowmem, i);
    }
    if (index_name) {
      snpq.reads = (ReadList **) malloc (db.n_kmers * sizeof (ReadList *));
      memset (snpq.reads, 0, db.n_kmers * sizeof (ReadList *));
    }
    gt4_queue_create_threads (&snpq.lmq.queue, process, &snpq);
    process (&snpq.lmq.queue, 0, &snpq);
    gt4_queue_lock (&snpq.lmq.queue);
    while (snpq.lmq.queue.nthreads_running > 1) {
      gt4_queue_wait (&snpq.lmq.queue);
    }
    gt4_queue_unlock (&snpq.lmq.queue);

    if (debug) {
      fprintf (stderr, "Reading time: %.1fs\n", get_time() - last_time);
    }
    last_time = get_time();

    fprintf (stdout, "#gmer_counter version %u.%u.%u (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_MICRO, VERSION_QUALIFIER);
    if (db_name) fprintf (stdout, "#TextDatabase\t%s\n", db_name);
    if (dbb) fprintf (stdout, "#BinaryDatabase\t%s\n", dbb);
        
    if (dm) {
      unsigned int med = get_pair_median (&db);
      fprintf (stdout, "#PairMedian\t%u\n", med);
    }

    if (stats) {
      fprintf (stdout, "LENGTH\tGC\tTOTAL_KMERS\tLIST_KMERS\tLIST_KMER_GC\n");
      fprintf (stdout, "%llu", snpq.n_nucl);
      fprintf (stdout, "\t%.3f", (double) snpq.n_gc / snpq.n_nucl);
      fprintf (stdout, "\t%llu", snpq.n_kmers_total);
      fprintf (stdout, "\t%llu", snpq.n_kmers);
      fprintf (stdout, "\t%.3f\n", (double) snpq.n_kmer_gc / (snpq.n_kmers * db.wordsize));
    }

    if (index_name) {
      write_index (&snpq, &db, seqnames, nseqs);
      if (debug) {
        fprintf (stderr, "Index writing time: %.1fs\n", get_time() - last_time);
      }
      last_time = get_time();
    }

    if (!silent) {
      print_counts (&snpq, &db);
    }
    /* Need queue for stats */
    az_instance_finalize (&snpq.lmq, GT4_TYPE_LISTMAKER_QUEUE);
  }
  if (debug) {
    fprintf (stderr, "Total time: %.1fs\n", get_time() - start_time);
  }
  
  return 0;
}

typedef struct _ReadBlock ReadBlock;
struct _ReadBlock {
  GT4Index *index;
  unsigned int start;
  unsigned int end;
  unsigned long long max_name_pos;
  unsigned long long max_kmer_pos;
  unsigned int n_reads;
};

static void
count_reads (GT4Queue *queue, unsigned int idx, void *arg)
{
  ReadBlock *rbs = (ReadBlock *) arg;
  ReadBlock *rb = rbs + idx;
  SNPQueue *snpq = (SNPQueue *) queue;
  unsigned int i;
  if (debug) fprintf (stderr, "counting %u-%u\n", rb->start, rb->end);
  for (i = rb->start; i < rb->end; i++) {
    unsigned int read_idx = 0;
    ReadList *rl;
    for (rl = snpq->reads[i]; rl; rl = rl->next) {
      GT4LMQSource *src = &snpq->lmq.sources[rl->read.source_idx];
      unsigned long long name_pos = src->start + src->subseqs[rl->read.subseq].name_pos;
      if (name_pos > rb->max_name_pos) rb->max_name_pos = name_pos;
      if (rl->read.kmer_pos > rb->max_kmer_pos) rb->max_kmer_pos = rl->read.kmer_pos;
      rb->n_reads += 1;
      read_idx += 1;
    }
    rb->index->read_blocks[i] = (unsigned long long) read_idx;
  }
  gt4_queue_lock (queue);
  gt4_queue_broadcast (queue);
  gt4_queue_unlock (queue);
}

static unsigned long long
write_reads (GT4Index *index, FILE *ofs, void *data)
{
  SNPQueue *snpq = (SNPQueue *) data;
  KMerDB *db = snpq->db;
  unsigned long long i;
  unsigned long long written = 0;
  unsigned long long buf[1024];
  unsigned int bp = 0;
  if (debug) fprintf (stderr, "Writing reads\n");
  for (i = 0; i < db->n_kmers; i++) {
    ReadList *rl;
    /* read_start = db->index.read_blocks[i]; */
    for (rl = snpq->reads[i]; rl; rl = rl->next) {
      GT4LMQSource *src = &snpq->lmq.sources[rl->read.source_idx];
      unsigned long long name_pos = src->start + src->subseqs[rl->read.subseq].name_pos;
      unsigned long long code = ((unsigned long long) rl->read.dir << (db->index.nbits_file + db->index.nbits_npos + db->index.nbits_kmer)) |
        ((unsigned long long) src->file_idx << (db->index.nbits_npos + db->index.nbits_kmer)) |
        (name_pos << db->index.nbits_kmer) |
        rl->read.kmer_pos;
      buf[bp++] = code;
      if (bp >= 1024) {
        fwrite (buf, 8, bp, ofs);
        written += bp;
        bp = 0;
      }
    }
  }
  if (bp) {
    fwrite (buf, 8, bp, ofs);
    written += bp;
  }
  return written;
}

static void
write_index (SNPQueue *snpq, KMerDB *db, const char *files[], unsigned int n_files)
{
  /* Build read index */
  unsigned long long max_name_pos = 0;
  unsigned int max_file_idx = 0, max_kmer_pos = 0;
  unsigned long long read_start;
  unsigned int i;
  double last_time;

  last_time = get_time();

  gt4_db_clear_index (db);

  /* Files */
  db->index.n_files = n_files;
  db->index.files = (char **) malloc (db->index.n_files * sizeof (char *));
  for (i = 0; i < n_files; i++) {
    db->index.files[i] = (char *) files[i];
  }
  db->index.n_kmers = db->n_kmers;
  db->index.read_blocks = (unsigned long long *) malloc (db->n_kmers * sizeof (unsigned long long));

  max_file_idx = n_files - 1;
  read_start = 0;
  if (debug) fprintf (stderr, "Calculate bitsizes\n");
  unsigned int n_threads = snpq->lmq.queue.nthreads_total;
  unsigned int bsize = (db->index.n_kmers + n_threads - 1) / n_threads;
  ReadBlock rbs[MAX_NUM_THREADS];
  memset (rbs, 0, sizeof (rbs));
  unsigned int pos = 0;
  for (i = 0; i < n_threads; i++) {
    rbs[i].index = &db->index;
    rbs[i].start = pos;
    pos += bsize;
    if (pos > db->index.n_kmers) pos = db->index.n_kmers;
    rbs[i].end = pos;
  }
  /* for (i = 0; i < n_threads; i++) {
    fprintf (stderr, "counting %u-%u\n", rbs[i].start, rbs[i].end);
    count_reads (&snpq->lmq.queue, i, rbs);
  }*/
  gt4_queue_create_threads (&snpq->lmq.queue, count_reads, rbs);
  count_reads (&snpq->lmq.queue, 0, rbs);
  gt4_queue_lock (&snpq->lmq.queue);
  while (snpq->lmq.queue.nthreads_running > 1) {
    gt4_queue_wait (&snpq->lmq.queue);
  }
  gt4_queue_unlock (&snpq->lmq.queue);
  for (i = 0; i < n_threads; i++) {
    if (rbs[i].max_name_pos > max_name_pos) max_name_pos = rbs[i].max_name_pos;
    if (rbs[i].max_kmer_pos > max_kmer_pos) max_kmer_pos = rbs[i].max_kmer_pos;
    db->index.n_reads += rbs[i].n_reads;
  }
  unsigned long long kms = 0;
  for (i = 0; i < db->index.n_kmers; i++) {
    unsigned long long nreads = db->index.read_blocks[i];
    db->index.read_blocks[i] = kms;
    kms += nreads;
  }
  if (debug) {
    fprintf (stderr, "Bitsize time: %.1fs\n", get_time() - last_time);
  }
  last_time = get_time();
  
  if (debug) fprintf (stderr, "Num files %u Max name pos %llu Max sequence pos %u\n", n_files, max_name_pos, max_kmer_pos);
  db->index.nbits_file = 1;
  while (max_file_idx > 1) {
    db->index.nbits_file += 1;
    max_file_idx /= 2;
  }
  db->index.nbits_npos = 1;
  while (max_name_pos > 1) {
    db->index.nbits_npos += 1;
    max_name_pos /= 2;
  }
  db->index.nbits_kmer = 1;
  while (max_kmer_pos > 1) {
    db->index.nbits_kmer += 1;
    max_kmer_pos /= 2;
  }
  if (debug) fprintf (stderr, "NBits file %u npos %u kmer %u\n", db->index.nbits_file, db->index.nbits_npos, db->index.nbits_kmer);

  /* Write database with index */
  FILE *ofs;
  if (debug) {
    fprintf (stderr, "Writing index database to %s\n", index_name);
  }
  ofs = (fopen (index_name, "w+"));
  if (!ofs) {
    fprintf (stderr, "Cannot open %s for writing\n", index_name);
    exit (1);
  }
  write_db_to_file_with_reads_callback (db, ofs, 0, write_reads, snpq);
  fclose (ofs);
  if (debug) {
    fprintf (stderr, "Done\n");
  }
  if (debug) {
    fprintf (stderr, "Writing time (reads): %.1fs\n", get_time() - last_time);
  }
  last_time = get_time();
}

static void
print_counts (SNPQueue *snpq, KMerDB *db)
{
  unsigned int i;
  if (header) {
    fprintf (stdout, "NODE\tN_KMERS");
    if (total) fprintf (stdout, "\tTOTAL");
    if (unique) fprintf (stdout, "\tUNIQUE");
    if (kmers) fprintf (stdout, "\tKMERS");
    if (distro) fprintf (stdout, "\tDISTRIBUTION");
    fprintf (stdout, "\n");
  }
  
  for (i = 0; i < db->n_nodes; i++) {
    unsigned int j;
    fprintf (stdout, "%s\t%u", db->names + db->nodes[i].name, db->nodes[i].nkmers);
    if (total) {
      unsigned long long total = 0;
      for (j = 0; j < db->nodes[i].nkmers; j++) {
        if (db->count_bits == 16) {
          total += db->kmers_16[db->nodes[i].kmers + j];
        } else {
          total += db->kmers_32[db->nodes[i].kmers + j];
        }
      }
      fprintf (stdout, "\t%llu", total);
    }
    if (unique) {
      unsigned int uniq = 0;
      for (j = 0; j < db->nodes[i].nkmers; j++) {
        if (db->count_bits == 16) {
          if (db->kmers_16[db->nodes[i].kmers + j]) uniq += 1;
        } else {
          if (db->kmers_16[db->nodes[i].kmers + j]) uniq += 1;
        }
      }
      fprintf (stdout, "\t%u", uniq);
    }
    if (kmers) {
      for (j = 0; j < db->nodes[i].nkmers; j++) {
        if (db->count_bits == 16) {
          fprintf (stdout, "\t%u", db->kmers_16[db->nodes[i].kmers + j]);
        } else {
          fprintf (stdout, "\t%u", db->kmers_32[db->nodes[i].kmers + j]);
        }
      }
    }
    if (distro) {
      static unsigned int c_len = 0;
      static unsigned int *c = NULL;
      if (c_len < db->nodes[i].nkmers) {
        c_len = c_len << 1;
        if (c_len < db->nodes[i].nkmers) c_len = db->nodes[i].nkmers;
        c = (unsigned int *) realloc (c, c_len * 4);
      }
      unsigned int current, count;
      if (db->count_bits == 16) {
        for (j = 0; j < db->nodes[i].nkmers; j++) c[j] = db->kmers_16[db->nodes[i].kmers + j];
      } else {
        memcpy (c, db->kmers_32 + db->nodes[i].kmers, db->nodes[i].nkmers * 4);
      }
      qsort (c, db->nodes[i].nkmers, 4, compare_counts);
      current = 0;
      j = 0;
      while (current <= distro) {
        count = 0;
        while ((j < db->nodes[i].nkmers) && (c[j] == current)) {
          count += 1;
          j += 1;
        }
        fprintf (stdout, "\t%u", count);
        current += 1;
      }
    }
    if (index_name) {
      for (j = 0; j < db->nodes[i].nkmers; j++) {
        unsigned int kmer_idx;
        ReadList *rl;
        kmer_idx = db->nodes[i].kmers + j;
        for (rl = snpq->reads[kmer_idx]; rl; rl = rl->next) {
          fprintf (stdout, " (%u/%u/%u)", rl->read.source_idx, rl->read.subseq, rl->read.kmer_pos);
        }
      }
    }
    fprintf (stdout, "\n");
  }
}

static unsigned int
read_file (SNPQueue *snpq, TaskRead *tr)
{
  unsigned int result;

  /* Read words from file */
  snpq->lmq.n_files_waiting -= 1;
  snpq->lmq.n_files_reading += 1;
  SNPTable *tbl = snpq->free_tables[--snpq->n_free_tables];
  gt4_queue_unlock (&snpq->lmq.queue);
  tbl->nwords = 0;
  tbl->n_nucl = 0;
  tbl->n_gc = 0;
  tr->data = tbl;
  /* if (debug > 0) fprintf (stderr, "Thread %d: reading file %s from %llu\n", idx, tt->_seqfile->path, tf->task_read.reader.cpos); */
  result = fasta_reader_read_nwords (&tr->reader, BLOCK_SIZE, start_sequence, end_sequence, NULL, (stats) ? read_nucleotide : NULL, read_word, tr);
  if (result) {
    fprintf (stderr, "read_file: Fasta reader %s returned %u\n", tr->reader.id, result);
    if (!recover) exit (1);
  }
  /* if (debug > 1) fprintf (stderr, "Thread %d: finished reading %s at %llu\n", idx, tt->_seqfile->path, tf->task_read.reader.cpos); */
  gt4_queue_lock (&snpq->lmq.queue);
  snpq->lmq.n_files_reading -= 1;
  if (result || tr->reader.in_eof) {
    task_read_delete (tr);
  } else {
    gt4_queue_add_task (&snpq->lmq.queue, &tr->task, 0);
    snpq->lmq.n_files_waiting += 1;
  }
  snpq->full_tables[snpq->n_full_tables++] = tbl;
  TaskTable *tt = task_table_new (&snpq->lmq.queue, index_name != NULL);
  tt->tbl = snpq->full_tables[--snpq->n_full_tables];
  gt4_queue_add_task (&snpq->lmq.queue, &tt->task, 0);
  return 0;
}

static unsigned int
process_table (SNPQueue *snpq, TaskTable *tt, unsigned int thread_idx)
{
  KMerDB *db = snpq->db;
  unsigned int i;

  gt4_queue_unlock (&snpq->lmq.queue);
  if (debug > 1) fprintf (stderr, "Thread %d: table lookup\n", thread_idx);
  SNPTable *tbl = tt->tbl;
  for (i = 0; i < tbl->nwords; i++) {
    if ((debug > 1) && ((i % 10000) == 0)) fprintf (stderr, ".");
    tbl->alleles[i] = trie_lookup (&snpq->db->trie, tbl->words[i]);
  }
  if (debug > 1) fprintf (stderr, "Thread %d: finished lookup\n", thread_idx);

  gt4_queue_lock (&snpq->lmq.queue);
  /* fixme: Create separate task / mutex */
  if (stats) {
    snpq->n_nucl += tbl->n_nucl;
    snpq->n_gc += tbl->n_gc;
    snpq->n_kmers_total += tbl->nwords;
  }
  for (i = 0; i <  tbl->nwords; i++) {
    unsigned int code, node, kmer, kmer_idx;
    code = tbl->alleles[i];
    if (!code) continue;
    code &= 0x7fffffff;
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
    kmer_idx = db->nodes[node].kmers + kmer;
    if (db->count_bits == 16) {
      if (db->kmers_16[kmer_idx] < 65535) db->kmers_16[kmer_idx] += 1;
    } else {
      if (db->kmers_32[kmer_idx] < 0xffffffff) db->kmers_32[kmer_idx] += 1;
    }
    if (stats) {
      unsigned int j;
      snpq->n_kmers += 1;
      for (j = 0; j < db->wordsize; j++) {
        unsigned long long word = tbl->words[i];
        snpq->n_kmer_gc += ((word ^ (word >> 1)) & 1);
        word = word >> 2;
      }
    }
    if (snpq->reads) {
      ReadList *rl = gm4_read_list_new ();
      rl->read = tbl->reads[i];
      rl->next = snpq->reads[kmer_idx];
      snpq->reads[kmer_idx] = rl;
    }
  }
  snpq->free_tables[snpq->n_free_tables++] = tbl;
  task_table_delete (tt);
  return 0;
}

static void
process (GT4Queue *queue, unsigned int idx, void *arg)
{
  SNPQueue *snpq = (SNPQueue *) arg;
  unsigned int finished = 0;

  if (debug > 1) {
    gt4_queue_lock (queue);
    fprintf (stderr, "Thread %d started (total %d)\n", idx, queue->nthreads_running);
    gt4_queue_unlock (queue);
  }
  /* Do work */
  while (!finished) {
    unsigned int wait = 0;
    /* Get exclusive lock on queue */
    gt4_queue_lock (queue);
    if (!queue->tasks) {
      /* No tasks left */
      if (!snpq->lmq.n_running) {
        /* No other threads running - finish */
        finished = 1;
      } else {
        /* Other tasks still working - Wait */
        wait = 1;
      }
    } else {
      GT4Task *task;
      for (task = queue->tasks; task; task = task->next) {
        if (task->type == TASK_READ) {
          if (snpq->n_free_tables) {
            gt4_queue_remove_task (&snpq->lmq.queue, task, 0);
            snpq->lmq.n_running += 1;
            wait = read_file (snpq, (TaskRead *) task);
            snpq->lmq.n_running -= 1;
            break;
          }
        } else if (task->type == TASK_TABLE) {
          gt4_queue_remove_task (&snpq->lmq.queue, task, 0);
          snpq->lmq.n_running += 1;
          wait = process_table (snpq, (TaskTable *) task, idx);
          snpq->lmq.n_running -= 1;
          break;
        }
      }
      if (!task) wait = 1;
    }
    if (wait) {
      gt4_queue_wait (&snpq->lmq.queue);
    } else {
      gt4_queue_broadcast (&snpq->lmq.queue);
    }
    gt4_queue_unlock (&snpq->lmq.queue);
  }
  /* Exit if everything is done */
  if (debug > 1) fprintf (stderr, "Thread %u exiting (remaining %d)\n", idx, snpq->lmq.queue.nthreads_running);
}


static int
start_sequence (GT4FastaReader *reader, void *data)
{
  TaskRead *tt = (TaskRead *) data;
  GT4ListMakerQueue *mq = (GT4ListMakerQueue *) tt->task.queue;
  GT4LMQSource *src = &mq->sources[tt->idx];
  maker_queue_add_subsequence (mq, tt->idx, reader->name_pos, reader->name_length);
  src->subseqs[src->n_subseqs - 1].sequence_pos = reader->cpos;
  return 0;
}
            
static int
end_sequence (GT4FastaReader *reader, void *data)
{
  TaskRead *tt = (TaskRead *) data;
  GT4ListMakerQueue *mq = (GT4ListMakerQueue *) tt->task.queue;
  GT4LMQSource *src = &mq->sources[tt->idx];
  src->subseqs[src->n_subseqs - 1].sequence_len = reader->cpos - src->subseqs[src->n_subseqs - 1].sequence_pos;
  return 0;
}


static int
read_word (GT4FastaReader *reader, unsigned long long word, void *data)
{
  TaskRead *tt = (TaskRead *) data;
  GT4ListMakerQueue *mq = (GT4ListMakerQueue *) tt->task.queue;
  GT4LMQSource *src = &mq->sources[tt->idx];
  SNPTable *tbl = (SNPTable *) tt->data;

  tbl->words[tbl->nwords] = word;
  if (tbl->reads) {
    tbl->reads[tbl->nwords].source_idx = tt->idx;
    tbl->reads[tbl->nwords].subseq = src->n_subseqs - 1;
    tbl->reads[tbl->nwords].kmer_pos = reader->seq_npos + 1 - reader->wordlength;
    /* tt->reads[tt->nwords].kmer_pos = reader->cpos - reader->name_pos + 1 - reader->wordlength; */
    tbl->reads[tbl->nwords].dir = (word != reader->wordfw);
  }
  tbl->nwords += 1;
  return 0;
}

static int
read_nucleotide (GT4FastaReader *reader, unsigned int nucl, void *data)
{
  TaskRead *tt = (TaskRead *) data;
  GT4ListMakerQueue *mq = (GT4ListMakerQueue *) tt->task.queue;
  GT4LMQSource *src = &mq->sources[tt->idx];
  SNPTable *tbl = (SNPTable *) tt->data;

  tbl->n_nucl += 1;
  tbl->n_gc += ((nucl ^ (nucl >> 1)) & 1);
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
task_table_new (GT4Queue *queue, unsigned int compile_index)
{
  TaskTable *tt = (TaskTable *) malloc (sizeof (TaskTable));
  memset (tt, 0, sizeof (TaskTable));
  tt->task.queue = queue;
  tt->task.type = TASK_TABLE;
  tt->task.priority = 10;
  return tt;
}

void
task_table_delete (TaskTable *tt)
{
  free (tt);
}

static SNPTable *
snp_table_new (unsigned int compile_index)
{
  SNPTable *tbl = (SNPTable *) malloc (sizeof (SNPTable));
  memset (tbl, 0, sizeof (SNPTable));
  tbl->words = (unsigned long long *) malloc (BLOCK_SIZE * 8);
  tbl->alleles = (unsigned int *) malloc (BLOCK_SIZE * 4);
  if (compile_index) tbl->reads = (Read *) malloc (BLOCK_SIZE * sizeof (Read));
  return tbl;
}

static void
snp_table_free (SNPTable *tbl)
{
  free (tbl->words);
  free (tbl->alleles);
  if (tbl->reads) free (tbl->reads);
}
