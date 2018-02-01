#define __GMER_COUNTER_C__

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <sys/mman.h>

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

#define BLOCK_SIZE (1024 * 1024 * 10)
#define DEFAULT_NUM_THREADS 24
#define DEFAULT_NUM_TABLES 24

#define MAX_KMERS 25

unsigned int debug = 0;

typedef struct _TaskTable TaskTable;
struct _TaskTable {
  TaskTable *next;
  GT4SequenceFile *seqfile;
  unsigned int nwords;
  unsigned long long *words;
  unsigned int *alleles;
  /* Stats */
  unsigned long long n_nucl;
  unsigned long long n_gc;
  /* Read indexing */
  unsigned int file_idx;
  unsigned long long name_npos;
  Read *reads;
};

TaskTable *task_table_new (unsigned int index);
void task_table_free (TaskTable *tt);

typedef struct _SNPQueue SNPQueue;
struct _SNPQueue {
  GT4Queue queue;
  /* Files to process */
  unsigned int nfiles;
  TaskFile *files;
  /* Free tables */
  TaskTable *free_tables;
  /* Filled tables */
  TaskTable *full_tables;
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
static void process (GT4Queue *queue, unsigned int idx, void *arg);
static int start_sequence (GT4FastaReader *reader, void *data);
static int end_sequence (GT4FastaReader *reader, void *data);
static int read_nucleotide (GT4FastaReader *reader, unsigned int nucleotide, void *data);
static int read_word (GT4FastaReader *reader, unsigned long long word, void *data);
static int compare_counts (const void *lhs, const void *rhs);
static unsigned int get_pair_median (KMerDB *db);

static void
print_usage (FILE *ofs) {
  fprintf (ofs, "gmer_counter version %u.%u (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_QUALIFIER);
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

static void
add_task (SNPQueue *snpq, const char *filename)
{
  unsigned int len = strlen (filename);
  if ((len > 3) && !strcmp (filename + len - 3, ".gz")) {
    fprintf (stderr, "Opening compressed stream %s\n", filename);
    GT4SequenceZStream *zstream = gt4_sequence_zstream_new (filename);
    TaskFile *tf = task_file_new_from_source (AZ_OBJECT (zstream), filename, 0);
    tf->next = snpq->files;
    tf->idx = snpq->nfiles++;
    snpq->files = tf;
  } else if (!strcmp (filename, "-")) {
    GT4SequenceStream *stream;
    stream = gt4_sequence_stream_new_from_stream (stdin, 0);
    TaskFile *tf = task_file_new_from_source (AZ_OBJECT (stream), "stdin", 0);
    tf->next = snpq->files;
    tf->idx = snpq->nfiles++;
    snpq->files = tf;
  } else {
    unsigned int nseqs, i;
    GT4SequenceBlock *seqs[32];
    GT4SequenceFile *seqf = gt4_sequence_file_new (filename, 1);
    gt4_sequence_file_map_sequence (seqf);
    nseqs = (unsigned int) (seqf->block.csize / 10000000000ULL) + 1;
    if (nseqs > 32) nseqs = 32;
    nseqs = gt4_sequence_block_split (&seqf->block, seqs, nseqs);
    for (i = 0; i < nseqs; i++) {
      TaskFile *tf = task_file_new_from_source (AZ_OBJECT (seqs[i]), "block", 0);
      if (debug) fprintf (stderr, "%s:%u from %llu to %llu\n", filename, i, (unsigned long long) seqs[i]->cdata - (unsigned long long) seqf->block.cdata, (unsigned long long) seqs[i]->cdata - (unsigned long long) seqf->block.cdata + seqs[i]->csize);
      tf->next = snpq->files;
      tf->idx = snpq->nfiles++;
      snpq->files = tf;
    }
    az_object_unref (AZ_OBJECT (seqf));
  }
}

int
main (int argc, const char *argv[])
{
  const char *db_name = NULL;
  const char *dbb = NULL;
  const char *wdb = NULL;
  const char *index = NULL;
  unsigned int max_kmers_per_node = 1000000000;
  unsigned int silent = 0, header = 0, total = 0, unique = 0, kmers = 0, distro = 0, big = 0, dm = 0;
  unsigned int lowmem = 1;
  unsigned int nseqs = 0;
  const char *seqnames[1024];
  GT4SequenceFile *seq_files[1024];
  unsigned long long i;

  unsigned int nthreads = DEFAULT_NUM_THREADS;
  SNPQueue snpq;

  KMerDB db;

  
  for (i = 1; i < argc; i++) {
    if (!strcmp (argv[i], "-v") || !strcmp (argv[i], "--version")) {
      fprintf (stdout, "gmer_counter version %d.%d (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_QUALIFIER);
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
      index = argv[i];
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
    if (debug) fprintf (stderr, "Loading text database %s\n", db_name);
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
    if (debug) fprintf (stderr, "Finished loading binary database (index = %u)\n", db.index.read_blocks != NULL);
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
  }

  if (nseqs > 0) {
    memset (&snpq, 0, sizeof (SNPQueue));
    az_instance_init (&snpq.queue, GT4_TYPE_QUEUE);
    gt4_queue_setup (&snpq.queue, nthreads);
    /* Read files */
    snpq.db = &db;
    for (i = 0; i < nseqs; i++) {
      TaskFile *tf;
      if (index) {
        if (!strcmp (seqnames[i], "-")) {
          GT4SequenceStream *stream;
          stream = gt4_sequence_stream_new_from_stream (stdin, 0);
          tf = task_file_new_from_source (AZ_OBJECT (stream), seqnames[i], 0);
          az_object_unref (AZ_OBJECT (stream));
        } else {
          tf = task_file_new (seqnames[i], !lowmem);
        }
        seq_files[i] = tf->seqfile;
        gt4_sequence_file_ref (seq_files[i]);
        tf->next = snpq.files;
        tf->idx = i;
        snpq.files = tf;
        snpq.nfiles += 1;
      } else {
        add_task (&snpq, seqnames[i]);
      }
    }
    for (i = 0; i < DEFAULT_NUM_TABLES; i++) {
      TaskTable *tt = task_table_new (index != NULL);
      tt->next = snpq.free_tables;
      snpq.free_tables = tt;
    }
    if (index) {
      snpq.reads = (ReadList **) malloc (db.n_kmers * sizeof (ReadList *));
      memset (snpq.reads, 0, db.n_kmers * sizeof (ReadList *));
    }
    gt4_queue_create_threads (&snpq.queue, process, &snpq);
    process (&snpq.queue, 0, &snpq);
    gt4_queue_lock (&snpq.queue);
    while (snpq.queue.nthreads_running > 1) {
      gt4_queue_wait (&snpq.queue);
    }
    gt4_queue_unlock (&snpq.queue);

    if (debug) {
      fprintf (stderr, "Finished reading files\n");
    }

    fprintf (stdout, "#gmer_counter version %u.%u (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_QUALIFIER);
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

    /* Need queue for stats */
    az_instance_finalize (&snpq.queue, GT4_TYPE_QUEUE);

    if (header) {
      fprintf (stdout, "NODE\tN_KMERS");
      if (total) fprintf (stdout, "\tTOTAL");
      if (unique) fprintf (stdout, "\tUNIQUE");
      if (kmers) fprintf (stdout, "\tKMERS");
      if (distro) fprintf (stdout, "\tDISTRIBUTION");
      fprintf (stdout, "\n");
    }
  
    if (index) {
      /* Build read index */
      unsigned long long max_name_pos = 0;
      unsigned int max_file_idx = 0, max_kmer_pos = 0;
      unsigned long long read_start = 0;
      
      gt4_db_clear_index (&db);

      /* Files */
      db.index.n_files = nseqs;
      db.index.files = (char **) malloc (db.index.n_files * sizeof (char *));
      for (i = 0; i < nseqs; i++) {
        db.index.files[i] = (char *) seqnames[i];
      }
      max_file_idx = nseqs - 1;

      if (debug) fprintf (stderr, "Calculate bitsizes\n");
      for (i = 0; i < db.n_nodes; i++) {
        unsigned int j;
        for (j = 0; j < db.nodes[i].nkmers; j++) {
          unsigned int kmer_idx;
          ReadList *rl;
          kmer_idx = db.nodes[i].kmers + j;
          for (rl = snpq.reads[kmer_idx]; rl; rl = rl->next) {
            unsigned long long name_pos = seq_files[rl->read.file_idx]->subseqs[rl->read.subseq].name_pos;
            if (name_pos > max_name_pos) max_name_pos = name_pos;
            if (rl->read.kmer_pos > max_kmer_pos) max_kmer_pos = rl->read.kmer_pos;
          }
        }
      }
      
      if (debug) fprintf (stderr, "Num files %u Max name pos %llu Max sequence pos %u\n", nseqs, max_name_pos, max_kmer_pos);
      db.index.nbits_file = 1;
      while (max_file_idx > 1) {
        db.index.nbits_file += 1;
        max_file_idx /= 2;
      }
      db.index.nbits_npos = 1;
      while (max_name_pos > 1) {
        db.index.nbits_npos += 1;
        max_name_pos /= 2;
      }
      db.index.nbits_kmer = 1;
      while (max_kmer_pos > 1) {
        db.index.nbits_kmer += 1;
        max_kmer_pos /= 2;
      }
      if (debug) fprintf (stderr, "NBits file %u npos %u kmer %u\n", db.index.nbits_file, db.index.nbits_npos, db.index.nbits_kmer);
      db.index.n_kmers = db.n_kmers;
      db.index.read_blocks = (unsigned long long *) malloc (db.n_kmers * sizeof (unsigned long long));
      if (debug) fprintf (stderr, "Calculate number of reads\n");
      for (i = 0; i < db.n_kmers; i++) {
        ReadList *rl;
        for (rl = snpq.reads[i]; rl; rl = rl->next) db.index.n_reads += 1;
      }
      db.index.reads = (unsigned long long *) malloc (db.index.n_reads * sizeof (unsigned long long));
      if (debug) fprintf (stderr, "Writing reads\n");
      for (i = 0; i < db.n_kmers; i++) {
        unsigned int read_idx = 0;
        ReadList *rl;
        for (rl = snpq.reads[i]; rl; rl = rl->next) {
          unsigned long long name_pos = seq_files[rl->read.file_idx]->subseqs[rl->read.subseq].name_pos;
          unsigned long long code = ((unsigned long long) rl->read.dir << (db.index.nbits_file + db.index.nbits_npos + db.index.nbits_kmer)) |
            ((unsigned long long) rl->read.file_idx << (db.index.nbits_npos + db.index.nbits_kmer)) |
            (name_pos << db.index.nbits_kmer) |
            rl->read.kmer_pos;
          db.index.reads[read_start + read_idx] = code;
#if 0
          unsigned long long _name_pos = (code >> db.index.nbits_kmer) & ((1ULL << db.index.nbits_npos) - 1);
          unsigned long long _file_idx = (code >> (db.index.nbits_npos + db.index.nbits_kmer)) & ((1ULL << db.index.nbits_file) - 1);
          unsigned long long _dir = (code >> (db.index.nbits_npos + db.index.nbits_kmer + db.index.nbits_file)) & 1;
          unsigned long long _kmer_pos = code & ((1ULL << db.index.nbits_kmer) - 1);
          assert (_name_pos == name_pos);
          assert (_file_idx == rl->read.file_idx);
          assert (_dir == rl->read.dir);
          assert (_kmer_pos == rl->read.kmer_pos);
#endif
          read_idx += 1;
        }
        db.index.read_blocks[i] = (unsigned long long) read_start << 24 | read_idx;
        read_start += read_idx;
      }
      /* Write database with index */
      FILE *ofs;
      if (debug) {
        fprintf (stderr, "Writing index database to %s\n", index);
      }
      ofs = (fopen (index, "w+"));
      if (!ofs) {
        fprintf (stderr, "Cannot open %s for writing\n", index);
        exit (1);
      }
      write_db_to_file (&db, ofs, 0);
      fclose (ofs);
      if (debug) {
        fprintf (stderr, "Done\n");
      }
    }

    if (!silent) {
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
        if (index) {
          for (j = 0; j < db.nodes[i].nkmers; j++) {
            unsigned int kmer_idx;
            ReadList *rl;
            kmer_idx = db.nodes[i].kmers + j;
            for (rl = snpq.reads[kmer_idx]; rl; rl = rl->next) {
              fprintf (stdout, " (%u/%u/%u)", rl->read.file_idx, rl->read.subseq, rl->read.kmer_pos);
            }
          }
        }
        fprintf (stdout, "\n");
      }
    }
  }
  
  return 0;
}

static void
process (GT4Queue *queue, unsigned int idx, void *arg)
{
  SNPQueue *snpq;
  KMerDB *db;
  unsigned int finished;

  snpq = (SNPQueue *) arg;
  db = snpq->db;

  finished = 0;
        
  if (debug > 1) fprintf (stderr, "Thread %d started (total %d)\n", idx, snpq->queue.nthreads_running);
  /* Do work */
  while (!finished) {
    /* Get exclusive lock on queue */
    gt4_queue_lock (queue);
    if ((snpq->files) && (snpq->free_tables)) {
      TaskFile *tf;
      TaskTable *tt;
      unsigned int result;
      /* Create new file reading task */
      tf = snpq->files;
      snpq->files = tf->next;
      tt = snpq->free_tables;
      snpq->free_tables = tt->next;
      pthread_mutex_unlock (&snpq->queue.mutex);
      /* Read words from file */
      if (tt->seqfile) gt4_sequence_file_unref (tt->seqfile);
      tt->seqfile = tf->seqfile;
      gt4_sequence_file_ref (tt->seqfile);
      tt->nwords = 0;
      tt->file_idx = tf->idx;
      tt->n_nucl = 0;
      tt->n_gc = 0;
      if (debug > 0) fprintf (stderr, "Thread %d: reading file %s from %llu\n", idx, tt->seqfile->path, tf->reader.cpos);
      result = task_file_read_nwords (tf, BLOCK_SIZE, snpq->db->wordsize, start_sequence, end_sequence, NULL, (stats) ? read_nucleotide : NULL, read_word, tt);
      if (result) {
        fprintf (stderr, "Error reading from file %s\n", tt->seqfile->path);
        if (!recover) exit (1);
      }
      if (debug > 1) fprintf (stderr, "Thread %d: finished reading %s at %llu\n", idx, tt->seqfile->path, tf->reader.cpos);
      pthread_mutex_lock (&snpq->queue.mutex);
      if (result || tf->reader.in_eof) {
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
      /* Create new lookup task */
      tt = (TaskTable *) snpq->full_tables;
      snpq->full_tables = tt->next;

      pthread_mutex_unlock (&snpq->queue.mutex);
      if (debug > 1) fprintf (stderr, "Thread %d: table lookup\n", idx);
      for (i = 0; i < tt->nwords; i++) {
        if ((debug > 1) && ((i % 10000) == 0)) fprintf (stderr, ".");
        tt->alleles[i] = trie_lookup (&snpq->db->trie, tt->words[i]);
      }
      if (debug > 1) fprintf (stderr, "Thread %d: finished lookup\n", idx);

      pthread_mutex_lock (&snpq->queue.mutex);
      /* fixme: Create separate task / mutex */
      if (stats) {
        snpq->n_nucl += tt->n_nucl;
        snpq->n_gc += tt->n_gc;
        snpq->n_kmers_total += tt->nwords;
      }
      for (i = 0; i <  tt->nwords; i++) {
        unsigned int code, node, kmer, kmer_idx;
        code = tt->alleles[i];
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
            unsigned long long word = tt->words[i];
            snpq->n_kmer_gc += ((word ^ (word >> 1)) & 1);
            word = word >> 2;
          }
        }
        if (snpq->reads) {
          ReadList *rl = gm4_read_list_new ();
          rl->read = tt->reads[i];
          rl->next = snpq->reads[kmer_idx];
          snpq->reads[kmer_idx] = rl;
        }
      }
      gt4_sequence_file_unref (tt->seqfile);
      tt->seqfile = NULL;
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
  if (debug > 1) fprintf (stderr, "Thread %u exiting (remaining %d)\n", idx, snpq->queue.nthreads_running);
}


static int
start_sequence (GT4FastaReader *reader, void *data)
{
  TaskTable *tt = (TaskTable *) data;
  if (debug > 2) fprintf (stderr, "%s\n", reader->name);
  gt4_sequence_file_lock (tt->seqfile);
  gt4_sequence_file_add_subsequence (tt->seqfile, reader->name_pos, reader->name_length);
  tt->seqfile->subseqs[tt->seqfile->n_subseqs - 1].sequence_pos = reader->cpos;
  gt4_sequence_file_unlock (tt->seqfile);
  return 0;
}

static int
end_sequence (GT4FastaReader *reader, void *data)
{
  TaskTable *tt = (TaskTable *) data;
  if (debug > 2) fprintf (stderr, "%s\n", reader->name);
  gt4_sequence_file_lock (tt->seqfile);
  tt->seqfile->subseqs[tt->seqfile->n_subseqs - 1].sequence_len = reader->cpos - tt->seqfile->subseqs[tt->seqfile->n_subseqs - 1].sequence_pos;
  gt4_sequence_file_unlock (tt->seqfile);
  return 0;
}

static int
read_word (GT4FastaReader *reader, unsigned long long word, void *data)
{
  TaskTable *tt = (TaskTable *) data;
  tt->words[tt->nwords] = word;
  if (tt->reads) {
    tt->reads[tt->nwords].file_idx = tt->file_idx;
    tt->reads[tt->nwords].subseq = tt->seqfile->n_subseqs - 1;
    tt->reads[tt->nwords].kmer_pos = reader->seq_npos + 1 - reader->wordlength;
    tt->reads[tt->nwords].dir = (word != reader->wordfw);
  }
  tt->nwords += 1;
  return 0;
}

static int
read_nucleotide (GT4FastaReader *reader, unsigned int nucl, void *data)
{
  TaskTable *tt = (TaskTable *) data;
  tt->n_nucl += 1;
  tt->n_gc += ((nucl ^ (nucl >> 1)) & 1);
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
task_table_new (unsigned int index)
{
  TaskTable *tt = (TaskTable *) malloc (sizeof (TaskTable));
  memset (tt, 0, sizeof (TaskTable));
  tt->words = (unsigned long long *) malloc (BLOCK_SIZE * 8);
  tt->alleles = (unsigned int *) malloc (BLOCK_SIZE * 4);
  if (index) {
    tt->reads = (Read *) malloc (BLOCK_SIZE * sizeof (Read));
  }
  return tt;
}

void
task_table_free (TaskTable *tt)
{
  if (tt->seqfile) gt4_sequence_file_unref (tt->seqfile);
  free (tt->words);
  free (tt->alleles);
  if (tt->reads) {
    free (tt->reads);
  }
  free (tt);
}
