#define __GASSEMBLER_C__

#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>

#include <binomial.h>
#include <database.h>
#include <math.h>
#include <matrix.h>
#include <queue.h>
#include <sequence.h>
#include <utils.h>
#include <version.h>

#define OUTPUT_FULL_STATS
#define noUSE_SUB

unsigned int debug = 0;
unsigned int debug_groups = 0;
unsigned int flush = 1;
static unsigned int alternative_calls = 0;

unsigned int max_regions = 1000000000;
unsigned int n_threads = 24;
static unsigned int min_coverage = 4;
static double min_p = 0.95;
static double min_pmut = 0.5;
#define SEX_AUTO 0
#define SEX_MALE 1
#define SEX_FEMALE 2
static unsigned int sex = 0;
#define OUTPUT_POLY_BEST 0
#define OUTPUT_ALL_BEST 1
#define OUTPUT_ALL 2
static unsigned int output = OUTPUT_POLY_BEST;
static unsigned int print_extra = 0;
static double error_prob = 0.001f;
static unsigned int exome = 0;
#define COVERAGE_IGNORE -2
#define COVERAGE_LOCAL -1
#define COVERAGE_MEDIAN 0
/*
 * -2 - per position
 * -1 - per block
 *  0  - median
 * >0  - value
 */

static float coverage = 0;

#define WORDLEN 25
#define MAX_THREADS 256
#define MAX_KMERS 1024
#define MAX_READS_PER_KMER 200
#define MAX_READS 4096
#define MIN_READS 10
#define MAX_ALIGNED_READS 1024
#define MAX_READ_LENGTH max_read_length
#define MAX_REFERENCE_LENGTH max_reference_length
#define MAX_GROUPS MAX_ALIGNED_READS
static unsigned int max_read_length = 200;
static unsigned int max_reference_length = 200;

/* For matrix indexing */
#define A_COLS (MAX_REFERENCE_LENGTH * 2)

typedef struct _GASMQueue GASMQueue;
typedef struct _GASMRead GASMRead;
typedef struct _ReadInfo ReadInfo;
typedef struct _SeqFile SeqFile;
typedef struct _SNV SNV;
typedef struct _SWCell SWCell;
typedef struct _AssemblyData AssemblyData;
typedef struct _CallBlock CallBlock;
typedef struct _Call Call;
typedef struct _CallExtra CallExtra;
typedef struct _Group Group;

/* Returns number of aligned positions */
static int align (AssemblyData *adata, const char *kmers[], unsigned int nkmers);
static int group (AssemblyData *adata, unsigned int print);
static void recalculate_and_call (AssemblyData *adata, Group *groups, unsigned int n_groups, unsigned int *good_groups, unsigned int n_included, unsigned int haploid, unsigned int print);
/* Returns true if last call is written */
static unsigned int call (AssemblyData *adata, CallBlock *cb, unsigned int a_pos, unsigned int sub, CallExtra *extra, unsigned int *call_alignment, unsigned int homozygote);
static int assemble (AssemblyData *adata, const char *kmers[], unsigned int nkmers, unsigned int print);
static int assemble_recursive (GT4GmerDB *db, SeqFile *files, unsigned int ref_chr, unsigned int ref_start, unsigned int ref_end, const char *ref, const char *kmers[], unsigned int nkmers);
unsigned int align_reads_to_reference (NSeq *ref_seq, GASMRead *reads[], unsigned int nreads, GASMRead *a_reads[], short *a, SWCell *sw_matrix);
unsigned int create_gapped_alignment (NSeq *ref_seq, unsigned int ref_start, GASMRead *a_reads[], unsigned int na, short *a, unsigned int aligned_ref[], int ref_pos[], short *ga);
static void test_alignment (const char *a, const char *b);
static float find_coverage (GT4Index *index);
static double calc_p_mdetect (Call *call, CallExtra *extra, unsigned int kmer_cov);
static double calc_p_qual_haploid (Call *call, CallExtra *extra, unsigned int kmer_cov);
static double calc_p_qual_diploid (Call *call, CallExtra *extra, unsigned int kmer_cov);
static double calc_p_select_haploid (Call *call, CallExtra *extra, unsigned int kmer_cov);
static double calc_p_select_diploid (Call *call, CallExtra *extra, unsigned int kmer_cov, unsigned int n0, unsigned int n1);

struct _GASMRead {
  char *name;
  char *seq;
  NSeq *nseq;
  unsigned long long tag;
  unsigned long long mask;
  unsigned long long unknown;
  unsigned short group;
  unsigned short dir;
};

struct _ReadInfo {
  unsigned long long name_pos;
  unsigned int kmer_pos;
  unsigned int file_idx;
  unsigned int dir;
};

static GASMRead *gasm_read_new (const char *name, const char *seq, unsigned int wlen, unsigned int dir);
static void gasm_read_delete (GASMRead *read);

struct _SeqFile {
  const char *name;
  const unsigned char *cdata;
  unsigned long long csize;
};

struct _SNV {
  unsigned int chr;
  unsigned long long pos;
  char *id;
  unsigned short ref_allele;
  unsigned short alt_allele;
  unsigned short genotype;
};

struct _SWCell {
  int16_t score;
  int16_t left_gap_score;
  int16_t top_gap_score;
  int8_t sx;
  int8_t sy;
  int8_t left_gap_len;
  int8_t top_gap_len;
};

struct _CallExtra {
  float prob;
  float rprob;
  float hzprob;
  unsigned short end_dist;
  unsigned short n_groups_total;
  unsigned short n_groups;
  unsigned short div_0;
  unsigned short div_1;
  unsigned short max_cov_0;
  unsigned short max_cov_1;
  unsigned short compat_0;
  unsigned short compat_1;
  unsigned short compat_both;
};

struct _Call {
  unsigned int pos;
  unsigned char sub;
  unsigned char ref;
  unsigned short cov;
  unsigned short counts[GAP + 1];
  unsigned short nucl[2];
  unsigned short poly;
  unsigned short prev_ref;
  float p;
  float q;
  float p_det;
  CallExtra extra;
};

struct _AssemblyData {
  GASMQueue *queue;
  GT4GmerDB *db;
  SeqFile *files;
  /* Reference data */
  unsigned int chr;
  unsigned int start;
  unsigned int end;
  const char *ref;
  NSeq *ref_seq;
  /* Extracted reads */
  GASMRead *reads[MAX_READS];
  unsigned int nreads;
  /* Smith-Waterman table */
  SWCell *sw_matrix;
  /* Gapped alignment (nucleotide values) */
  //short (*alignment)[MAX_REFERENCE_LENGTH * 2];
  short *alignment;

  GASMRead *aligned_reads[MAX_ALIGNED_READS];
  unsigned int *aligned_ref;
  int *ref_pos;
  unsigned int na, p_len;

  /* Number of reads per position */
  short *coverage;
  /* Nucleotide counts */
  short (*nucl_counts)[GAP + 1];
  /* Group compatibility */
  unsigned char (*is_compat)[MAX_GROUPS];
  unsigned short (*n_common)[MAX_GROUPS];
  /* Calls */
  CallBlock *cblock;
};

struct _CallBlock {
  CallBlock *next;
  unsigned int chr;
  unsigned int start;
  unsigned int end;
  unsigned int n_calls;
  unsigned int chr_cov;
  unsigned int haploid;
  /* Calls */
  Call *calls;
};

struct _GASMQueue {
  GT4Queue gt4_queue;
  /* Source */
  const unsigned char *cdata;
  unsigned long long csize;
  unsigned long long cpos;
  unsigned int line;
  unsigned int nrunning;
  unsigned int finished;
  AssemblyData *adata[MAX_THREADS];
  CallBlock *free_blocks;
  CallBlock *processing_blocks;
  CallBlock *finished_blocks;
  /* Last printed position */
  unsigned int last_chr;
  unsigned int last_pos;
};

static void
queue_setup (GASMQueue *queue, unsigned int nthreads)
{
  gt4_queue_setup (&queue->gt4_queue, nthreads);
  queue->cpos = 0;
  queue->line = 0;
  queue->nrunning = 0;
  queue->finished = 0;
}

static CallBlock *
queue_get_call_block (GASMQueue *gq, unsigned int chr, unsigned int start, unsigned int end)
{
  CallBlock *cb;
  if (gq->free_blocks) {
    cb = gq->free_blocks;
    gq->free_blocks = cb->next;
  } else {
    cb = (CallBlock *) malloc (sizeof (CallBlock));
    memset (cb, 0, sizeof (CallBlock));
    cb->calls = (Call *) malloc (MAX_REFERENCE_LENGTH * 2 * sizeof (Call));
  }
  memset (cb->calls, 0, MAX_REFERENCE_LENGTH * 2 * sizeof (Call));
  cb->next = gq->processing_blocks;
  cb->chr = chr;
  cb->start = start;
  cb->end = end;
  cb->n_calls = 0;
  cb->chr_cov = 0;
  cb->haploid = (((sex == SEX_MALE) && ((chr == CHR_X) || (chr == CHR_Y))) || (chr == CHR_MT));
  gq->processing_blocks = cb;
  return cb;
}

static void
queue_finish_call_block (GASMQueue *gq, CallBlock *cb)
{
  CallBlock *prev = NULL;
  CallBlock *cur = gq->processing_blocks;
  while (cur != cb) {
    prev = cur;
    cur = cur->next;
  }
  if (!prev) {
    gq->processing_blocks = cb->next;
  } else {
    prev->next = cb->next;
  }
  cb->next = gq->finished_blocks;
  gq->finished_blocks = cb;
}

static void
queue_free_call_block (GASMQueue *gq, CallBlock *cb)
{
  CallBlock *prev = NULL;
  CallBlock *cur = gq->finished_blocks;
  while (cur != cb) {
    prev = cur;
    cur = cur->next;
  }
  if (!prev) {
    gq->finished_blocks = cb->next;
  } else {
    prev->next = cb->next;
  }
  cb->next = gq->free_blocks;
  gq->free_blocks = cb;
}

static AssemblyData *
assembly_data_new (GT4GmerDB *db, SeqFile *files)
{
  AssemblyData *adata = (AssemblyData *) malloc (sizeof (AssemblyData));
  memset (adata, 0, sizeof (AssemblyData));
  adata->db = db;
  adata->files = files;
  adata->sw_matrix = (SWCell *) malloc ((MAX_REFERENCE_LENGTH + 1) * (MAX_READ_LENGTH + 1) * sizeof (SWCell));
  adata->alignment = (short *) malloc (MAX_ALIGNED_READS * MAX_REFERENCE_LENGTH * 2 * 2);
  adata->aligned_ref = (unsigned int *) malloc (MAX_REFERENCE_LENGTH * 2 * 4);
  memset (adata->aligned_ref, 0, MAX_REFERENCE_LENGTH * 2 * 4);
  adata->ref_pos = (int *) malloc (MAX_REFERENCE_LENGTH * 2 * 4);
  memset (adata->ref_pos, 0, MAX_REFERENCE_LENGTH * 2 * 4);
  adata->coverage = (short *) malloc (MAX_REFERENCE_LENGTH * 2 * 2);
  adata->nucl_counts = (short (*)[GAP + 1]) malloc (MAX_REFERENCE_LENGTH * 2 * (GAP + 1) * 2);
  adata->is_compat = (unsigned char (*)[MAX_GROUPS]) malloc (MAX_GROUPS * MAX_GROUPS);
  adata->n_common = (unsigned short (*)[MAX_GROUPS]) malloc (MAX_GROUPS * MAX_GROUPS * 2);
  return adata;
}

static void
assembly_data_clear (AssemblyData *adata)
{
  unsigned int i;
  for (i = 0; i < adata->nreads; i++) gasm_read_delete (adata->reads[i]);
  adata->nreads = 0;
  if (adata->ref_seq) {
    n_seq_delete (adata->ref_seq);
    adata->ref_seq = NULL;
  }
  adata->cblock = NULL;
}

static void
print_header (FILE *ofs) {
  fprintf (ofs, "CHR\tPOS\tSUB\tREF\tCOV\tCALL\tCLASS\tP\tPMUT");
  if (print_extra > 1) fprintf (ofs, "\tPREV");
  if (print_extra > 0) {
    fprintf (ofs, "\tA\tC\tG\tT\tGAP");
  }
  if (print_extra > 1) {
    fprintf (ofs, "\tPROB\tRPROB\tHZPROB\tEDIST\tGRP_ALL\tGRP\tDIV0\tDIV1\tG0\tG1\tG0_COMP\tG1_COMP\tCOMP_2");
  }
}

static void
print_call (CallBlock *cb, unsigned int pos)
{
  Call *call = &cb->calls[pos];
  /* CHR POS REF COV */
  fprintf (stdout, "%s\t%u\t%u\t%c\t%u", chr_names[cb->chr], call->pos, call->sub, n2c[call->ref], call->cov);
  /* CALL */
  if ((call->ref != N) && (call->cov >= min_coverage) && (call->q >= min_p) && (call->poly || (call->p_det >= min_pmut)) && (call->nucl[0] != NONE)) {
    fprintf (stdout, "\t%c%c", n2c[call->nucl[0]], n2c[call->nucl[1]]);
  } else {
    fprintf (stdout, "\tNC");
  }
  /* CLASS */
  if (call->ref == GAP) {
    fprintf (stdout, "\tI");
  } else if (call->nucl[1] == GAP) {
    fprintf (stdout, "\tD");
  } else if (call->poly) {
    fprintf (stdout, "\tS");
  } else {
    fprintf (stdout, "\t0");
  }
  /* PVALUE */
  fprintf (stdout, "\t%.3f", call->q);
  fprintf (stdout, "\t%.3f", call->p_det);
  if (print_extra > 1) fprintf (stdout, "\t%c", call->prev_ref);
  if (print_extra > 0) {
    /* A C G T GAP */
    fprintf (stdout, "\t%u\t%u\t%u\t%u\t%u", call->counts[A], call->counts[C], call->counts[G], call->counts[T], call->counts[GAP]);
  }
#ifdef OUTPUT_FULL_STATS
  if (print_extra > 1) {
    fprintf (stdout, "\t%.5f\t%.5f\t%.5f", call->extra.prob, call->extra.rprob, call->extra.hzprob);
    fprintf (stdout, "\t%2u", call->extra.end_dist);
    fprintf (stdout, "\t%2u\t%2u\t%2u\t%2u", call->extra.n_groups_total, call->extra.n_groups, call->extra.div_0, call->extra.div_1);
    fprintf (stdout, "\t%2u\t%2u\t%2u\t%2u\t%2u", call->extra.max_cov_0, call->extra.max_cov_1, call->extra.compat_0, call->extra.compat_1, call->extra.compat_both);
  }
#endif
  if (flush) fflush (stdout);
}

static void
print_calls_poly_best (CallBlock *cb_f, GASMQueue *queue, unsigned int only_poly)
{
  unsigned int pos;
    /* Print all relevant positions */
    for (pos = cb_f->start; pos < cb_f->end; pos++) {
      CallBlock *best_cb, *ccb;
      float best_p;
      unsigned int j, has_poly;
      if ((cb_f->chr == queue->last_chr) && (pos <= queue->last_pos)) continue;
      /* Find call for given position with best aggregate value */
      best_cb = cb_f;
      best_p = 0;
      has_poly = 0;
      for (ccb = queue->finished_blocks; ccb; ccb = ccb->next) {
        unsigned int local_poly = 0;
        if (ccb->chr > cb_f->chr) continue;
        if (ccb->start > pos) continue;
        for (j = 0; j < ccb->n_calls; j++) {
          if (ccb->calls[j].pos > pos) break;
          if (ccb->calls[j].pos != pos) continue;
          /* We have the same position (there can be many in ccb) */
          if (ccb->calls[j].poly) local_poly = 1;
          if (ccb->calls[j].p < best_p) continue;
          best_cb = ccb;
          best_p = ccb->calls[j].p;
        }
        if (best_cb == ccb) has_poly = local_poly;
      }
      if (only_poly) {
        if (has_poly) {
          /* Iterate over best call block & print all calls for this position */
          for (j = 0; j < best_cb->n_calls; j++) {
            double p;
            if (best_cb->calls[j].pos > pos) break;
            if (best_cb->calls[j].pos != pos) continue;
            p = best_cb->calls[j].q;
            if (p >= min_p) {
              /* p >= cutoff, print all changed positions (changed nucleotides) */
              if (best_cb->calls[j].poly) {
                print_call (best_cb, j);
                fprintf (stdout, "\n");
              }
            } else {
              /* p < cutoff, print single position (will be NC) */
              print_call (best_cb, j);
              fprintf (stdout, "\n");
              break;
            }
          }
        } else {
          for (j = 0; j < best_cb->n_calls; j++) {
            double pmut;
            if (best_cb->calls[j].pos > pos) break;
            if (best_cb->calls[j].pos != pos) continue;
            /* Best is reference */
            pmut = best_cb->calls[j].p_det;
            if (pmut < min_pmut) {
              print_call (best_cb, j);
              fprintf (stdout, "\n");
            }
          }
        }
      } else {
        for (j = 0; j < best_cb->n_calls; j++) {
          if (best_cb->calls[j].pos > pos) break;
          if (best_cb->calls[j].pos != pos) continue;
          print_call (best_cb, j);
          fprintf (stdout, "\n");
        }
      }
      queue->last_chr = cb_f->chr;
      queue->last_pos = pos;
    }
}

static void
print_calls_all (CallBlock *cb_f, GASMQueue *queue)
{
  unsigned int pos;
  /* Print all relevant positions */
  for (pos = cb_f->start; pos < cb_f->end; pos++) {
    CallBlock *ccb;
    if ((cb_f->chr == queue->last_chr) && (pos <= queue->last_pos)) continue;
    /* Find call for given position with best aggregate value */
    for (ccb = queue->finished_blocks; ccb; ccb = ccb->next) {
      unsigned int j;
      if (ccb->chr != cb_f->chr) continue;
      for (j = 0; j < ccb->n_calls; j++) {
        if (ccb->calls[j].pos > pos) break;
        if (ccb->calls[j].pos != pos) continue;
        print_call (ccb, j);
        fprintf (stdout, "\n");
      }
    }
    queue->last_chr = cb_f->chr;
    queue->last_pos = pos;
  }
}

static void
print_calls (GASMQueue *queue)
{
  unsigned int min_chr_p = 0xffffffff;
  unsigned int min_start_p = 0xffffffff;
  CallBlock *cb;
  /* Find smallest processing */
  for (cb = queue->processing_blocks; cb; cb = cb->next) {
    if ((cb->chr < min_chr_p) || ((cb->chr == min_chr_p) && (cb->start < min_start_p))) {
      min_chr_p = cb->chr;
      min_start_p = cb->start;
    }
  }
  /* fprintf (stderr, "Smallest processing: %u %u\n", min_chr_p, min_start_p); */
  while (queue->finished_blocks) {
    CallBlock *cb_f = NULL;
    unsigned int min_chr_f = 0xffffffff;
    unsigned int min_start_f = 0xffffffff;
    /* Find block with smallest start address */
    for (cb = queue->finished_blocks; cb; cb = cb->next) {
      if ((cb->chr < min_chr_f) || ((cb->chr == min_chr_f) && (cb->start < min_start_f))) {
        min_chr_f = cb->chr;
        min_start_f = cb->start;
        cb_f = cb;
      }
    }
    if (!cb_f) return;
    if (cb_f->chr > min_chr_p) return;
    if ((cb_f->chr == min_chr_p) && (cb_f->end > min_start_p)) return;
    if (output == OUTPUT_POLY_BEST) {
      print_calls_poly_best (cb_f, queue, 1);
    } else if (output == OUTPUT_ALL_BEST) {
      print_calls_poly_best (cb_f, queue, 0);
    } else if (output == OUTPUT_ALL) {
      print_calls_all (cb_f, queue);
    }
    queue_free_call_block (queue, cb_f);
  }
}

static void
process (GT4Queue *queue, unsigned int idx, void *data)
{
  GASMQueue *gasm_queue = (GASMQueue *) queue;
  AssemblyData *adata = gasm_queue->adata[idx];
  gt4_queue_lock (queue);
  while (!gasm_queue->finished) {
    if ((gasm_queue->cpos < gasm_queue->csize) && (gasm_queue->line < max_regions)) {
      const unsigned char *tokenz[MAX_KMERS + 4];
      unsigned int lengths[MAX_KMERS + 4];
      unsigned int ntokenz;
      //fprintf (stderr, "process: Thread %u reading file from %llu\n", idx, gasm_queue->cpos);
      ntokenz = split_line (gasm_queue->cdata + gasm_queue->cpos, gasm_queue->csize - gasm_queue->cpos, tokenz, lengths, MAX_KMERS + 4);
      while ((gasm_queue->cpos < gasm_queue->csize) && (gasm_queue->cdata[gasm_queue->cpos] != '\n')) gasm_queue->cpos += 1;
      while ((gasm_queue->cpos < gasm_queue->csize) && (gasm_queue->cdata[gasm_queue->cpos] <= ' ')) gasm_queue->cpos += 1;
      gasm_queue->line += 1;
      //fprintf (stderr, "process: Thread %u line %u\n", idx, gasm_queue->line);
      if (ntokenz < 5) {
        fprintf (stderr, "process: Too few tokens at line %u\n", gasm_queue->line);
      } else {
        char *kmers[MAX_KMERS];
        unsigned int nkmers;
        unsigned int i;
        char chr[32];
        gasm_queue->nrunning += 1;
        if (lengths[0] > 31) lengths[0] = 31;
        memcpy (chr, tokenz[0], lengths[0]);
        chr[lengths[0]] = 0;
        nkmers = 0;
        for (i = 4; i < ntokenz; i++) {
          kmers[nkmers++] = strndup ((const char *) tokenz[i], lengths[i]);
        }
        adata->chr = gt4_chr_from_string (chr);
        adata->start = strtol ((const char *) tokenz[1], NULL, 10);
        adata->end = strtol ((const char *) tokenz[2], NULL, 10);
        adata->ref = (const char *) tokenz[3];
        adata->cblock = queue_get_call_block (gasm_queue, adata->chr, adata->start, adata->end);
        /* Now the new block is allocated so it is safe to print all completed blocks */
        print_calls (gasm_queue);
        gt4_queue_unlock (queue);
        assemble (adata, (const char **) kmers, nkmers, 0);
        gt4_queue_lock (queue);
        queue_finish_call_block (gasm_queue, adata->cblock);
        assembly_data_clear (adata);

        for (i = 0; i < nkmers; i++) free (kmers[i]);
        gasm_queue->nrunning -= 1;
      }
      gt4_queue_broadcast (queue);
    } else {
      if (!gasm_queue->nrunning) {
        gasm_queue->finished = 1;
      } else {
        gt4_queue_wait (queue);
      }
    }
  }
  //fprintf (stderr, "Thread %u exiting\n", idx);
  gt4_queue_broadcast (queue);
  gt4_queue_unlock (queue);
}

static SNV *read_snvs (const char *filename, unsigned int *n_snvs);
static SNV *read_fps (const char *filename, unsigned int *n_snvs);
/* Get index of SNV at or next to pos */
static unsigned int lookup_snv (SNV *snvs, unsigned int n_snvs, unsigned int chr, unsigned long long pos);

static GT4GmerDB *load_db_or_die (const char *db_name, const char *seq_dir, const char *id);
static SeqFile *map_sequences (GT4GmerDB *db, const char *seq_dir);

/* Get list of unique reads containing at lest one kmer from list */
static unsigned int get_unique_reads (ReadInfo reads[], unsigned int nreads, GT4GmerDB *db, SeqFile files[], const char *kmers[], unsigned int nkmers, unsigned int max_reads_per_kmer);
/* Get proper read sequences (forward direction) */
static unsigned int get_read_sequences (GASMRead *reads[], const ReadInfo read_info[], unsigned int nreads, SeqFile files[]);
static void print_db_reads (GT4Index *index, SeqFile *files, unsigned long long kmer_idx, unsigned int kmer_dir, FILE *ofs);

unsigned int smith_waterman_seq (unsigned int a_pos[], unsigned int b_pos[], const NSeq *a, const NSeq *b, SWCell *t, unsigned int debug);
static void print_alignment (FILE *ofs, unsigned int a_pos[], unsigned int b_pos[], unsigned int len, NSeq *a, NSeq *b);


/* Basic parameters */
const char *db_name = NULL;
const char *snv_db_name = NULL;
const char *fp_db_name = NULL;
const char *seq_dir = NULL;
static unsigned int print_reads = 0;
static unsigned int prefetch_db = 1;
static unsigned int prefetch_seq = 0;
/* Advanced parameters */
SNV *snvs = NULL;
SNV *fps = NULL;
unsigned int n_snvs = 0;
unsigned int n_fps = 0;
static unsigned int single_cutoff = 10;
static unsigned int min_confirming = 2;
static unsigned int min_group_coverage = 1;
static unsigned int max_divergent = 4;
static unsigned int min_align_len = 25;
static unsigned int min_group_size = 3;
static float min_group_rsize = 0.0f;
static unsigned int max_group_divergence = 3;
static unsigned int max_group_rdivergence = 3;
static unsigned int skip_end_align = 10;
static unsigned int skip_end_call = 10;
static unsigned int require_both_dirs = 1;

static void
print_usage (FILE *ofs, unsigned int advanced, int exit_value)
{
  fprintf (ofs, "gassembler version %u.%u.%u (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_MICRO, VERSION_QUALIFIER);
  fprintf (ofs, "Usage: gassembler --dbi FILENAME --region_file FILENAME [ARGUMENTS]\n");
  fprintf (ofs, "Common options:\n");
  fprintf (ofs, "    -v, --version                    - print version information and exit\n");
  fprintf (ofs, "    -h, --help                       - print this usage screen and exit\n");
  fprintf (ofs, "    --dbi FILENAME                   - index of sequenced reads (mandatory)\n");
  fprintf (ofs, "    --region_file FILENAME           - reference and kmer database (mandatory)\n");
  fprintf (ofs, "    --sex male|female|auto           - sex of the individual (default auto)\n");
  fprintf (ofs, "    --coverage FLOAT | median | local | ignore - average sequencing depth (default - median, local - use local number of reads)\n");
  fprintf (ofs, "    --num_threads                    - number of threads to use (default %u)\n", n_threads);
  fprintf (ofs, "    --min_p FLOAT                    - minimum call quality (default %.2f)\n", min_p);
  fprintf (ofs, "    --min_pmut FLOAT                 - minimum reference call quality (default %.2f)\n", min_pmut);
  fprintf (ofs, "    --exome                          - Disable quality models (needed if coverage variability is high)\n");
  fprintf (ofs, "    --advanced                       - print advanced usage options\n");
  if (advanced) {
    fprintf (ofs, "Advanced options:\n");
    fprintf (ofs, "    --seq_dir DIRECTORY              - directory of fastq files (overrides location in index)\n");
    fprintf (ofs, "    --region CHR START END SEQ       - call single reference region\n");
    fprintf (ofs, "    --min_coverage INTEGER           - minimum coverage for a call (default %u)\n", min_coverage);
    fprintf (ofs, "    --output poly | best | all       - output type (only polymorphisms, best calls for positon, all calls) (default poly)\n");
    fprintf (ofs, "    --counts                         - output nucleotide counts\n");
    fprintf (ofs, "    --extra                          - output extra information about call\n");
    fprintf (ofs, "    --min_confirming INTEGER         - minimum confirming nucleotide count for a call (default %u)\n", min_confirming);
    fprintf (ofs, "    --min_group_coverage INTEGER     - minimum coverage of group (default %u)\n", min_group_coverage);
    fprintf (ofs, "    --max_divergent INTEGER          - maximum number of mismatches per read (default %u)\n", max_divergent);
    fprintf (ofs, "    --min_align_len INTEGER          - minimum alignment length (default %u)\n", min_align_len);
    fprintf (ofs, "    --min_group_size INTEGER         - minimum group size (default %u)\n", min_group_size);
    fprintf (ofs, "    --min_group_rsize FLOAT          - minimum relative group size (default %.2f)\n", min_group_rsize);
    fprintf (ofs, "    --max_group_divergence INTEGER   - maximum divergence in group (default %u)\n", max_group_divergence);
    fprintf (ofs, "    --max_group_rdivergence INTEGER  - maximum relative divergence in group (default %u)\n", max_group_rdivergence);
    fprintf (ofs, "    --skip_end_align INTEGER         - skip nucleotides at region ends during alignment (default %u)\n", skip_end_align);
    fprintf (ofs, "    --skip_end_call INTEGER          - skip nucleotides at alignment ends (default %u)\n", skip_end_call);
    fprintf (ofs, "    --allow_one_dir                  - Allow calling if all confirming reads have the same dir\n");
    fprintf (ofs, "    --alternatives                   - output also homozygous variant for each heterozygous position\n");
    fprintf (ofs, "    --max_read_length INTEGER        - maximum length of reads (default %u)\n", max_read_length);
    fprintf (ofs, "    --max_reference_length INTEGER   - maximum length of reference region (default %u)\n", max_reference_length);
    fprintf (ofs, "    --error_prob FLOAT               - Probability of error (default %f)\n", error_prob);
    fprintf (ofs, "    --prefetch_seq                   - Prefetch FastQ sequences (slightly faster but uses more virtual memory/IO)\n");
    fprintf (ofs, "    --dont_prefetch_db               - Do not prefetch index (much slower but uses less memory/IO)\n");
    fprintf (ofs, "    -D                               - increase debug level\n");
    fprintf (ofs, "    -DG                              - increase group debug level\n");
    /*
     * fprintf (ofs, "    --snvs FILENAME                  - gmer_caller called SNVs\n");
     * fprintf (ofs, "    --fp FILENAME                    - List of known false positives\n");
     */
  }
  exit (exit_value);
}

static unsigned int only_chr = CHR_1;
static unsigned int only_pos = 0;

GT4Scout db_scout;
GT4Scout seq_scout;

int
main (int argc, const char *argv[])
{
  unsigned int i;
  const char *input_name = NULL;
  const char *kmers[MAX_KMERS];
  unsigned int nkmers = 0;
  unsigned int ref_chr = CHR_NONE;
  unsigned int ref_start = 0, ref_end = 0;
  const char *ref = NULL;

  GT4GmerDB *db;
  SeqFile *files;

  srand(1);

  //double sum = 0;
  //for (i = 0; i <= 57; i++) sum += dbinom (i, 48 + 57, 0.5);
  //fprintf (stderr, "%g\n", sum);
  //exit (0);
    
  for (i = 1; i < argc; i++) {
    if (!strcmp (argv[i], "-v") || !strcmp (argv[i], "--version")) {
      fprintf (stdout, "gassembler version %u.%u.%u (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_MICRO, VERSION_QUALIFIER);
      exit (0);
    } else if (!strcmp (argv[i], "-h") || !strcmp (argv[i], "--help")) {
      print_usage (stdout, 0, 0);
    } else if (!strcmp (argv[i], "--advanced")) {
      print_usage (stdout, 1, 0);
    } else if (!strcmp (argv[i], "-dbi") || !strcmp (argv[i], "-dbb") || !strcmp (argv[i], "-db") || !strcmp (argv[i], "--dbi")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      db_name = argv[i];
    } else if (!strcmp (argv[i], "--reference") || !strcmp (argv[i], "--region")) {
      if ((i + 4) >= argc) print_usage (stderr, 0, 1);
      ref_chr = gt4_chr_from_string (argv[i + 1]);
      if (!ref_chr) print_usage (stderr, 0, 1);
      ref_start = atoi (argv[i + 2]);
      ref_end = atoi (argv[i + 3]);
      ref = (const char *) argv[i + 4];
      i += 4;
    } else if (!strcmp (argv[i], "--snvs")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      snv_db_name = argv[i];
    } else if (!strcmp (argv[i], "--fp")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      fp_db_name = argv[i];
    } else if (!strcmp (argv[i], "--region_file") || !strcmp (argv[i], "--file")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      input_name = argv[i];
    } else if (!strcmp (argv[i], "--pos")) {
      char *p;
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      p = strchr (argv[i], ':');
      if (p) {
        char c[4] = { 0 };
        memcpy (c, argv[i], p - argv[i]);
        only_chr = gt4_chr_from_string (c);
        only_pos = strtol (p + 1, NULL, 10);
      } else {
        only_pos = strtol (argv[i], NULL, 10);
      }
    } else if (!strcmp (argv[i], "--max_regions")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      max_regions = strtol (argv[i], NULL, 10);
    } else if (!strcmp (argv[i], "--min_coverage")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      min_coverage = strtol (argv[i], NULL, 10);
    } else if (!strcmp (argv[i], "--sex")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      if (!strcmp (argv[i], "male")) {
        sex = SEX_MALE;
      } else if (!strcmp (argv[i], "female")) {
        sex = SEX_FEMALE;
      } else if (!strcmp (argv[i], "auto")) {
        sex = SEX_AUTO;
      } else {
        print_usage (stderr, 0, 1);
      }
    } else if (!strcmp (argv[i], "--error_prob")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      error_prob = atof (argv[i]);
    } else if (!strcmp (argv[i], "--min_confirming")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      min_confirming = strtol (argv[i], NULL, 10);
    } else if (!strcmp (argv[i], "--min_group_coverage")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      min_group_coverage = strtol (argv[i], NULL, 10);
    } else if (!strcmp (argv[i], "--max_divergent")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      max_divergent = strtol (argv[i], NULL, 10);
    } else if (!strcmp (argv[i], "--min_align_len")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      min_align_len = strtol (argv[i], NULL, 10);
    } else if (!strcmp (argv[i], "--min_group_size")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      min_group_size = strtol (argv[i], NULL, 10);
    } else if (!strcmp (argv[i], "--min_group_rsize")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      min_group_rsize = atof (argv[i]);
    } else if (!strcmp (argv[i], "--max_group_divergence")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      max_group_divergence = strtol (argv[i], NULL, 10);
    } else if (!strcmp (argv[i], "--max_group_rdivergence")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      max_group_rdivergence = strtol (argv[i], NULL, 10);
    } else if (!strcmp (argv[i], "--skip_end_align")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      skip_end_align = strtol (argv[i], NULL, 10);
    } else if (!strcmp (argv[i], "--skip_end_call")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      skip_end_call = strtol (argv[i], NULL, 10);
    } else if (!strcmp (argv[i], "--allow_one_dir")) {
      require_both_dirs = 0;
    } else if (!strcmp (argv[i], "--coverage")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      if (!strcmp (argv[i], "ignore")) {
        coverage = COVERAGE_IGNORE;
      } else if (!strcmp (argv[i], "local")) {
        coverage = COVERAGE_LOCAL;
      } else if (!strcmp (argv[i], "median")) {
        coverage = COVERAGE_MEDIAN;
      } else {
        coverage = (float) atof (argv[i]);
        if (!coverage) {
          fprintf (stderr, "Coverage has to be positive real value\n");
          exit (1);
        }
      }
    } else if (!strcmp (argv[i], "--min_p")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      min_p = atof (argv[i]);
    } else if (!strcmp (argv[i], "--min_pmut")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      min_pmut = atof (argv[i]);
    } else if (!strcmp (argv[i], "--exome")) {
      exome = 1;
    } else if (!strcmp (argv[i], "--num_threads")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      n_threads = strtol (argv[i], NULL, 10);
    } else if (!strcmp (argv[i], "--print_reads")) {
      print_reads = 1;
    } else if (!strcmp (argv[i], "--seq_dir")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      seq_dir = argv[i];
    } else if (!strcmp (argv[i], "--output")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      if (!strcmp (argv[i], "poly")) {
        output = OUTPUT_POLY_BEST;
      } else if (!strcmp (argv[i], "best")) {
        output = OUTPUT_ALL_BEST;
      } else if (!strcmp (argv[i], "all")) {
        output = OUTPUT_ALL;
      } else {
        print_usage (stderr, 0, 1);
      }
    } else if (!strcmp (argv[i], "--counts")) {
      print_extra = 1;
    } else if (!strcmp (argv[i], "--extra")) {
      print_extra = 2;
    } else if (!strcmp (argv[i], "--alternatives")) {
      alternative_calls = 1;
    } else if (!strcmp (argv[i], "--max_read_length")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      max_read_length = strtol (argv[i], NULL, 10);
    } else if (!strcmp (argv[i], "--max_reference_length")) {
      i += 1;
      if (i >= argc) print_usage (stderr, 0, 1);
      max_reference_length = strtol (argv[i], NULL, 10);
    } else if (!strcmp (argv[i], "--prefetch_seq")) {
      prefetch_seq = 1;
    } else if (!strcmp (argv[i], "--dont_prefetch_db")) {
      prefetch_db = 0;
    } else if (!strcmp (argv[i], "-D")) {
      debug += 1;
    } else if (!strcmp (argv[i], "-DG")) {
      debug_groups += 1;
    } else if (!strcmp (argv[i], "-ta")) {
      test_alignment (argv[i + 1], argv[i + 2]);
      exit (0);
    } else  {
      if (!isalpha (argv[i][0])) {
        fprintf (stderr, "Invalid argument %s\n", argv[i]);
        print_usage (stderr, 0, 1);
      }
      if (nkmers < MAX_KMERS) kmers[nkmers++] = argv[i];
    }
  }

  /* if (debug > debug_groups) debug_groups = debug; */
  if (sex == SEX_AUTO) prefetch_db = 0;

  /* Check arguments */
  if (!db_name) {
    print_usage (stderr, 0, 1);
  }
  if (!input_name && !ref) {
    print_usage (stderr, 0, 1);
  }

  /* Read databases */
  db = load_db_or_die (db_name, seq_dir, "reads");

  if (coverage == COVERAGE_MEDIAN) coverage = find_coverage (&db->index);

  if (snv_db_name) {
    fprintf (stderr, "Loading SNV database\n");
    snvs = read_snvs (snv_db_name, &n_snvs);
    fprintf (stderr, "Num SNVs %u\n", n_snvs);
  }

  if (fp_db_name) {
    fprintf (stderr, "Loading known false positives\n");
    fps = read_fps (fp_db_name, &n_fps);
    fprintf (stderr, "Num false positives %u\n", n_fps);
  }

  /* Set up file lists */
  if (debug) fprintf (stderr, "Loading read sequences\n");
  files = map_sequences (db, seq_dir);
  if (!files) {
    fprintf (stderr, "Cannot read sequences: terminating\n");
    exit (1);
  }

  if (sex == SEX_AUTO) {
    uint64_t j;
    unsigned long long sum[3] = { 0 };
    unsigned int count[3] = { 0 };
    double avg[3];
    fprintf (stderr, "Determine sex\n");
    for (j = 0; j < db->n_nodes; j++) {
      Node *node = &db->nodes[j];
      const char *name = db->names + node->name;
      unsigned int first_kmer = node->kmers;
      unsigned int nkmers = node->nkmers;
      /* 0 - autosome, 1 - X, 2 - Y */
      unsigned int klass = 0;
      if (name[0] == 'X') {
        klass = 1;
      } else if (name[0] == 'Y') {
        klass = 2;
      }
      for (i = 0; i < nkmers; i++) {
        unsigned int kmer_count;
        gt4_index_get_kmer_info (&db->index, first_kmer + i, &kmer_count);
        sum[klass] += kmer_count;
        count[klass] += 1;
      }
    }
    if (!count[1]) {
      fprintf (stderr, "No X kmers found, cannot determine sex (use --sex)\n");
      exit (1);
    }
    for (i = 0; i < 3; i++) {
      avg[i] = (double) sum[i] / count[i];
      fprintf (stderr, "Klass %u kmers %u sum %llu avg %.3f\n", i, count[i], sum[i], avg[i]);
    }
    if ((100 * avg[2] / avg[1]) < (avg[1] / avg[0])) {
      sex = SEX_FEMALE;
    } else {
      sex = SEX_MALE;
    }
    fprintf (stderr, "Sex: %s\n", (sex == SEX_MALE) ? "Male" : "Female");
  }

  if (input_name) {
    if (!only_pos) {
      GASMQueue queue;
      unsigned int i;
      memset (&queue, 0, sizeof (GASMQueue));
      queue.cdata = gt4_mmap (input_name, &queue.csize);
      if (!queue.cdata) {
        fprintf (stderr, "Cannot mmap input file %s\n", input_name);
        exit (1);
      }
      /* Print info */
      fprintf (stdout, "#KATK version: %u.%u.%u\n", VERSION_MAJOR, VERSION_MINOR, VERSION_MICRO);
      fprintf (stdout, "#KMer Database: %s\n", db_name);
      if (coverage >= 0) {
        fprintf (stdout, "#Coverage: %.2f\n", coverage);
      } else {
        fprintf (stdout, "#Coverage: local\n");
      }
      print_header (stdout);
      fprintf (stdout, "\n");
      queue_setup (&queue, n_threads);
      queue.cpos = 0;
      queue.line = 0;
      queue.last_chr = 0;
      queue.last_pos = 0;
      for (i = 0; i < n_threads; i++) {
        queue.adata[i] = assembly_data_new (db, files);
      }
      gt4_queue_create_threads (&queue.gt4_queue, process, NULL);
      process (&queue.gt4_queue, 0, NULL);
      gt4_queue_lock (&queue.gt4_queue);
      while (queue.gt4_queue.nthreads_running > 1) {
        gt4_queue_wait (&queue.gt4_queue);
      }
      print_calls (&queue);
      gt4_queue_unlock (&queue.gt4_queue);
    } else {
      /* Only single position */
      const unsigned char *cdata;
      unsigned long long csize, cpos;
      cdata = gt4_mmap (input_name, &csize);
      if (!cdata) {
        fprintf (stderr, "Cannot mmap input file %s\n", input_name);
        exit (1);
      }
      cpos = 0;
      while (cpos < csize) {
        const unsigned char *tokenz[MAX_KMERS + 4];
        unsigned int lengths[MAX_KMERS + 4];
        unsigned int ntokenz;
        ntokenz = split_line (cdata + cpos, csize - cpos, tokenz, lengths, MAX_KMERS + 4);
        while ((cpos < csize) && (cdata[cpos] != '\n')) cpos += 1;
        while ((cpos < csize) && (cdata[cpos] <= ' ')) cpos += 1;
        if (ntokenz < 5) {
          fprintf (stderr, "process: Too few tokens at line\n");
        } else {
          const char *kmers[MAX_KMERS];
          unsigned int chr, start, end, nkmers;
          const char *ref;
          unsigned int i;
          char chrc[32];
          if (lengths[0] > 31) lengths[0] = 31;
          memcpy (chrc, tokenz[0], lengths[0]);
          chrc[lengths[0]] = 0;
          chr = gt4_chr_from_string (chrc);
          if (chr != only_chr) continue;
          start = strtol ((const char *) tokenz[1], NULL, 10);
          if (start > only_pos) break;
          end = strtol ((const char *) tokenz[2], NULL, 10);
          if (end <= only_pos) continue;
          if ((end - start) > max_reference_length) {
            fprintf (stderr, "WARNING: Region %u-%u is longer than maximum allowed length (%u), skipping\n", start, end, max_reference_length);
            continue;
          }
          ref = (const char *) tokenz[3];
          nkmers = 0;
          for (i = 4; i < ntokenz; i++) {
            kmers[nkmers++] = strndup ((const char *) tokenz[i], lengths[i]);
          }
          assemble_recursive (db, files, chr, start, end, ref, kmers, nkmers);
        }
      }
    }
  } else {
    assemble_recursive (db, files, ref_chr, ref_start, ref_end, ref, kmers, nkmers);
  }

  if (prefetch_db) {
    gt4_delete_scout (&db_scout);
  }
  if (prefetch_seq) {
    //gt4_delete_scout (&seq_scout);
  }

  return 0;
}

static int
assemble_recursive (GT4GmerDB *db, SeqFile *files, unsigned int ref_chr, unsigned int ref_start, unsigned int ref_end, const char *ref, const char *kmers[], unsigned int nkmers)
{
  AssemblyData *adata = assembly_data_new (db, files);
  int result;
  unsigned int len;
  char *dup;
  adata->chr = ref_chr;
  adata->start = ref_start;
  adata->end = ref_end;
  len = ref_end - ref_start;
  dup = (char *) malloc (len + 1);
  strncpy (dup, ref, len);
  adata->ref = dup;
  adata->cblock = (CallBlock *) malloc (sizeof (CallBlock));
  memset (adata->cblock, 0, sizeof (CallBlock));
  adata->cblock->calls = (Call *) malloc (MAX_REFERENCE_LENGTH * 2 * sizeof (Call));
  memset (adata->cblock->calls, 0, MAX_REFERENCE_LENGTH * 2 * sizeof (Call));
  adata->cblock->chr = adata->chr;
  adata->cblock->start = adata->start;
  adata->cblock->end = adata->end;
  adata->cblock->haploid = (((sex == SEX_MALE) && ((ref_chr == CHR_X) || (ref_chr == CHR_Y))) || (ref_chr == CHR_MT));
  result = align (adata, kmers, nkmers);
  if (result > 0) {
    result = group (adata, 1);
  } else if (result == 0) {
    unsigned int mid;
    mid = (ref_start + ref_end) / 2;
    result = 0;
    result += assemble_recursive (db, files, ref_chr, ref_start, mid, ref, kmers, nkmers);
    result += assemble_recursive (db, files, ref_chr, mid, ref_end, ref + (mid - ref_start), kmers, nkmers);
  }
  assembly_data_clear (adata);
  free (adata->cblock);
  free (dup);
  return result;
}

static double
gt1_prob (const unsigned short counts[], unsigned int n0, unsigned int coverage)
{
  unsigned int i;
  double log_p = lgamma (coverage);
  for (i = A; i <= GAP; i++) {
    log_p -= lgamma (counts[i] + 1);
    if (i == n0) {
      log_p += (log (1 - error_prob) * counts[i]);
    } else {
      log_p += (log (error_prob / 4) * counts[i]);
    }
  }
  return exp (log_p);
}

static float
gt2_prob (const unsigned short counts[], unsigned int n0, unsigned int n1, unsigned int coverage)
{
  unsigned int i;
  double log_p = lgamma (coverage);
  for (i = A; i <= GAP; i++) {
    log_p -= lgamma (counts[i] + 1);
    if ((i == n0) || (i == n1)) {
      log_p += (log (0.5 - error_prob / 2) * counts[i]);
    } else {
      log_p += (log (error_prob / 3) * counts[i]);
    }
  }
  return exp (log_p);
}

static unsigned int
count_divergent_from_alignment (NSeq *a, NSeq *b, unsigned int a_p[], unsigned int b_p[], unsigned int align_len, unsigned int *n_gaps, unsigned int *s_gap, unsigned int *e_gap, unsigned int *gaps_total)
{
  unsigned int n_divergent, i;
  *n_gaps = 0;
  *gaps_total = 0;
  if ((a_p[0] > 0) && (b_p[0] > 0)) {
    /* Starts are unaligned */
    unsigned int gap_a = a_p[0];
    unsigned int gap_b = b_p[0];
    unsigned int min = (gap_a < gap_b) ? gap_a : gap_b;
    *n_gaps += 1;
    *s_gap = min;
    *gaps_total = *gaps_total + min;
  }
  if ((a_p[align_len - 1] < (a->len - 1)) && (b_p[align_len - 1] < (b->len - 1))) {
    /* Ends are unaligned */
    unsigned int gap_a = a->len - 1 - a_p[align_len - 1];
    unsigned int gap_b = b->len - 1 - b_p[align_len - 1];
    unsigned int min = (gap_a < gap_b) ? gap_a : gap_b;
    *n_gaps += 1;
    *e_gap = min;
    *gaps_total = *gaps_total + min;
  }
  n_divergent = *n_gaps;
  for (i = 0; i < align_len; i++) {
    if (a->pos[a_p[i]].nucl != b->pos[b_p[i]].nucl) n_divergent += 1;
  }
  return n_divergent;
}

struct _Group {
  unsigned long long tag;
  unsigned long long mask;
  unsigned int size;
  unsigned int included;
  unsigned int compat;
  unsigned int min_cov;
  unsigned int max_cov;
  unsigned short divergent;
  unsigned short dirs;
  unsigned int *consensus;
};

#define MAX_UNALIGNED_SIZE 5

static int
align (AssemblyData *adata, const char *kmers[], unsigned int nkmers)
{
  ReadInfo read_info[MAX_READS];
  short *alignment;
  unsigned int i;
  adata->ref_seq = n_seq_new_length (adata->ref, adata->end - adata->start, WORDLEN);
  /* Check reference length */
  if ((adata->end - adata->start) > MAX_REFERENCE_LENGTH) {
    fprintf (stderr, "align: reference length (%u) too big (max %u)\n", adata->end - adata->start, MAX_REFERENCE_LENGTH);
    return 0;
  }
  /* Get all unique reads */
  unsigned int max_reads_per_kmer = (adata->chr == CHR_MT) ? 2000 : MAX_READS_PER_KMER;
  adata->nreads = get_unique_reads (read_info, MAX_READS, adata->db, adata->files, kmers, nkmers, max_reads_per_kmer);
  if (debug > 1) fprintf (stderr, "Got %u unique reads\n", adata->nreads);
  /* Create actual sequences in correct direction */
  get_read_sequences (adata->reads, read_info, adata->nreads, adata->files);
  if (print_reads) {
    for (i = 0; i < adata->nreads; i++) {
      fprintf (stdout, ">Read_%u\n", i);
      fprintf (stdout, "%s\n", adata->reads[i]->seq);
    }
  }
  /* Sanitize */
  if (debug > 1) fprintf (stderr, "Number of usable reads: %u\n", adata->nreads);
  if (print_reads) {
    for (i = 0; i < adata->nreads; i++) {
      fprintf (stdout, ">Read_%u\n", i);
      fprintf (stdout, "%s\n", adata->reads[i]->seq);
    }
  }
  if (debug == 1) {
    fprintf (stderr, "Block: %s %u %u Reads: %u\n", chr_names[adata->chr], adata->start, adata->end, adata->nreads);
  }
  if (adata->nreads < MIN_READS) {
    if (debug) fprintf (stderr, "Final number of reads (%u) too low (min %u)\n", adata->nreads, MIN_READS);
    return -1;
  }
  /* Align all reads to reference */
  if (debug > 1) fprintf (stderr, "Aligning reads to reference...");
  alignment = (short *) malloc (MAX_ALIGNED_READS * MAX_REFERENCE_LENGTH * 2);
  adata->na = align_reads_to_reference (adata->ref_seq, adata->reads, adata->nreads, adata->aligned_reads, alignment, adata->sw_matrix);
  if (debug > 1) fprintf (stderr, "\n");
  /* Generate gapped alignment */
  adata->p_len = create_gapped_alignment (adata->ref_seq, adata->start, adata->aligned_reads, adata->na, alignment, adata->aligned_ref, adata->ref_pos, adata->alignment);
  /* Calculate totals */
  memset (adata->coverage, 0, adata->p_len * 2);
  memset (adata->nucl_counts, 0, adata->p_len * (GAP + 1) * 2);
  for (i = 0; i < adata->p_len; i++) {
    unsigned int j;
    for (j = 0; j < adata->na; j++) {
      if (adata->alignment[j * A_COLS + i] <= GAP) {
        int nucl = adata->alignment[j * A_COLS + i];
        adata->nucl_counts[i][nucl] += 1;
        adata->coverage[i] += 1;
      }
    }
  }
  /* Tag all reads by divergent positions */
  unsigned int n_divergent = 0;
  for (i = 0; i < adata->p_len; i++) {
    unsigned int j;
    unsigned int diverges = 0;
    unsigned int cutoff = (adata->coverage[i] >= single_cutoff) ? 2 : 1;
    for (j = 0; j <= GAP; j++) {
      if (j == adata->aligned_ref[i]) continue;
      if (j == N) continue;
      if (adata->nucl_counts[i][j] >= cutoff) diverges = 1;
    }
    if (diverges) {
      unsigned int known = 0;
      unsigned int ref_allele = 0, alt_allele = 0;
      if (n_divergent >= 21) {
        fprintf (stderr, "assemble: Too many divergent positions (max 21), ignoring the rest\n");
        break;
      }
      if (debug > 1) fprintf (stderr, "Divergent position: %u\n", adata->ref_pos[i]);
      if (snvs) {
        unsigned int snv = lookup_snv (snvs, n_snvs, adata->chr, adata->start + i);
        if ((snv < n_snvs) && (snvs[snv].chr == adata->chr) && (snvs[snv].pos == (adata->start + i))) {
          if (debug > 1) fprintf (stderr, "Known SNV %s (%c/%c)\n", snvs[snv].id, n2c[snvs[snv].ref_allele], n2c[snvs[snv].alt_allele]);
          known = 1;
          ref_allele = snvs[snv].ref_allele;
          alt_allele = snvs[snv].alt_allele;
        } else {
          if (debug > 1) fprintf (stderr, "Potential DeNovo\n");
        }
      }
      for (j = 0; j < adata->na; j++) {
        unsigned int ref = adata->aligned_ref[i];
        unsigned int nucl = adata->alignment[j * A_COLS + i];
        unsigned int mask = 7;
        /* Do not count single nucleotides */
        if ((nucl <= GAP) && (adata->nucl_counts[i][nucl] < cutoff)) mask = 0;
        /* N-s are counted same as reference */
        if (nucl == N) nucl = ref;
        /* Uncovered positions are counted as reference but masked */
        if (nucl > GAP) {
          nucl = ref;
          mask = 0;
        }
        adata->aligned_reads[j]->unknown = adata->aligned_reads[j]->unknown << 3;
        if (!known || ((nucl != ref_allele) && (nucl != alt_allele))) {
          adata->aligned_reads[j]->unknown |= 7;
        }
        /* Make 000 reference variant */
        nucl = nucl ^ ref;
        adata->aligned_reads[j]->tag = (adata->aligned_reads[j]->tag << 3) | nucl;
        adata->aligned_reads[j]->mask = (adata->aligned_reads[j]->mask << 3) | mask;
      }
      n_divergent += 1;
    }
  }
  free (alignment);

  return adata->nreads;
}

static int
group (AssemblyData *adata, unsigned int print)
{
  Group groups[MAX_ALIGNED_READS];
  unsigned int i, j, k;

  /* Recalculate totals */
  memset (adata->coverage, 0, adata->p_len * 2);
  memset (adata->nucl_counts, 0, adata->p_len * (GAP + 1) * 2);
  for (i = 0; i < adata->p_len; i++) {
    unsigned int j;
    for (j = 0; j < adata->na; j++) {
      int nucl = adata->alignment[j * A_COLS + i];
      if (nucl <= GAP) {
        adata->nucl_counts[i][nucl] += 1;
        adata->coverage[i] += 1;
      }
    }
  }
  /* Distribute reads into singleto groups */
  memset (groups, 0, sizeof (groups));
  unsigned int n_groups = adata->na;
  for (i = 0; i < adata->na; i++) {
    adata->aligned_reads[i]->group = i;
    groups[i].size = 1;
    groups[i].tag = adata->aligned_reads[i]->tag & adata->aligned_reads[i]->mask;
    groups[i].mask = adata->aligned_reads[i]->mask;
    groups[i].dirs = adata->aligned_reads[i]->dir;
  }
  if (debug > 1) {
    for (i = 0; i < n_groups; i++) fprintf (stderr, "%llx\t", groups[i].tag);
    fprintf (stderr, "\n");
    for (i = 0; i < n_groups; i++) fprintf (stderr, "%llx\t", groups[i].mask);
    fprintf (stderr, "\n");
  }
  memset (adata->is_compat, 0, n_groups * MAX_GROUPS);
  memset (adata->n_common, 0, n_groups * MAX_GROUPS * 2);

  while (n_groups > 1) {
    for (i = 0; i < n_groups; i++) {
      for (j = 0; j < n_groups; j++) {
        if (j == i) {
          adata->is_compat[i][j] = 0;
          adata->n_common[i][j] = 0;
          continue;
        }
        unsigned long long common = groups[i].mask & groups[j].mask;
        if ((groups[i].tag & common) == (groups[j].tag & common)) {
          /* Groups are compatible */
          adata->is_compat[i][j] = 1;
        } else {
          adata->is_compat[i][j] = 0;
        }
        adata->n_common[i][j] = 0;
        while (common != 0) {
          if (common & 7) adata->n_common[i][j] += 1;
          common = common >> 3;
        }
      }
    }
    unsigned int max_i = 0;
    unsigned int max_j = 0;
    unsigned int found = 0;
    for (i = 0; i < n_groups; i++) {
      for (j = i + 1; j < n_groups; j++) {
        if (adata->is_compat[i][j]) {
          if (!found) {
            max_i = i;
            max_j = j;
            found = 1;
          } else {
            if (adata->n_common[i][j] > adata->n_common[max_i][max_j]) {
              max_i = i;
              max_j = j;
            } else if (adata->n_common[i][j] == adata->n_common[max_i][max_j]) {
              if ((groups[i].size + groups[j].size) > (groups[max_i].size + groups[max_j].size)) {
                max_i = i;
                max_j = j;
              }
            }
          }
        }
      }
    }
    if (max_i == max_j) break;
    /* Can merge groups i and j */
    if (debug_groups) fprintf (stderr, "Merging groups %u (size %u) and %u (size %u) (common %u): %llx %llx %llx %llx -> ", max_i, groups[max_i].size, max_j, groups[max_j].size, adata->n_common[max_i][max_j], groups[max_i].tag, groups[max_i].mask, groups[max_j].tag, groups[max_j].mask);
    groups[max_i].tag = (groups[max_i].tag & groups[max_i].mask) | (groups[max_j].tag & groups[max_j].mask);
    groups[max_i].mask = groups[max_i].mask | groups[max_j].mask;
    groups[max_i].size += groups[max_j].size;
    groups[max_i].dirs |= groups[max_j].dirs;
    if (debug_groups) fprintf (stderr, "%llx %llx\n", groups[max_i].tag, groups[max_j].mask);
    for (i = 0; i < adata->na; i++) if (adata->aligned_reads[i]->group == max_j) adata->aligned_reads[i]->group = max_i;
    n_groups -= 1;
    groups[max_j].tag = groups[n_groups].tag;
    groups[max_j].mask = groups[n_groups].mask;
    groups[max_j].size = groups[n_groups].size;
    groups[max_j].dirs = groups[n_groups].dirs;
    for (i = 0; i < adata->na; i++) if (adata->aligned_reads[i]->group == n_groups) adata->aligned_reads[i]->group = max_j;
  }
  if (debug_groups) fprintf (stderr, "Num remaining groups: %u\n", n_groups);
  
  /* Calculate group min and max */
  for (i = 0; i < n_groups; i++) {
    groups[i].min_cov = adata->na;
    for (j = 0; j < adata->p_len; j++) {
      unsigned int cov = 0;
      for (k = 0; k < adata->na; k++) {
        if (adata->aligned_reads[k]->group != i) continue;
        if (adata->alignment[k * A_COLS + j] <= GAP) cov += 1;
      }
      if (cov < groups[i].min_cov) groups[i].min_cov = cov;
      if (cov > groups[i].max_cov) groups[i].max_cov = cov;
    }
    /* Calculate compat */
    for (j = 0; j < adata->na; j++) {
      unsigned long long common = groups[i].mask & adata->aligned_reads[j]->mask;
      if ((groups[i].tag & common) == (adata->aligned_reads[j]->tag & common)) {
        /* Read is compatible with this group */
        groups[i].compat += 1;
      }
    }
  }

  /* Create group consensus and calculate divergences */
  unsigned int *g_cons = (unsigned int *) malloc (n_groups * adata->p_len * 4);
  unsigned int last_aligned_ref = N;
  unsigned int last_consensus = N;
  for (j = 0; j < n_groups; j++) {
    groups[j].consensus = &g_cons[j * adata->p_len];
    for (i = 0; i < adata->p_len; i++) {
      unsigned int c[10] = { 0 };
      for (k = 0; k < adata->na; k++) {
        if (adata->aligned_reads[k]->group == j) c[adata->alignment[k * A_COLS + i]] += 1;
      }
      unsigned int best = adata->aligned_ref[i];
      for (k = 0; k <= GAP; k++) {
        if (k == N) continue;
        if ((adata->nucl_counts[i][k] > 1) && (c[k] > c[best])) best = k;
      }
      g_cons[j * adata->p_len + i] = best;
      if (best != adata->aligned_ref[i]) {
        unsigned int snv;
        if (debug_groups) fprintf (stderr, "Divergent position in group %u %u:%u\n", j, adata->chr, adata->ref_pos[i]);
        snv = lookup_snv (snvs, n_snvs, adata->chr, adata->start + i);
        if ((snv < n_snvs) && (snvs[snv].chr == adata->chr) && (snvs[snv].pos == (adata->start + i))) {
          if (debug_groups) fprintf (stderr, "Known SNV (%c/%c)\n", n2c[snvs[snv].ref_allele], n2c[snvs[snv].alt_allele]);
        } else {
          if (debug_groups) fprintf (stderr, "Potential DeNovo\n");
          if (((last_aligned_ref != GAP) || (adata->aligned_ref[i] != GAP)) && ((last_consensus != GAP) || (best != GAP))) {
            groups[j].divergent += 1;
          }
        }
      }
      last_aligned_ref = adata->aligned_ref[i];
      last_consensus = best;
    }
  }

  /* Sort groups by divergence/size */
  for (i = 0; i < n_groups; i++) {
    for (j = i + 1; j < n_groups; j++) {
      if ((groups[j].divergent < groups[i].divergent) || ((groups[j].divergent == groups[i].divergent) && (groups[j].size > groups[i].size))) {
        Group t = groups[i];
        groups[i] = groups[j];
        groups[j] = t;
        for (k = 0; k < adata->na; k++) {
          if (adata->aligned_reads[k]->group == i) {
            adata->aligned_reads[k]->group = j;
          } else {
            if (adata->aligned_reads[k]->group == j) adata->aligned_reads[k]->group = i;
          }
        }
      }
    }
  }

  if (debug_groups) {
    for (i = 0; i < n_groups; i++) fprintf (stderr, "%llu\t", groups[i].tag);
    fprintf (stderr, "\n");
    for (i = 0; i < n_groups; i++) fprintf (stderr, "%llu\t", groups[i].mask);
    fprintf (stderr, "\n");
  }

  if (debug_groups) {
    fprintf (stderr, "Read groups:");
    for (i = 0; i < adata->na; i++) fprintf (stderr, " %u:%u", i, adata->aligned_reads[i]->group);
    fprintf (stderr, "\n");
  }

  if (debug_groups > 0) {
    for (i = 0; i < n_groups; i++) {
      unsigned int j;
      if (debug_groups > 0) fprintf (stderr, "Group %u size %u divergent %u, min %u max %u tag %llx mask %llx\n", i, groups[i].size, groups[i].divergent, groups[i].min_cov, groups[i].max_cov, groups[i].tag, groups[i].mask);
      if (debug_groups > 1) {
        for (j = 0; j < adata->p_len; j++) fprintf (stderr, "%c", n2c[groups[i].consensus[j]]);
        fprintf (stderr, "\n");
        for (j = 0; j < adata->na; j++) {
          if (adata->aligned_reads[j]->group == i) {
            fprintf (stderr, "%s\n", adata->aligned_reads[j]->name);
          }
        }
      }
    }
  }

  /* Discard groups with too many divergent positions */
  unsigned int max_groups = 2;
  if (sex == SEX_MALE) {
    if ((adata->chr == CHR_X) || (adata->chr == CHR_Y)) {
      max_groups = 1;
    }
  }
  if (adata->chr == CHR_MT) max_groups = 1;


  unsigned int min_div = groups[0].divergent;
  for (i = 0; i < n_groups; i++) if (groups[i].divergent < min_div) min_div = groups[i].divergent;
  unsigned int good_groups[2];
  unsigned int n_included = 0;
  for (i = 0; i < n_groups; i++) {
    groups[i].included = (n_included < max_groups);
    if (require_both_dirs && (groups[i].dirs != 3)) {
      groups[i].included = 0;
      if (debug_groups > 0) fprintf (stderr, "Discarded group %u (%u): All reads have the same dir (%s)\n", i, groups[i].size, (groups[i].dirs == 2) ? "rev" : "fwd");
    }
    if (groups[i].min_cov < min_group_coverage) {
      groups[i].included = 0;
      if (debug_groups > 0) fprintf (stderr, "Discarded group %u (%u): Minimum coverage is 0\\n", i, groups[i].size);
    }
    if (groups[i].size < min_group_size) {
      groups[i].included = 0;
      if (debug_groups > 0) fprintf (stderr, "Discarded group %u (%u): size too small (%u < %u)\n", i, groups[i].size, groups[i].size, min_group_size);
    }
    if (groups[i].divergent > max_group_divergence) {
      groups[i].included = 0;
      if (debug_groups > 0) fprintf (stderr, "Discarded group %u (%u): too big divergence (%u > %u)\n", i, groups[i].size, groups[i].divergent, max_group_divergence);
    }
    if (groups[i].divergent > (min_div + max_group_rdivergence)) {
      groups[i].included = 0;
      if (debug_groups > 0) fprintf (stderr, "Discarded group %u (%u): too big relative divergence (%u > %u)\n", i, groups[i].size, groups[i].divergent, min_div + max_group_rdivergence);
    }
    if ((float) groups[i].size < (groups[0].size * min_group_rsize)) {
      groups[i].included = 0;
      if (debug_groups > 0) fprintf (stderr, "Discarded group %u (%u): relative size too small (%.2f < %.2f)\n", i, groups[i].size, (double) groups[i].size / groups[0].size, min_group_rsize);
    }
    if (groups[i].included) good_groups[n_included++] = i;
  }

  if (n_included < 1) {
    free (g_cons);
    return 0;
  }

  recalculate_and_call (adata, groups, n_groups, good_groups, n_included, max_groups == 1, print);
#if 0
  if (alternative_calls && (n_included == 2)) {
    /* Call the same position with one group */
    recalculate_and_call (adata, groups, n_groups, good_groups, 1, 0, print);
  }
#endif
  free (g_cons);
  
  return adata->p_len;
}

static void
recalculate_and_call (AssemblyData *adata, Group groups[], unsigned int n_groups, unsigned int good_groups[], unsigned int n_included, unsigned int haploid, unsigned int print)
{
  unsigned int max_cov_0 = groups[good_groups[0]].max_cov;
  unsigned int div_0 = groups[good_groups[0]].divergent;
  unsigned int compat_0 = groups[good_groups[0]].compat;
  unsigned int max_cov_1 = 0;
  unsigned int div_1 = 0;
  unsigned int compat_1 = 0;
  unsigned int compat_both = 0;
  unsigned int i, j;
  if (n_included > 1) {
    max_cov_1 = groups[good_groups[1]].max_cov;
    div_1 = groups[good_groups[1]].divergent;
    compat_1 = groups[good_groups[1]].compat;
    /* Calculate compatible common */
    for (j = 0; j < adata->na; j++) {
      unsigned long long common = groups[good_groups[0]].mask & adata->aligned_reads[j]->mask;
      if ((groups[good_groups[0]].tag & common) != (adata->aligned_reads[j]->tag & common)) continue;
      common = groups[good_groups[1]].mask & adata->aligned_reads[j]->mask;
      if ((groups[good_groups[1]].tag & common) != (adata->aligned_reads[j]->tag & common)) continue;
      /* Read is compatible with both groups */
      compat_both += 1;
    }
  }
    
  if (debug_groups > 0) {
    for (i = 0; i < n_groups; i++) {
      unsigned int j;
      if (debug_groups > 0) fprintf (stderr, "Group %u size %u divergent %u, min %u max %u, included %u\n", i, groups[i].size, groups[i].divergent, groups[i].min_cov, groups[i].max_cov, groups[i].included);
      if (debug_groups > 1) {
        for (j = 0; j < adata->p_len; j++) fprintf (stderr, "%c", n2c[groups[i].consensus[j]]);
        fprintf (stderr, "\n");
        for (j = 0; j < adata->na; j++) {
          if (adata->aligned_reads[j]->group == i) {
            fprintf (stderr, "%s\n", adata->aligned_reads[j]->name);
          }
        }
      }
    }
  }

  /* Recalculate totals */
  unsigned int max_coverage = 0;
  memset (adata->coverage, 0, adata->p_len * 2);
  memset (adata->nucl_counts, 0, adata->p_len * (GAP + 1) * 2);
  for (i = 0; i < adata->p_len; i++) {
    unsigned int j;
    for (j = 0; j < adata->na; j++) {
      unsigned int grp = adata->aligned_reads[j]->group;
      if (!groups[grp].included) continue;
      if (adata->alignment[j * A_COLS + i] <= GAP) {
        int nucl = adata->alignment[j * A_COLS + i];
        if (nucl != groups[grp].consensus[i]) continue;
        adata->nucl_counts[i][nucl] += 1;
        adata->coverage[i] += 1;
      }
    }
    if (adata->coverage[i] > max_coverage) max_coverage = adata->coverage[i];
  }
  unsigned int chr_coverage = max_coverage;
  if ((coverage > 0) && (adata->chr != CHR_MT)) {
    chr_coverage = coverage;
    if (sex == SEX_MALE) {
      if ((adata->chr == CHR_X) || (adata->chr == CHR_Y)) {
        chr_coverage /= 2;
      }
    }
  }

  /* Call */
  unsigned int last_call_pos = 0;
  unsigned int sub = 0;
  unsigned int call_alignment[1024];
  CallBlock *cb = adata->cblock;
  cb->n_calls = 0;
  cb->chr_cov = chr_coverage;
  for (i = skip_end_call; i < (adata->p_len - skip_end_call); i++) {
    CallExtra extra = { 0 };
    unsigned int hz = 0;
    if (adata->ref_pos[i] == last_call_pos) {
      sub += 1;
    } else {
      sub = 0;
    }
    last_call_pos = adata->ref_pos[i];
    extra.n_groups_total = n_groups;
    extra.n_groups = n_included;
    extra.div_0 = div_0;
    extra.div_1 = div_1;
    extra.max_cov_0 = max_cov_0;
    extra.max_cov_1 = max_cov_1;
    extra.compat_0 = compat_0;
    extra.compat_1 = compat_1;
    extra.compat_both = compat_both;
    extra.end_dist = (i < (adata->p_len - 1 - i)) ? i : adata->p_len - 1 - i;
    if (call (adata, cb, i, sub, &extra, call_alignment, 0)) {
      cb->calls[cb->n_calls].extra = extra;
      call_alignment[cb->n_calls] = i;
      hz = (cb->calls[cb->n_calls].nucl[0] != cb->calls[cb->n_calls].nucl[1]);
      cb->n_calls += 1;
      if (alternative_calls && hz) {
        /* Call also homozygote */
        if (call (adata, cb, i, sub, &extra, call_alignment, 1)) {
          cb->calls[cb->n_calls].extra = extra;
          call_alignment[cb->n_calls] = i;
          cb->n_calls += 1;
        }
      }
    } else {
      cb->n_calls += 1;
    }
  }

  if (print) {
    /* Output alignment */
    print_header (stdout);
    if (debug) {
      fprintf (stdout, "\t ");
      for (i = 0; i < n_included; i++) {
        unsigned int j;
        fprintf (stdout, "       ");
        for (j = 0; j < adata->na; j++) {
          if (adata->aligned_reads[j]->group == good_groups[i]) fprintf (stdout, "%c", 'A' + i);
        }
      }
    }
    fprintf (stdout, "\n");
    for (i = 0; i < cb->n_calls; i++) {
      print_call (cb, i);
      if (debug_groups) {
        /* Aligned reads */
        unsigned int a_i = call_alignment[i];
        fprintf (stdout, "\t%c", n2c[adata->aligned_ref[a_i]]);
        for (j = 0; j < n_groups; j++) {
          unsigned int k;
          fprintf (stdout, "  [%c%c] ", n2c[groups[j].consensus[a_i]], (groups[j].consensus[a_i] == adata->aligned_ref[a_i]) ? ' ' : '*');
          for (k = 0; k < adata->na; k++) {
            if (adata->aligned_reads[k]->group == j) fprintf (stdout, "%c", n2c[adata->alignment[k * A_COLS + a_i]]);
          }
        }
      }
      fprintf (stdout, "\n");
    }
  }
}

static unsigned int
call (AssemblyData *adata, CallBlock *cb, unsigned int a_pos, unsigned int sub, CallExtra *extra, unsigned int *call_alignment, unsigned int force_homozygote)
{
    Call *call = &cb->calls[cb->n_calls];
    unsigned int j;

    memset (call, 0, sizeof (Call));
    call->nucl[0] = call->nucl[1] = NONE;

    call->pos = adata->ref_pos[a_pos];
    call->sub = sub;
    call->ref = adata->aligned_ref[a_pos];

    if (call->ref == GAP) {
      /* Insertion, prev is reference nucleotide */
      call->prev_ref = adata->ref[call->pos - adata->start];
    } else {
      if (call->pos > adata->start) {
        call->prev_ref = adata->ref[call->pos - adata->start - 1];
      } else {
        call->prev_ref = '!';
      }
    }
    call->cov = adata->coverage[a_pos];
    for (j = A; j <= GAP; j++) {
      call->counts[j] = adata->nucl_counts[a_pos][j];
    }

    if (fps) {
      /* NC if known false positive */
      unsigned int fp = lookup_snv (fps, n_fps, adata->chr, adata->start + a_pos);
      if ((fp < n_fps) && (fps[fp].chr == adata->chr) && (fps[fp].pos == call->pos)) return 1;
    }

    /* Pick best and second best allele */
    unsigned int best0 = 0, best1 = 0, n;
    unsigned int best_n0 = A, best_n1 = A;
    for (n = A; n <= GAP; n++) {
      if (n == N) continue;
      if (call->counts[n] > best0) {
        best1 = best0;
        best_n1 = best_n0;
        best0 = call->counts[n];
        best_n0 = n;
      } else if (call->counts[n] > best1) {
        best1 = call->counts[n];
        best_n1 = n;
      }
    }
    /* NC if too small coverage of most numerous allele */
    if (best0 < min_confirming) return 1;

    unsigned int local_cov = cb->chr_cov;
    double p_hom, p_het;
    if (!exome) {
      p_hom = calc_p_select_diploid (call, extra, local_cov, best_n0, best_n0);
      p_het = calc_p_select_diploid (call, extra, local_cov, best_n0, best_n1);
    } else {
      p_hom = gt1_prob (call->counts, best_n0, adata->coverage[a_pos] - call->counts[N]);
      p_het = (best1 >= min_confirming) ? gt2_prob (call->counts, best_n0, best_n1, adata->coverage[a_pos] - call->counts[N]) : 0;
    }
    double sum_probs = p_hom + p_het;
    if (!sum_probs) sum_probs = 1;
    p_hom /= sum_probs;
    p_het /= sum_probs;
    double hzp = 1;
    if (coverage == COVERAGE_IGNORE) local_cov = call->cov;
    if (cb->haploid) {
      /* Haploid */
      best_n1 = best_n0;
      call->nucl[0] = call->nucl[1] = best_n0;
      if (!exome) {
        call->p = calc_p_select_haploid (call, extra, local_cov);
        call->q = calc_p_qual_haploid (call, extra, local_cov);
      } else {
        call->p = call->q = p_hom;
      }
    } else if (!best1 || force_homozygote) {
      /* Diploid homozygote */
      best_n1 = best_n0;
      call->nucl[0] = call->nucl[1] = best_n0;
      if (!exome) {
        call->p = calc_p_select_diploid (call, extra, local_cov, best_n0, best_n0);
        call->q = calc_p_qual_diploid (call, extra, local_cov);
      } else {
        call->p = call->q = p_hom;
      }
    } else {
      /* Heterozygote */
      if (p_het >= p_hom) {
        call->nucl[0] = (best_n0 < best_n1) ? best_n0 : best_n1;
        call->nucl[1] = (best_n0 < best_n1) ? best_n1 : best_n0;
        call->p = p_het;
      } else {
        assert (best0 >= best1);
        assert (call->counts[best_n0] >= call->counts[best_n1]);
        best_n1 = best_n0;
        call->nucl[0] = call->nucl[1] = best_n0;
        call->p = p_hom;
      }
      if (!exome) {
        call->q = calc_p_qual_diploid (call, extra, local_cov);
      } else {
        call->q = call->p;
      }
    }
    call->p_det = calc_p_mdetect (call, extra, local_cov);
    call->poly = ((call->nucl[0] != adata->aligned_ref[a_pos]) || (call->nucl[1] != adata->aligned_ref[a_pos]));
    extra->prob = 1;//best_prob;
    extra->rprob = call->q / sum_probs;
    extra->hzprob = hzp;
    
    return 1;
}

static int
assemble (AssemblyData *adata, const char *kmers[], unsigned int nkmers, unsigned int print)
{
  unsigned int i;
  int result;

  /* Print virtual command line to simplify debugging */
  if (debug > 1) {
    fprintf (stderr, "Arguments: -db %s --reference %s %u %u ", db_name, chr_names[adata->chr], adata->start, adata->end);
    for (i = adata->start; i < adata->end; i++) fprintf (stderr, "%c", adata->ref[i - adata->start]);
    for (i = 0; i < nkmers; i++) fprintf (stderr, " %s", kmers[i]);
    fprintf (stderr, "\n");
  }
  result = align (adata, kmers, nkmers);
  if (result > 0) {
    result = group (adata, print);
  }
  if (result <= 0) {
    /* Fill call block with NC-s */
    CallBlock *cb = adata->cblock;
    cb->n_calls = adata->end - adata->start - 2 * skip_end_align - 2 * skip_end_call;
    for (i = 0; i < cb->n_calls; i++) {
      memset (&cb->calls[i], 0, sizeof (Call));
      cb->calls[i].pos = adata->start + skip_end_align + skip_end_call + i;
      cb->calls[i].ref = adata->ref_seq->pos[skip_end_align + skip_end_call + i].nucl;
      cb->calls[i].prev_ref = '.';
    }
  }
  
  return result;
}

/*
 * Align all reads to reference, discard too divergent
 * 
 * Return number of remaining reads
 * a is an alignment table (num_reads * reference_length)
 */

static const unsigned int max_endgap = 1;
static const unsigned int max_gaps = 10;

static void
test_alignment (const char *a, const char *b)
{
  NSeq *a_seq, *b_seq;
  unsigned int ref_p[MAX_REFERENCE_LENGTH], read_p[MAX_READ_LENGTH];
  unsigned int align_len;
  SWCell *sw_matrix = (SWCell *) malloc (1000 * 1000 * sizeof (SWCell));
  a_seq = n_seq_new (a, 25);
  b_seq = n_seq_new (b, 25);
  align_len = smith_waterman_seq (ref_p, read_p, a_seq, b_seq, sw_matrix, 1);
  print_alignment (stdout, ref_p, read_p, align_len, a_seq, b_seq);
}


unsigned int
align_reads_to_reference (NSeq *ref_seq, GASMRead *reads[], unsigned int nreads, GASMRead *a_reads[], short *a, SWCell *sw_matrix)
{
  unsigned int i;
  unsigned int na = 0;
  assert (ref_seq->len <= MAX_REFERENCE_LENGTH);
  assert (nreads <= MAX_READS);
  for (i = 0; i < nreads; i++) {
    unsigned int ref_p[MAX_REFERENCE_LENGTH], read_p[MAX_READ_LENGTH];
    unsigned int align_len, n_divergent, n_gaps, s_gap = 0, e_gap = 0, gaps_total;
    int j, r_p;
    unsigned int last;
    
    /* Align seq to reference */
    align_len = smith_waterman_seq (ref_p, read_p, ref_seq, reads[i]->nseq, sw_matrix, 0 /* (i == 14) && (na == 9)*/);
    /* Calculate divergence */
    n_divergent = count_divergent_from_alignment (ref_seq, reads[i]->nseq, ref_p, read_p, align_len, &n_gaps, &s_gap, &e_gap, &gaps_total);
    if (debug > 1) fprintf (stderr, "Read %u: %u divergen %u gaps %u gap length start %u end %u\n", i, n_divergent, n_gaps, gaps_total, s_gap, e_gap);

    if (debug > 2) {
      fprintf (stderr, ">%u/%u\n", i, na);
      print_alignment (stderr, ref_p, read_p, align_len, ref_seq, reads[i]->nseq);
    }

    /* Ignore too divergent */
    if (n_divergent > max_divergent) {
      if (debug > 1) {
        fprintf (stderr, "Read %u: %s\n", i, reads[i]->seq);
        fprintf (stderr, "  has too many divergences: %u total, %u gaps (len = %u)\n", n_divergent, n_gaps, gaps_total);
      }
      continue;
    }
    /* Ignore too short alignments */
    if (align_len < min_align_len) {
      if (debug > 1) {
        fprintf (stderr, "Read %u: %s\n", i, reads[i]->seq);
        fprintf (stderr, "  has too short alignment: %u\n", align_len);
      }
      continue;
    }
    /* Ignore end and start gaps */
    if ((s_gap > max_endgap) || (e_gap > max_endgap)) {
      if (debug > 1) {
        fprintf (stderr, "Read %u: %s\n", i, reads[i]->seq);
        fprintf (stderr, "  has too long endgaps: %u/%u\n", s_gap, e_gap);
      }
      continue;
    }
    /* Ignore too many gaps */
    if (gaps_total > max_gaps) {
      /* fixme: Potential long indel */
      if (debug > 1) {
        fprintf (stderr, "Read %u: %s\n", i, reads[i]->seq);
        fprintf (stderr, "  has too long gaps: %u\n", gaps_total);
      }
      continue;
    }
    a_reads[na] = reads[i];

    /* fixme: Remove this */
    for (j = 0; j < ref_seq->len; j++) a[na * A_COLS + j] = -1000;
    /* Initial part */
    for (j = 0; j < ref_p[0]; j++) {
      int d = j - (int) ref_p[0];
      r_p = (int) read_p[0] + d;
      a[na * A_COLS + j] = (r_p < 0) ? BEFORE : UNKNOWN;
    }
    /* Aligned part */
    a[na * A_COLS + ref_p[0]] = read_p[0];
    last = ref_p[0];
    for (j = 1; j < align_len; j++) {
      unsigned int k;
      for (k = last + 1; k < ref_p[j]; k++) a[na * A_COLS + k] = a[na * A_COLS + last];
      if (ref_p[j] > ref_p[j - 1]) a[na * A_COLS + ref_p[j]] = read_p[j];
      last = ref_p[j];
    }
    /* Final part */
    for (j = ref_p[align_len - 1] + 1; j < ref_seq->len; j++) {
      int d = j - (int) ref_p[align_len - 1];
      r_p = (int) read_p[align_len - 1] + d;
      a[na * A_COLS + j] = (r_p >= a_reads[na]->nseq->len) ? AFTER : UNKNOWN;
    }
    for (j = 0; j < ref_seq->len; j++) {
      assert (a[na * A_COLS + j] >= -3);
      assert (a[na * A_COLS + j] < 1000);
    }
    na += 1;
    if (na >= MAX_ALIGNED_READS) {
      fprintf (stderr, "align_reads_to_reference: maximum number of aligned reads (%u) achieved\n", MAX_ALIGNED_READS);
      break;
    }
  }
  assert (na <= MAX_ALIGNED_READS);
  return na;
}

unsigned int
create_gapped_alignment (NSeq *ref_seq, unsigned int ref_start, GASMRead *a_reads[], unsigned int na, short *a, unsigned int aligned_ref[], int ref_pos[], short *ga)
{
  int ref_p, last_ref_p;
  int read_p[1024], last_read_p[1024];
  unsigned int p_len = 0;
  unsigned int i;
  /* Set read positions to alignment start (+ skip) */
  for (i = 0; i < na; i++) {
    read_p[i] = a[i * A_COLS + skip_end_align];
    last_read_p[i] = UNKNOWN;
  }
  ref_p = skip_end_align;
  last_ref_p = UNKNOWN;
  /* Iterate over reference (+ skip) */
  while (ref_p < (ref_seq->len - skip_end_align)) {
    /* Write alignment */
    /* Ref */
    if ((last_ref_p < 0) || (ref_p > last_ref_p)) {
      aligned_ref[p_len] = ref_seq->pos[ref_p].nucl;
      ref_pos[p_len] = ref_start + ref_p;
      last_ref_p = ref_p;
    } else {
      aligned_ref[p_len] = GAP;
      ref_pos[p_len] = ref_start + ref_p;
    }
    /* Reads */
    for (i = 0; i < na; i++) {
      if ((read_p[i] >= 0) && ((last_read_p[i] < 0) || (read_p[i] > last_read_p[i]))) {
        ga[i * A_COLS + p_len] = a_reads[i]->nseq->pos[read_p[i]].nucl;
        last_read_p[i] = read_p[i];
      } else if (read_p[i] >= 0) {
        ga[i * A_COLS + p_len] = GAP;
      } else {
        ga[i * A_COLS + p_len] = NONE;
      }
    }
    /* Advance */
    int rgap = 1;
    if (ref_p < (ref_seq->len - skip_end_align - 1)) {
      int next_ref_p = ref_p + 1;
      for (i = 0; i < na; i++) {
        int next_read_p = a[i * A_COLS + next_ref_p];
        if ((read_p[i] >= 0) && (next_read_p >= 0)) {
          int gap = next_read_p - read_p[i];
          if (gap > rgap) rgap = gap;
        }
      }
    }
    if (ref_p < (ref_seq->len - skip_end_align - 1)) {
      int next_ref_p = ref_p + 1;
      for (i = 0; i < na; i++) {
        int next_read_p = a[i * A_COLS + next_ref_p];
        if (next_read_p >= 0) {
          if (read_p[i] < 0) {
            if (rgap == 1) read_p[i] = next_read_p;
          } else if (read_p[i] < next_read_p) {
            int delta = next_read_p - read_p[i];
            if (delta == rgap) read_p[i] += 1;
          }
        } else {
          read_p[i] = next_read_p;
        }
      }
    }
    if (rgap == 1) ref_p += 1;
    p_len += 1;
  }
  return p_len;
}

static void
print_alignment (FILE *ofs, unsigned int a_pos[], unsigned int b_pos[], unsigned int len, NSeq *a, NSeq *b)
{
  int left, i, last_a, last_b;
  /* A */
  /* Left unaligned part */
  left = (a_pos[0] > b_pos[0]) ? a_pos[0] : b_pos[0];
  for (i = 0; i < left; i++) {
    int a_p = (int) a_pos[0] - (left - i);
    if (a_p >= 0) {
      fprintf (ofs, "%c", n2c[a->pos[a_p].nucl]);
    } else {
      fprintf (ofs, " ");
    }
  }
  /* Aligned part */
  last_a = a_pos[0];
  last_b = b_pos[0];
  for (i = 0; i < len; i++) {
    while ((int) b_pos[i] > last_b) {
      fprintf (ofs, "-");
      last_b += 1;
    }
    while (last_a <= (int) a_pos[i]) {
      fprintf (ofs, "%c", n2c[a->pos[last_a].nucl]);
      last_a += 1;
    }
    last_b = b_pos[i] + 1;
  }
  /* Final unaligned part */
  for (i = a_pos[len - 1] + 1; i < a->len; i++) {
    fprintf (ofs, "%c", n2c[a->pos[i].nucl]);
  }
  fprintf (ofs, "\n");
  
  /* Alignment */
  /* Left */
  for (i = 0; i < left; i++) {
    fprintf (ofs, " ");
  }
  /* Aligned part */
  last_a = a_pos[0];
  last_b = b_pos[0];
  for (i = 0; i < len; i++) {
    while ((int) b_pos[i] > last_b) {
      fprintf (ofs, " ");
      last_b += 1;
    }
    while ((int) a_pos[i] > last_a) {
      fprintf (ofs, " ");
      last_a += 1;
    }
    if (a->pos[a_pos[i]].nucl == b->pos[b_pos[i]].nucl) {
      fprintf (ofs, "|");
    } else {
      fprintf (stderr, " ");
    }
    last_a = a_pos[i] + 1;
    last_b = b_pos[i] + 1;
  }
  fprintf (ofs, "\n");

  /* B */
  /* Left unaligned part */
  for (i = 0; i < left; i++) {
    int b_p = (int) b_pos[0] - (left - i);
    if (b_p >= 0) {
      fprintf (ofs, "%c", n2c[b->pos[b_p].nucl]);
    } else {
      fprintf (ofs, " ");
    }
  }
  /* Aligned part */
  last_a = a_pos[0];
  last_b = b_pos[0];
  for (i = 0; i < len; i++) {
    while ((int) a_pos[i] > last_a) {
      fprintf (ofs, "-");
      last_a += 1;
    }
    while (last_b <= (int) b_pos[i]) {
      fprintf (ofs, "%c", n2c[b->pos[last_b].nucl]);
      last_b += 1;
    }
    last_a = a_pos[i] + 1;
  }
  /* Final unaligned part */
  for (i = b_pos[len - 1] + 1; i < b->len; i++) {
    fprintf (ofs, "%c", n2c[b->pos[i].nucl]);
  }
  fprintf (ofs, "\n");
}

/* Smith-Waterman algorithm */

#define M_SCORE 2
#define N_SCORE 0
#define MM_SCORE -3
#define GAP_OPEN_SCORE -4
#define GAP_SCORE -2

#define NROWS (n + 1)
#define NCOLS (m + 1)

#define TBACK 10

unsigned int
smith_waterman_seq (unsigned int a_pos[], unsigned int b_pos[], const NSeq *a, const NSeq *b, SWCell *t, unsigned int debug)
{
  unsigned int i, j;
  unsigned int n = a->len, m = b->len;
  /* Fill first column with zeroes */
  /* Fill first row with zeroes */
  /* Fill table starting from 1-st row and pick maximum score */
  int max_i = 0, max_j = 0;
  memset (t, 0, NCOLS * sizeof (SWCell));
  for (j = 0; j <= m; j++) {
    t[0 * NCOLS + j].left_gap_score = -1000;
    t[0 * NCOLS + j].top_gap_score = -1000;
  }
  for (i = 1; i <= n; i++) {
    memset (t + i * NCOLS, 0, sizeof (SWCell));
    t[i * NCOLS + 0].left_gap_score = -1000;
    t[i * NCOLS + 0].top_gap_score = -1000;
    for (j = 1; j <= m; j++) {
      memset (t + i * NCOLS + j, 0, sizeof (SWCell));
      int score = ((a->pos[i - 1].nucl >= N) || (b->pos[j - 1].nucl >= N)) ? N_SCORE : (a->pos[i - 1].nucl == b->pos[j - 1].nucl) ? M_SCORE : MM_SCORE;
      t[i * NCOLS + j].score = 0;
      /* Diagonal */
      if ((t[(i - 1) * NCOLS + (j - 1)].score + score) > 0) {
        t[i * NCOLS + j].score = t[(i - 1) * NCOLS + (j - 1)].score + score;
        t[i * NCOLS + j].sx = -1;
        t[i * NCOLS + j].sy = -1;
      }
#if 1
      /* Left */
      t[i * NCOLS + j].left_gap_score = t[i * NCOLS + j].score + GAP_OPEN_SCORE;
      t[i * NCOLS + j].left_gap_len = 0;
      if ((t[i * NCOLS + (j - 1)].left_gap_score + GAP_SCORE) > t[i * NCOLS + j].left_gap_score) {
        t[i * NCOLS + j].left_gap_score = t[i * NCOLS + (j - 1)].left_gap_score + GAP_SCORE;
        t[i * NCOLS + j].left_gap_len = t[i * NCOLS + (j - 1)].left_gap_len + 1;
      }
      if (t[i * NCOLS + j].left_gap_score >= t[i * NCOLS + j].score) {
        t[i * NCOLS + j].score = t[i * NCOLS + j].left_gap_score;
        t[i * NCOLS + j].sx = -t[i * NCOLS + j].left_gap_len;
        t[i * NCOLS + j].sy = 0;
      }
      t[i * NCOLS + j].top_gap_score = t[i * NCOLS + j].score + GAP_OPEN_SCORE;
      t[i * NCOLS + j].top_gap_len = 0;
      if ((t[(i - 1) * NCOLS + j].top_gap_score + GAP_SCORE) > t[i * NCOLS + j].top_gap_score) {
        t[i * NCOLS + j].top_gap_score = t[(i - 1) * NCOLS + j].top_gap_score + GAP_SCORE;
        t[i * NCOLS + j].top_gap_len = t[(i - 1) * NCOLS + j].top_gap_len + 1;
      }
      if (t[i * NCOLS + j].top_gap_score >= t[i * NCOLS + j].score) {
        t[i * NCOLS + j].score = t[i * NCOLS + j].top_gap_score;
        t[i * NCOLS + j].sx = 0;
        t[i * NCOLS + j].sy = -t[i * NCOLS + j].top_gap_len;
      }
#else
      /* Left */
      if ((t[i * NCOLS + (j - 1)].left_gap_score + GAP_SCORE) > 0) {
        t[i * NCOLS + j].left_gap_score = t[i * NCOLS + (j - 1)].left_gap_score + GAP_SCORE;
        t[i * NCOLS + j].left_gap_len = t[i * NCOLS + (j - 1)].left_gap_len + 1;
        if (t[i * NCOLS + j].left_gap_score >= t[i * NCOLS + j].score) {
          t[i * NCOLS + j].score = t[i * NCOLS + j].left_gap_score;
          t[i * NCOLS + j].sx = -t[i * NCOLS + j].left_gap_len;
          t[i * NCOLS + j].sy = 0;
        } else if ((t[i * NCOLS + j].score + GAP_OPEN_SCORE) >= t[i * NCOLS + j].left_gap_score) {
          t[i * NCOLS + j].left_gap_score = t[i * NCOLS + j].score + GAP_OPEN_SCORE;
          t[i * NCOLS + j].left_gap_len = 0;
        }
      } else {
        t[i * NCOLS + j].left_gap_score = t[i * NCOLS + j].score + GAP_OPEN_SCORE;
        t[i * NCOLS + j].left_gap_len = 0;
      }
      if ((t[(i - 1) * NCOLS + j].top_gap_score + GAP_SCORE) > 0) {
        t[i * NCOLS + j].top_gap_score = t[(i - 1) * NCOLS + j].top_gap_score + GAP_SCORE;
        t[i * NCOLS + j].top_gap_len = t[(i - 1) * NCOLS + j].top_gap_len + 1;
        if (t[i * NCOLS + j].top_gap_score >= t[i * NCOLS + j].score) {
          t[i * NCOLS + j].score = t[i * NCOLS + j].top_gap_score;
          t[i * NCOLS + j].sx = 0;
          t[i * NCOLS + j].sy = -t[i * NCOLS + j].top_gap_len;
        } else if ((t[i * NCOLS + j].score + GAP_OPEN_SCORE) >= t[i * NCOLS + j].top_gap_score) {
          t[i * NCOLS + j].left_gap_score = t[i * NCOLS + j].score + GAP_OPEN_SCORE;
          t[i * NCOLS + j].top_gap_len = 0;
        }
      } else {
        t[i * NCOLS + j].top_gap_score = t[i * NCOLS + j].score + GAP_OPEN_SCORE;
        t[i * NCOLS + j].top_gap_len = 0;
      }
#endif
      if (t[i * NCOLS + j].score > t[max_i * NCOLS + max_j].score) {
        max_i = i;
        max_j = j;
      }
    }
  }
  if (debug > 2) {
    fprintf (stderr, "    ");
    for (j = 0; j < m; j++) fprintf (stderr, "%c          ", n2c[b->pos[j].nucl]);
    fprintf (stderr, "\n");
    for (i = 0; i < n; i++) {
      fprintf (stderr, "%c ", n2c[a->pos[i].nucl]);
      for (j = 0; j < m; j++) {
        fprintf (stderr, "%3d(%2d/%2d)[%2d/%2d/%2d/%2d] ", t[(i + 1) * NCOLS + (j + 1)].score, t[(i + 1) * NCOLS + (j + 1)].sx, t[(i + 1) * NCOLS + (j + 1)].sy,
          t[(i + 1) * NCOLS + (j + 1)].left_gap_score, t[(i + 1) * NCOLS + (j + 1)].left_gap_len, t[(i + 1) * NCOLS + (j + 1)].top_gap_score, t[(i + 1) * NCOLS + (j + 1)].top_gap_len);
      }
      fprintf (stderr, "  %c\n", n2c[a->pos[i].nucl]);
    }
    fprintf (stderr, "    ");
    for (j = 0; j < m; j++) fprintf (stderr, "%c          ", n2c[b->pos[j].nucl]);
    fprintf (stderr, "\n");
  }
  unsigned int len = 0;
  while ((max_i > 0) && (max_j > 0)) {
    int sx = t[max_i * NCOLS + max_j].sx;
    int sy = t[max_i * NCOLS + max_j].sy;
    if (!sx && !sy) break;
    if (t[max_i * NCOLS + max_j].score < 1) break;
    if (sx && sy) {
      a_pos[len] = max_i - 1;
      b_pos[len] = max_j - 1;
      len += 1;
    }
    max_i += sy;
    max_j += sx;
  }
  /* Write result in reverse */
  for (i = 0; i < len / 2; i++) {
    unsigned int t = a_pos[i];
    a_pos[i] = a_pos[len - 1 - i];
    a_pos[len - 1 - i] = t;
    t = b_pos[i];
    b_pos[i] = b_pos[len - 1 - i];
    b_pos[len - 1 - i] = t;
  }
  if (debug > 2) {
    for (i = 0; i < len; i++) fprintf (stderr, "%c", n2c[a->pos[a_pos[i]].nucl]);
    fprintf (stderr, "\n");
    for (i = 0; i < len; i++) fprintf (stderr, "%c", n2c[b->pos[b_pos[i]].nucl]);
    fprintf (stderr, "\n");
  }
  return len;
}

static SNV *
read_snvs (const char *filename, unsigned int *n_snvs)
{
  const unsigned char *cdata;
  unsigned long long csize, cpos;
  unsigned int n_lines;
  SNV *snvs;
  *n_snvs = 0;
  if (!filename) return NULL;
  cdata = gt4_mmap (filename, &csize);
  if (!cdata) return NULL;
  /* Count lines */
  n_lines = 1;
  for (cpos = 0; cpos < csize; cpos++) {
    if (cdata[cpos] == '\n') n_lines += 1;
  }
  snvs = (SNV *) malloc (n_lines * sizeof (SNV));
  cpos = 0;
  while (cpos < csize) {
    if (cdata[cpos] != '#') {
      const unsigned char *tokenz[5];
      unsigned int lengths[5];
      unsigned int ntokenz;
      ntokenz = split_line (cdata + cpos, csize - cpos, tokenz, lengths, 5);
      if (ntokenz < 2) {
        fprintf (stderr, "read_snvs: too few tokens at line %u\n", *n_snvs);
      } else {
        const unsigned char *stok[4];
        unsigned int slen[4];
        char chr[32];
        split_line_chr (tokenz[0], lengths[0], stok, slen, 5, ':');
        if (slen[0] > 31) slen[0] = 31;
        memcpy (chr, stok[0], slen[0]);
        chr[slen[0]] = 0;
        snvs[*n_snvs].chr = gt4_chr_from_string (chr);
        if (!snvs[*n_snvs].chr) {
          static unsigned int warned = 0;
          if (!warned) {
            fprintf (stderr, "read_snvs: invalid chromosome name %s\n", chr);
            warned = 1;
          }
        } else {
          snvs[*n_snvs].pos = strtol ((const char *) stok[1], NULL, 10) - 1;
          /*snvs[*n_snvs].id = strndup ((const char *) tokenz[2], lengths[2]);*/
          snvs[*n_snvs].id = "*";
          snvs[*n_snvs].ref_allele = c2n (stok[3][0]);
          snvs[*n_snvs].alt_allele = c2n (stok[3][2]);
          snvs[*n_snvs].genotype = (tokenz[1][0] != 'A') || (tokenz[1][1] != 'A');
          /* if (debug) {
            fprintf (stderr, "SNV %llu %u%u\n", snvs[*n_snvs].pos, snvs[*n_snvs].ref_allele, snvs[*n_snvs].alt_allele);
          } */
          *n_snvs += 1;
        }
      }
    }
    while ((cpos < csize) && (cdata[cpos] != '\n')) cpos += 1;
    while ((cpos < csize) && (cdata[cpos] <= ' ')) cpos += 1;
  }
  return snvs;
}

static SNV *
read_fps (const char *filename, unsigned int *n_snvs)
{
  const unsigned char *cdata;
  unsigned long long csize, cpos;
  unsigned int n_lines;
  SNV *snvs;
  *n_snvs = 0;
  if (!filename) return NULL;
  cdata = gt4_mmap (filename, &csize);
  if (!cdata) return NULL;
  /* Count lines */
  n_lines = 1;
  for (cpos = 0; cpos < csize; cpos++) {
    if (cdata[cpos] == '\n') n_lines += 1;
  }
  snvs = (SNV *) malloc (n_lines * sizeof (SNV));
  memset (snvs, 0, n_lines * sizeof (SNV));
  cpos = 0;
  while (cpos < csize) {
    if (cdata[cpos] != '#') {
      const unsigned char *tokenz[5];
      unsigned int lengths[5];
      unsigned int ntokenz;
      ntokenz = split_line (cdata + cpos, csize - cpos, tokenz, lengths, 5);
      if (ntokenz < 2) {
        fprintf (stderr, "read_fps: too few tokens at line %u\n", *n_snvs);
      } else {
        const unsigned char *stok[4];
        unsigned int slen[4];
        char chr[32];
        split_line_chr (tokenz[0], lengths[0], stok, slen, 5, ':');
        if (slen[0] > 31) slen[0] = 31;
        memcpy (chr, stok[0], slen[0]);
        chr[slen[0]] = 0;
        snvs[*n_snvs].chr = gt4_chr_from_string (chr);
        if (!snvs[*n_snvs].chr) {
          static unsigned int warned = 0;
          if (!warned) {
            fprintf (stderr, "read_fps: invalid chromosome name %s\n", chr);
            warned = 1;
          }
        } else {
          snvs[*n_snvs].pos = strtol ((const char *) stok[1], NULL, 10);
          if (debug > 2) fprintf (stderr, "FP: %u %llu\n", snvs[*n_snvs].chr, snvs[*n_snvs].pos);
          *n_snvs += 1;
        }
      }
    }
    while ((cpos < csize) && (cdata[cpos] != '\n')) cpos += 1;
    while ((cpos < csize) && (cdata[cpos] <= ' ')) cpos += 1;
  }
  return snvs;
}

/* Get index of SNV at or next to pos */

static unsigned int
lookup_snv (SNV *snvs, unsigned int n_snvs, unsigned int chr, unsigned long long pos)
{
  unsigned int min, max, mid;
  min = 0;
  max = n_snvs;
  mid = (min + max) / 2;
  /* fprintf (stderr, "Looking up %u %llu\n", chr, pos); */
  while ((mid != min) && (mid != max)) {
    /* fprintf (stderr, "%u\t%u\t%u\t%u\t%llu\n", min, mid, max, snvs[mid].chr, snvs[mid].pos); */
    if (mid >= n_snvs) {
      break;
    } else if (snvs[mid].chr < chr) {
      min = mid;
    } else if (snvs[mid].chr > chr) {
      max = mid;
    } else if (snvs[mid].pos < pos) {
      min = mid;
    } else if (snvs[mid].pos > pos) {
      max = mid;
    } else {
      break;
    }
    mid = (min + max) / 2;
  }
  return mid;
}

static GT4GmerDB *
load_db_or_die (const char *db_name, const char *seq_dir, const char *id)
{
  const unsigned char *cdata;
  unsigned long long csize;
  GT4GmerDB *db;
  /* Read database */
  if (debug) fprintf (stderr, "Loading %s database %s... ", id, db_name);
  cdata = gt4_mmap (db_name, &csize);
  if (!cdata) {
    fprintf (stderr, "cannot mmap (no such file?)\n");
    exit (1);
  }
  if (prefetch_db) {
    db_scout.cdata = cdata;
    db_scout.csize = csize;
    gt4_scout_mmap (&db_scout);
  }
  db = gt4_gmer_db_new_from_binary (cdata, csize);
  if (!db) {
    fprintf (stderr, "cannot read (wrong file format?)\n");
    exit (1);
  }
  if (!db->index.read_blocks) {
    fprintf (stderr, "no index\n");
    exit (1);
  }
  if (debug) fprintf (stderr, "done\n");
  return db;
}

static char *
get_seq_name (const char *in_name, const char *seq_dir)
{
  if (!seq_dir) {
    return strdup (in_name);
  } else {
    unsigned int dir_len = (unsigned int) strlen (seq_dir);
    unsigned int len = (unsigned int) strlen (in_name);
    unsigned int i;
    for (i = 0; i < len; i++) {
      if (in_name[len - 1 - i] == '/') break;
    }
    if (i >= len) {
      char *c = (char *) malloc (dir_len + 1 + len + 1);
      memcpy (c, seq_dir, dir_len);
      c[dir_len] = '/';
      memcpy (c + dir_len + 1, in_name, len);
      c[dir_len + 1 + len] = 0;
      return c;
    } else {
      char *c = (char *) malloc (dir_len + 1 + i + 1);
      memcpy (c, seq_dir, dir_len);
      c[dir_len] = '/';
      memcpy (c + dir_len + 1, in_name + len - i, i);
      c[dir_len + 1 + i] = 0;
      return c;
    }
  }
}

static SeqFile *
map_sequences (GT4GmerDB *db, const char *seq_dir)
{
  unsigned int i;
  SeqFile *files = (SeqFile *) malloc (db->index.n_files * sizeof (SeqFile));
  memset (files, 0, db->index.n_files * sizeof (SeqFile));
  for (i = 0; i < db->index.n_files; i++) {
    files[i].name = get_seq_name (db->index.files[i], seq_dir);
    if (!files[i].cdata) {
      files[i].cdata = gt4_mmap (files[i].name, &files[i].csize);
      if (!files[i].cdata) {
        fprintf (stderr, "Cannot memory map %s\n", files[i].name);
        free (files);
        return NULL;
      }
      if (prefetch_seq) {
        //scout_mmap (files[i].cdata, files[i].csize);
      }
    }
  }
  return files;
}

#define max_reads_per_region 200

static unsigned int
get_unique_reads (ReadInfo reads[], unsigned int max_reads, GT4GmerDB *db, SeqFile files[], const char *kmers[], unsigned int nkmers, unsigned int max_reads_per_kmer)
{
  unsigned int nreads = 0, i;

  for (i = 0; i < nkmers; i++) {
    unsigned long long first_read;
    unsigned int n_reads, n_new_reads, j, kmer_dir;
    unsigned long long word, rword;
    unsigned int code, node_idx, node_kmer, kmer_idx;
  
    word = string_to_word (kmers[i], strlen (kmers[i]));
    rword = get_reverse_complement (word, strlen (kmers[i]));
    if (rword < word) word = rword;
    code = trie_lookup (&db->trie, word);
    if (!code) {
      fprintf (stderr, "No such kmer: %s\n", kmers[i]);
      exit (0);
    }
    kmer_dir = ((code & 0x80000000) != 0);
    if (debug > 1) fprintf (stderr, "Kmer %s word %llu code %u\n", kmers[i], word, code);
    code &= 0x7fffffff;
    node_idx = (code >> db->kmer_bits) - 1;
    node_kmer = code & ((1 << db->kmer_bits) - 1);
    kmer_idx = db->nodes[node_idx].kmers + node_kmer;
    if (debug > 1) fprintf (stderr, "Node %u kmer %u idx %u dir %u\n", node_idx, node_kmer, kmer_idx, kmer_dir);
    if (debug > 2) print_db_reads (&db->index, files, kmer_idx, kmer_dir, stderr);
    first_read = gt4_index_get_kmer_info (&db->index, kmer_idx, &n_reads);
    if (n_reads > max_reads_per_kmer) {
      if (debug > 1) fprintf (stderr, "Kmer %u has too many reads: %u\n", i, n_reads);
      continue;
    }
    if (debug > 1) fprintf (stderr, "Num reads %u\n", n_reads);
    n_new_reads = 0;
    for (j = 0; j < n_reads; j++) {
      unsigned long long name_pos;
      unsigned int file_idx, dir;
      unsigned int k;
      gt4_index_get_read_info (&db->index, first_read + j, &file_idx, &name_pos, &dir);
      for (k = 0; k < nreads; k++) {
        if ((reads[k].file_idx == file_idx) && (reads[k].name_pos == name_pos)) break;
      }
      if (k >= nreads) {
        n_new_reads += 1;
        if (debug > 1) fprintf (stderr, "Adding read %u dir %u\n", nreads, dir);
        reads[nreads].name_pos = name_pos;
        /* Reads are shared between k-kmers */
        reads[nreads].kmer_pos = 0;
        reads[nreads].file_idx = file_idx;
        /* fixme: What to do if two kmers have conflicting read directions? */
        reads[nreads].dir = (dir != kmer_dir);
        nreads += 1;
        if (nreads >= max_reads) {
          fprintf (stderr, "get_unique_reads: Maximum number of reads (%u) reached, ignoring the rest\n", max_reads);
          break;
        }
      } else {
        if (debug > 2) fprintf (stderr, "  Already registered as %u\n", k);
      }
    }
    if (debug > 1) fprintf (stderr, "Kmer %u %s reads %u new %u\n", i, kmers[i], n_reads, n_new_reads);
    if (nreads >= max_reads) break;
  }
  if (nreads > max_reads_per_region) {
    for (i = 0; i < max_reads_per_region; i++) {
      unsigned int p = (unsigned int) (rand() / (1.0 + RAND_MAX));
      ReadInfo ri = reads[p];
      reads[p] = reads[i];
      reads[i] = ri;
    }
    nreads = max_reads_per_region;
  }
  return nreads;
}

static unsigned int
get_read_sequences (GASMRead *seqs[], const ReadInfo reads[], unsigned int nreads, SeqFile files[])
{
  unsigned int i;
  /* Create actual sequences */
  /* Change directionality of reads if needed */
  for (i = 0; i < nreads; i++) {
    char name[2048], seq[2048];
    if (!files[reads[i].file_idx].cdata) {
      files[reads[i].file_idx].cdata = gt4_mmap (files[reads[i].file_idx].name, &files[reads[i].file_idx].csize);
      if (!files[reads[i].file_idx].cdata) {
        fprintf (stderr, "Cannot mmap %s\n", files[reads[i].file_idx].name);
        return 0;
      }
    }
    const unsigned char *p = files[reads[i].file_idx].cdata + reads[i].name_pos;
    unsigned int len = 0;
    while (p[len] != '\n') len += 1;
    memcpy (name, p, len);
    name[len] = 0;
    p += len;
    p += 1;
    len = 0;
    while (p[len] >= 'A') len += 1;
    if (len > max_read_length) {
      fprintf (stderr, "WARNING: Read is longer than maximum allowed length (%u, max %u), truncating\n", len, max_read_length);
      len = max_read_length;
    }
    memcpy (seq, p, len);
    seq[len] = 0;
    if (reads[i].dir) gt4_string_revcomp_inplace (seq, len);
    seqs[i] = gasm_read_new (name, seq, WORDLEN, reads[i].dir);
    if (debug > 1) fprintf (stderr, "Read %2u(%u): >%s\n%s\n", i, reads[i].dir, seqs[i]->name, seqs[i]->seq);
  }
  return 1;
}

static void
print_db_reads (GT4Index *index, SeqFile *files, unsigned long long kmer_idx, unsigned int kmer_dir, FILE *ofs)
{
  unsigned long long first_read;
  unsigned int i, num_reads;
  char c[2048];

  first_read = gt4_index_get_kmer_info (index, kmer_idx, &num_reads);
  fprintf (stderr, "Reads %llu first %llu num %u\n", index->read_blocks[kmer_idx], first_read, num_reads);

  for (i = 0; i < num_reads; i++) {
    unsigned long long kmer_pos, name_pos;
    unsigned int file_idx, dir;
    unsigned long long cpos;
    unsigned int j, k, len;

    kmer_pos = gt4_index_get_read_info (index, first_read + i, &file_idx, &name_pos, &dir);
    fprintf (stderr, "%u %s %u %llu %llu (dir %u)\n", i, index->files[file_idx], file_idx, name_pos, kmer_pos, dir);
    if (!files[file_idx].cdata) {
      files[file_idx].cdata = gt4_mmap (files[file_idx].name, &files[file_idx].csize);
    }
    const unsigned char *cdata = files[file_idx].cdata;

    /* Name */
    fprintf (ofs, ">");
    cpos = name_pos + 1;
    j = 0;
    while (cdata[cpos + j] >= ' ') {
      fprintf (ofs, "%c", cdata[cpos + j]);
      j += 1;
    }
    fprintf (ofs, "\n");

    while (cdata[cpos + j] < ' ') j++;
    len = 0;
    while (cdata[cpos + j] >= ' ') {
      c[len++] = cdata[cpos + j];
      j += 1;
    }
    c[len] = 0;
    if (kmer_dir != dir) {
      for (k = 0; k < (len / 2); k++) {
        char t = c[k];
        c[k] = c[len - 1 - k];
        c[len - 1 - k] = t;
      }
      for (k = 0; k < len; k++) {
        if (c[k] == 'A') c[k] = 'T';
        else if (c[k] == 'C') c[k] = 'G';
        else if (c[k] == 'G') c[k] = 'C';
        else if (c[k] == 'T') c[k] = 'A';
      }
    }
    fprintf (ofs, "%s", c);
    fprintf (ofs, "\n");
  }
}

static float
find_coverage (GT4Index *index)
{
  /* Find coverage */
  unsigned int MEDIAN_KMERS = 10000;
  unsigned int *counts = (unsigned int *) malloc (MEDIAN_KMERS * 4);
  unsigned int max = 0;
  unsigned int min = 0xffffffff;
  unsigned int i, ci = 0;
  while(ci < MEDIAN_KMERS) {
    unsigned int kmer_idx = rand () % index->n_kmers;
    gt4_index_get_kmer_info (index, kmer_idx, &counts[ci]);
    if (!counts[ci]) continue;
    if (counts[ci] < min) min = counts[ci];
    if (counts[ci] > max) max = counts[ci];
    ci += 1;
  }
  if (debug) fprintf (stderr, "Sample min %u max %u\n", min, max);
  unsigned int med = (min + max) / 2;
  while (max > min) {
    unsigned int above = 0, below = 0, equal;
    med = (min + max) / 2;
    for (i = 0; i < MEDIAN_KMERS; i++) {
      if (counts[i] < med) below += 1;
      if (counts[i] > med) above += 1;
    }
    equal = MEDIAN_KMERS - above - below;
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
  }
  if (debug) fprintf (stderr, "Sample median %u\n", med);
  if (debug > 1) {
    int bins[100] = { 0 };
    for (i = 0; i < MEDIAN_KMERS; i++) {
      bins[(counts[i] < 100) ? counts[i] : 99] += 1;
    }
    for (i = 0; i < 2 * med; i++) fprintf (stderr, "%u\t%u\n", i, bins[i]);
  }
  free (counts);
  return med;
}

static GASMRead *
gasm_read_new (const char *name, const char *seq, unsigned int wlen, unsigned int dir)
{
  GASMRead *read = (GASMRead *) malloc (sizeof (GASMRead));
  memset (read, 0, sizeof (GASMRead));
  read->name = strdup (name);
  read->seq = strdup (seq);
  read->nseq = n_seq_new (seq, wlen);
  read->dir = (1U << dir);
  return read;
}

static void
gasm_read_delete (GASMRead *read)
{
  free (read->name);
  free (read->seq);
  n_seq_delete (read->nseq);
  free (read);
}

static double
calc_p_select_diploid (Call *call, CallExtra *extra, unsigned int kmer_cov, unsigned int n0, unsigned int n1)
{
  if (exome) return call->cov / (call->cov + 0.25);
  /* Diloid selection model */
  double COMP_2 = extra->compat_both;
  double G0_COMP = extra->compat_0;
  double katvus = kmer_cov;
  double EDIST = extra->end_dist;
  double EDIST0 = (extra->end_dist == 0);
  double EDIST1 = (extra->end_dist == 1);
  double EDIST2 = (extra->end_dist == 2);
  double alternatiiv = (extra->n_groups_total > 1);
  double ignoreeri = (extra->n_groups_total != extra->n_groups);
  double max = (call->counts[n0] >= call->counts[n1]) ? call->counts[n0] : call->counts[n1];
  double all = call->counts[A] + call->counts[C] + call->counts[G] + call->counts[T] + call->counts[GAP];
  double kaugus1 = (call->cov - katvus) / sqrt(katvus);
  double kaugus2 = ((max - 0.5 * all) / sqrt(call->cov)) * (extra->n_groups >= 2);
  double suhe = max / (call->counts[A] + call->counts[C] + call->counts[G] + call->counts[T] + call->counts[GAP]) * (extra->n_groups != 1);
  double deletsioon2 = ((n0 == GAP) && (n1 == GAP));
  double deletsioon1 = (((n0 != GAP) && (n1 == GAP)) || ((n0 == GAP) && (n1 != GAP)));
  double HET = ((n0 == n1) && (n0 != GAP));
  
  double p = 1.549817e+01 +
    COMP_2/G0_COMP * 3.214268e+00 +
    HET * -1.603723e+01 +
    deletsioon1 * 4.057173e+00 +
    deletsioon2 * -1.295838e+01 +
    katvus * 3.327203e-01 +
    EDIST0 * -2.055305e+00 +
    EDIST1 * -1.914959e+00 +
    EDIST2 * -5.105844e-01 +
    EDIST * 5.987854e-02 +
    alternatiiv * -7.634908e-01 +
    kaugus1 * 1.563516e+00 +
    kaugus2 * -1.233070e+01 +
    (kaugus1 + 0.5) * (kaugus1 > (-0.5)) * -3.456876e-01 +
    (kaugus1 - 2) * (kaugus1 > 2) * -1.089758e-01 +
    (kaugus1 - 3) * (kaugus1 > 3) * -8.686674e-01 +
    kaugus2 * kaugus2 * -6.547970e-01 +
    G0_COMP / katvus * -1.655326e+00 +
    G0_COMP * G0_COMP / (katvus * katvus) * 2.113226e-01 +
    (EDIST - 40) * (EDIST - 40) * (EDIST - 40) * (EDIST > 40) * 2.992796e-03 +
    (EDIST - 45) * (EDIST - 45) * (EDIST - 45) * (EDIST > 45) * -6.197973e-03 +
    ignoreeri * -2.224370e-01 +
    suhe * -1.255600e+02 +
    suhe * suhe * 3.233437e+02 +
    suhe * suhe * suhe * -2.755079e+02 +
    suhe * suhe * suhe * suhe * 7.897496e+01 +
    EDIST * EDIST * -8.887499e-04 +
    HET * (EDIST <= 5) * -2.998684e-01 +
    COMP_2 / G0_COMP * katvus * -1.062955e-01 +
    HET * katvus * -2.855130e-01 +
    deletsioon1 * katvus * -9.098014e-02 +
    deletsioon2 * katvus * -2.018754e-01 +
    deletsioon2 * EDIST * 7.388170e-02 +
    deletsioon2 * alternatiiv * -4.950726e+00 +
    deletsioon2 * kaugus1 * -6.573440e-01 +
    deletsioon2 * kaugus2 * 1.337017e+01 +
    HET * kaugus2 * kaugus2 * 2.234410e+00 +
    HET * G0_COMP / katvus * 2.994476e+00 +
    HET * G0_COMP * G0_COMP / (katvus * katvus) * -4.286640e-01 +
    HET * kaugus1 * -8.026551e-01 +
    HET * kaugus2 * 9.614824e+00 +
    deletsioon1 * EDIST * -1.301157e-01 +
    EDIST * kaugus1 * -1.017782e-02 +
    kaugus1 * EDIST * EDIST * 1.413317e-04 +
    deletsioon1 * EDIST * EDIST * 2.472375e-03;
  p = exp(p);
  return (isfinite (p)) ? p / (1 + p) : 1;
}

static double
calc_p_select_haploid (Call *call, CallExtra *extra, unsigned int kmer_cov)
{
  if (exome) return call->cov / (call->cov + 0.25);
  /* Haploid selection model */
  double katvus = kmer_cov;
  double EDIST = extra->end_dist;
  double EDIST0 = (extra->end_dist == 0);
  double kaugus1 = (call->cov - katvus) / sqrt(katvus);

  double p = 2.734031375 +
    EDIST0 * -8.395304525 +
    ((EDIST == 1) || (EDIST == 2)) * -2.292773866 +
    (EDIST - 45) * (EDIST > 45) * 1.502826728 +
    kaugus1 * 0.617528244 +
    EDIST * kaugus1 * -0.009752782;
  p = exp(p);
  return (isfinite (p)) ? p / (1 + p) : 1;
}

static double
calc_p_qual_diploid (Call *call, CallExtra *extra, unsigned int kmer_cov)
{
  if (exome) return 1.0 + call->cov / (call->cov + 0.25);
#ifdef USE_SUB
  double SUB = call->sub;
#else
  double SUB = 0;
#endif
  double COMP_2 = extra->compat_both;
  double G0_COMP = extra->compat_0;
  double katvus = kmer_cov;
  double EDIST = extra->end_dist;
  double EDIST0 = (extra->end_dist == 0);
  double EDIST1 = (extra->end_dist == 1);
  double EDIST2 = (extra->end_dist == 2);
  double alternatiiv = (extra->n_groups_total > 1);
  double mitualternatiivi = (extra->n_groups_total > 2);
  double ignoreeri = (extra->n_groups_total != extra->n_groups);
  double max = (call->counts[call->nucl[0]] >= call->counts[call->nucl[1]]) ? call->counts[call->nucl[0]] : call->counts[call->nucl[1]];
  double all = call->counts[A] + call->counts[C] + call->counts[G] + call->counts[T] + call->counts[GAP];
  double kaugus1 = (call->cov - katvus) / sqrt(katvus);
  double kaugus2 = ((max - 0.5 * all) / sqrt(call->cov)) * (extra->n_groups >= 2);
  double suhe = max / (call->counts[A] + call->counts[C] + call->counts[G] + call->counts[T] + call->counts[GAP]) * (extra->n_groups != 1);
  double deletsioon2 = ((call->nucl[0] == GAP) && (call->nucl[1] == GAP));
  double deletsioon1 = (((call->nucl[0] != GAP) && (call->nucl[1] == GAP)) || ((call->nucl[0] == GAP) && (call->nucl[1] != GAP)));
  double HET = ((call->nucl[0] == call->nucl[1]) && (call->nucl[0] != GAP));

  double p = 5.625990e+00 +
    HET * -1.926639e+00 +
    deletsioon2 * -4.149465e+00 +
    kaugus1 * 1.976799e+00 +
    deletsioon1 * -3.674773e-01 +
    katvus * 2.505259e-01 +
    COMP_2 / G0_COMP * 3.530792e+00 +
    mitualternatiivi * 2.384205e-01 +
    alternatiiv * -1.893987e+00 +
    EDIST0 * 2.488365e+00 +
    EDIST1 * 3.614451e+00 +
    EDIST2 * -8.343540e-01 +
    (kaugus1 + 2) * (kaugus1 > (-2)) * -3.608020e-01 +
    (kaugus1 - 2) * (kaugus1 > (2)) * -1.369033e+00 +
    kaugus2 * -8.717219e-01 +
    (kaugus1 + 1) * (kaugus1 > (-1)) * -5.990449e-01 +
    G0_COMP / katvus * -5.090870e-01 +
    (EDIST - 35) * (EDIST > 35) * 7.200000e-02 +
    (EDIST - 30) * (EDIST > 30) * -6.277709e-02 +
    (EDIST - 45) * (EDIST > 45) * 1.407460e-01 +
    katvus * katvus * -3.807892e-03 +
    ignoreeri * -5.524936e-01 +
    1.0 * (SUB > 0) * -1.085515e+00 +
    HET * (EDIST < 5) * 1.155368e+00 +
    suhe * -1.489082e+02 +
    suhe * suhe * 6.542650e+02 +
    suhe * suhe * suhe * -9.392902e+02 +
    suhe * suhe * suhe * suhe * 4.360459e+02 +
    kaugus1 * deletsioon1 * -2.069432e-01 +
    HET * katvus * 1.598539e-01 +
    deletsioon2 * katvus * 2.304383e-01 +
    kaugus1 * katvus * -1.981619e-02 +
    deletsioon1 * katvus * 5.554233e-02 +
    deletsioon2 * alternatiiv * -5.609686e-01 +
    deletsioon2 * kaugus1 * 7.001617e-01 +
    deletsioon1 * kaugus1 * kaugus1 * 1.859963e-01 +
    HET * kaugus1 * 6.971654e-01 +
    HET * kaugus2 * -1.003972e-01 +
    HET * G0_COMP / katvus * -6.196470e-01 +
    HET * G0_COMP * G0_COMP / (katvus * katvus) * 1.267673e-01 +
    katvus * 1.0 * (SUB > 0) * -1.475575e-01 +
    kaugus1 * 1.0 * (SUB > 0) * -7.022790e-01 +
    HET * kaugus2 * kaugus2 * 2.281341e-01 +
    kaugus1 * deletsioon1 * katvus * 1.536606e-02 +
    HET * kaugus2 * 1.0 * (SUB > 0) * -5.997786e-01;
  p = exp(p);
  return (isfinite (p)) ? p / (1 + p) : 1;
}

static double
calc_p_qual_haploid (Call *call, CallExtra *extra, unsigned int kmer_cov)
{
  if (exome) return 1.0 + call->cov / (call->cov + 0.25);
  /* Homozygote quality model */
#ifdef USE_SUB
  double SUB = call->sub;
#else
  double SUB = 0;
#endif
  double EDIST = extra->end_dist;
  double HET = ((call->nucl[0] == call->nucl[1]) && (call->nucl[0] != GAP));

  double p = 7.7911387 +
    (EDIST - 45) * (EDIST > 45) * 0.7390936 +
    (SUB > 0) * -5.7026205 +
    HET * (EDIST < 5) * -0.9447409;
  p = exp(p);
  return (isfinite (p)) ? p / (1 + p) : 1;
}

static double
calc_p_mdetect (Call *call, CallExtra *extra, unsigned int kmer_cov)
{
  if (exome) return call->cov / (call->cov + 8.0);
  /* Mutation detectability model */
  double katvus = kmer_cov;
  double EDIST = extra->end_dist;
  double EDIST2 = (extra->end_dist == 2);
  double COV = call->cov;
  double kaugus1 = (call->cov - katvus) / sqrt(katvus);

  double p = -7.339851e+00 +
    kaugus1 * 2.457963e+00 +
    kaugus1 * kaugus1 * -2.092731e-01 +
    kaugus1 * kaugus1 * kaugus1 * 1.757365e-02 +
    EDIST * 1.174253e+00 +
    COV * 2.189787e-01 +
    katvus * 7.489705e-01 +
    COV * COV * -1.873808e-02 +
    COV * COV * COV * 2.716039e-04 +
    (kaugus1 + 0.5) * (kaugus1 > (-0.5)) * -5.814003e-01 +
    (kaugus1 - 3) * (kaugus1 > 3) * -8.967198e-02 +
    EDIST2 * 1.881940e+00 +
    EDIST * EDIST * -1.146688e-01 +
    EDIST * EDIST * EDIST * 4.807719e-03 +
    EDIST * EDIST * EDIST * EDIST * -9.036972e-05 +
    EDIST * EDIST * EDIST * EDIST * EDIST * 6.263128e-07 +
    kaugus1 * EDIST * -4.384856e-03 +
    COV * katvus * -3.309976e-02 +
    katvus * COV * COV * 9.086561e-04 +
    katvus * COV * COV * COV * -9.727565e-06 +
    EDIST * katvus * -9.141201e-05;
  p = exp(p);
  return (isfinite (p)) ? p / (1 + p) : 1;
}

