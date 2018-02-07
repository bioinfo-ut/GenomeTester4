#define __GASSEMBLER_C__

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <binomial.h>
#include <database.h>
#include <math.h>
#include <matrix.h>
#include <queue.h>
#include <sequence.h>
#include <utils.h>
#include <version.h>

int debug = 1;

#define WORDLEN 25
#define MAX_KMERS 1024
#define MAX_READS 4096
#define MAX_READS_PER_KMER 100

typedef struct _GASMRead GASMRead;

struct _GASMRead {
  char *name;
  char *seq;
  NSeq *nseq;
  unsigned long long tag;
  unsigned long long mask;
  unsigned long long unknown;
};

typedef struct _ReadInfo ReadInfo;
struct _ReadInfo {
  unsigned long long name_pos;
  unsigned int kmer_pos;
  unsigned int file_idx;
  unsigned int dir;
};

static GASMRead *gasm_read_new (const char *name, const char *seq, unsigned int wlen);
static void gasm_read_delete (GASMRead *read);
static unsigned int gasm_read_is_compatible (GASMRead *lhs, unsigned long long tag, unsigned long long mask);

static GASMRead *
gasm_read_new (const char *name, const char *seq, unsigned int wlen)
{
  GASMRead *read = (GASMRead *) malloc (sizeof (GASMRead));
  memset (read, 0, sizeof (GASMRead));
  read->name = strdup (name);
  read->seq = strdup (seq);
  read->nseq = n_seq_new (seq, wlen);
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

static unsigned int
gasm_read_is_compatible (GASMRead *lhs, unsigned long long tag, unsigned long long mask)
{
  unsigned long long common = lhs->mask & mask;
  return (lhs->tag & common) == (tag & common);
}

typedef struct _SeqFile SeqFile;

struct _SeqFile {
  const char *name;
  const unsigned char *cdata;
  unsigned long long csize;
};

typedef struct _SNV SNV;

struct _SNV {
  unsigned int chr;
  unsigned long long pos;
  char *id;
  unsigned short ref_allele;
  unsigned short alt_allele;
  unsigned short genotype;
};

static SNV *read_snvs (const char *filename, unsigned int *n_snvs);
/* Get index of SNV at or next to pos */
static unsigned int lookup_snv (SNV *snvs, unsigned int n_snvs, unsigned int chr, unsigned long long pos);

static void load_db_or_die (KMerDB *db, const char *db_name, const char *id);
static SeqFile *map_sequences (KMerDB *db);

static void assemble (KMerDB *db, SeqFile *files, KMerDB *gdb, SeqFile *g_files, const char *kmers[], unsigned int nkmers, unsigned int ref_chr, unsigned int ref_start, unsigned int ref_end, const char *ref);
static void assemble_by_matrix (KMerDB *db, KMerDB *gdb, const char *seqs[], unsigned int nreads, const char *ref_chr, unsigned int ref_start, unsigned int ref_end, const char *ref);

/* Get list of unique reads containing at lest one kmer from list */
static unsigned int get_unique_reads (ReadInfo reads[], unsigned int nreads, KMerDB *db, SeqFile files[], const char *kmers[], unsigned int nkmers);
/* Get proper read sequences (forward direction) */
static unsigned int get_read_sequences (GASMRead *reads[], const ReadInfo read_info[], unsigned int nreads, SeqFile files[]);
/* Remove reads that have multiple instances of same k-mer, too long gaps or too comon k-mers */
unsigned int remove_bad_reads (GASMRead *reads[], unsigned int nseqs, KMerDB *gdb);
static unsigned long long get_kmer_location (KMerDB *gdb, unsigned long long word, unsigned int *num_seqs, unsigned int *file_idx, unsigned int *dir);
static void print_db_reads (GT4Index *index, SeqFile *files, unsigned long long kmer_idx, unsigned int kmer_dir, FILE *ofs);
static void print_kmer_statistics (NMatrix *mat);

static unsigned int smith_waterman_seq (unsigned int a_pos[], unsigned int b_pos[], const NSeq *a, const NSeq *b);
static unsigned int smith_waterman (unsigned int a_pos[], unsigned int b_pos[], const unsigned int a[], unsigned int na, const unsigned int b[], unsigned int nb);
static void print_alignment (FILE *ofs, unsigned int a_pos[], unsigned int b_pos[], unsigned int len, NSeq *a, NSeq *b);
static unsigned int alignment_from_matrix (NCell *alignment[], unsigned int alignment_size, KMerDB *db, KMerDB *gdb, NMatrix *mat, unsigned int min_cov);
static void align_seq (NMatrix *mat, NCell *a[], unsigned int alen, unsigned int seq_idx);

enum { CHR_NONE, CHR_1, CHR_2, CHR_3, CHR_4, CHR_5, CHR_6, CHR_7, CHR_8, CHR_9, CHR_10, CHR_11, CHR_12, CHR_13, CHR_14, CHR_15, CHR_16, CHR_17, CHR_18, CHR_19, CHR_20, CHR_21, CHR_22, CHR_X, CHR_Y };
static const char *chr_names[] = { "INVALID", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y" };

static unsigned int
chr_from_text (const char *name)
{
  unsigned int val;
  char *e;
  if (!strcmp (name, "X")) return CHR_X;
  if (!strcmp (name, "Y")) return CHR_Y;
  val = strtol (name, &e, 10);
  if (*e) return CHR_NONE;
  if (val > CHR_22) return CHR_NONE;
  return val;
}

const char *db_name = NULL;
const char *gdb_name = NULL;
const char *snv_db_name = NULL;
static unsigned int print_chains = 0;
static unsigned int print_reads = 0, print_kmers = 0, analyze_kmers = 0;
static unsigned int min_coverage = 8;
SNV *snvs = NULL;
unsigned int n_snvs = 0;

static void
print_usage (int exit_value)
{
  fprintf (stderr, "gassembler version %u.%u (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_QUALIFIER);
  fprintf (stderr, "Usage: gassembler [OPTIONS] [KMERS...]\n");
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "    -v, --version           - print version information and exit\n");
  fprintf (stderr, "    -h, --help              - print this usage screen and exit\n");
  fprintf (stderr, "    -dbb, -db FILENAME      - name of read index file\n");
  fprintf (stderr, "    -gdb FILENAME           - name of genome index file\n");
  fprintf (stderr, "    --reference CHR START END SEQ - reference position and sequence\n");
  fprintf (stderr, "    --snvs FILENAME         - gmer_caller called SNVs\n");
  fprintf (stderr, "    --file FILENAME         - read reference and kmers from file (one line at time)\n");
  fprintf (stderr, "    -D                      - increase debug level\n");
  exit (exit_value);
}

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
  unsigned int max_regions = 1000000000;

  KMerDB db, gdb;
  SeqFile *files, *g_files;
    
  for (i = 1; i < argc; i++) {
    if (!strcmp (argv[i], "-v") || !strcmp (argv[i], "--version")) {
      fprintf (stdout, "gassembler version %d.%d (%s)\n", VERSION_MAJOR, VERSION_MINOR, VERSION_QUALIFIER);
    } else if (!strcmp (argv[i], "-h") || !strcmp (argv[i], "--help")) {
      print_usage (0);
    } else if (!strcmp (argv[i], "-dbb") || !strcmp (argv[i], "-db")) {
      i += 1;
      if (i >= argc) exit (1);
      db_name = argv[i];
    } else if (!strcmp (argv[i], "-gdb")) {
      i += 1;
      if (i >= argc) exit (1);
      gdb_name = argv[i];
    } else if (!strcmp (argv[i], "--reference")) {
      if ((i + 4) >= argc) exit (1);
      ref_chr = chr_from_text (argv[i + 1]);
      if (!ref_chr) exit (1);
      ref_start = atoi (argv[i + 2]);
      ref_end = atoi (argv[i + 3]);
      ref = (const char *) argv[i + 4];
      i += 4;
    } else if (!strcmp (argv[i], "--snvs")) {
      i += 1;
      if (i >= argc) exit (1);
      snv_db_name = argv[i];
    } else if (!strcmp (argv[i], "--file")) {
      i += 1;
      if (i >= argc) exit (1);
      input_name = argv[i];
    } else if (!strcmp (argv[i], "--max_regions")) {
      i += 1;
      if (i >= argc) exit (1);
      max_regions = strtol (argv[i], NULL, 10);
    } else if (!strcmp (argv[i], "--print_reads")) {
      print_reads = 1;
    } else if (!strcmp (argv[i], "--print_kmers")) {
      print_kmers = 1;
    } else if (!strcmp (argv[i], "--print_chains")) {
      print_chains = 1;
    } else if (!strcmp (argv[i], "--analyze_kmers")) {
      analyze_kmers = 1;
    } else if (!strcmp (argv[i], "-D")) {
      debug += 1;
    } else  {
      if (nkmers < MAX_KMERS) kmers[nkmers++] = argv[i];
    }
  }

  /* Check arguments */
  if (!db_name || !gdb_name) {
    print_usage (1);
  }

  /* Read databases */
  load_db_or_die (&db, db_name, "reads");
  load_db_or_die (&gdb, gdb_name, "genome");

  if (snv_db_name) {
    fprintf (stderr, "Loading SNV database\n");
    snvs = read_snvs (snv_db_name, &n_snvs);
    fprintf (stderr, "Num SNVs %u\n", n_snvs);
  }

  /* Set up file lists */
  fprintf (stderr, "Loading read sequences\n");
  files = map_sequences (&db);
  fprintf (stderr, "Loading genome sequence\n");
  g_files = map_sequences (&gdb);
  if (!files || !g_files) {
    fprintf (stderr, "Terminating\n");
    exit (1);
  }

  if (input_name) {
    const unsigned char *cdata;
    unsigned long long csize, cpos;
    unsigned int line;
    cdata = gt4_mmap (input_name, &csize);
    if (!cdata) {
      fprintf (stderr, "Cannot mmap input file %s\n", input_name);
      exit (1);
    }
    cpos = 0;
    line = 0;
    while (cpos < csize) {
      const unsigned char *tokenz[MAX_KMERS + 4];
      unsigned int lengths[MAX_KMERS + 4];
      unsigned int ntokenz;
      ntokenz = split_line (cdata + cpos, csize - cpos, tokenz, lengths, MAX_KMERS + 4);
      if (ntokenz < 5) {
        fprintf (stderr, "Too few tokens at line %u\n", line);
      } else {
        unsigned int i;
        char chr[32];
        if (lengths[0] > 31) lengths[0] = 31;
        memcpy (chr, tokenz[0], lengths[0]);
        chr[lengths[0]] = 0;
        ref_chr = chr_from_text (chr);
        ref_start = strtol ((const char *) tokenz[1], NULL, 10);
        ref_end = strtol ((const char *) tokenz[2], NULL, 10);
        ref = strndup ((const char *) tokenz[3], lengths[3]);
        nkmers = 0;
        for (i = 4; i < ntokenz; i++) {
          kmers[nkmers++] = strndup ((const char *) tokenz[i], lengths[i]);
        }
        assemble (&db, files, &gdb, g_files, kmers, nkmers, ref_chr, ref_start, ref_end, ref);
      }
      while ((cpos < csize) && (cdata[cpos] != '\n')) cpos += 1;
      while ((cpos < csize) && (cdata[cpos] <= ' ')) cpos += 1;
      line += 1;
      if (line >= max_regions) break;
    }
  } else {
    assemble (&db, files, &gdb, g_files, kmers, nkmers, ref_chr, ref_start, ref_end, ref);
  }

  return 0;
}

#define ERROR_PROB 0.01f

static float
gt1_prob (unsigned int gt_count, unsigned int total_count)
{
  double q0, q1;
  unsigned int err_count = total_count - gt_count;
  q0 = poisson (err_count, ERROR_PROB);
  q1 = poisson (gt_count, total_count);
  return (float) (q0 * q1);
}

static float
gt2_prob (unsigned int gt1_count, unsigned int gt2_count, unsigned int total_count)
{
  double q0, q1, q2;
  unsigned int err_count = total_count - gt1_count - gt2_count;
  q0 = poisson (err_count, ERROR_PROB);
  q1 = poisson (gt1_count, total_count / 2.0);
  q2 = poisson (gt2_count, total_count / 2.0);
  return (float) (q0 * q1 * q2);
}

#define MAX_UNALIGNED_SIZE 5

static void
assemble (KMerDB *db, SeqFile *files, KMerDB *gdb, SeqFile *g_files, const char *kmers[], unsigned int nkmers, unsigned int ref_chr, unsigned int ref_start, unsigned int ref_end, const char *ref)
{
  ReadInfo read_info[MAX_READS];
  GASMRead *reads[MAX_READS];
  unsigned int nreads, i, j, k;

  /* Print virtual command line to simplify debugging */
  if (debug) {
    fprintf (stderr, "Arguments: -db %s -gdb %s -reference %s %u %u %s", db_name, gdb_name, chr_names[ref_chr], ref_start, ref_end, ref);
    for (i = 0; i < nkmers; i++) {
      fprintf (stderr, " %s", kmers[i]);
    }
    fprintf (stderr, "\n");
  }

  /* Get all unique reads */
  nreads = get_unique_reads (read_info, MAX_READS, db, files, kmers, nkmers);
  if (debug) fprintf (stderr, "Got %u unique reads\n", nreads);
  /* Create actual sequences in correct direction */
  get_read_sequences (reads, read_info, nreads, files);

  /* Now we have bunch of reads */
  if (print_reads) {
    for (i = 0; i < nreads; i++) {
      fprintf (stdout, ">Read_%u\n", i);
      fprintf (stdout, "%s\n", reads[i]->seq);
    }
  }
  if (debug) fprintf (stderr, "Number of unique reads: %u\n", nreads);

  /* Sanitize */
  nreads = remove_bad_reads (reads, nreads, gdb);

  if (debug) fprintf (stderr, "Number of usable reads: %u\n", nreads);
  /* Now we have bunch of reads */
  if (print_reads) {
    for (i = 0; i < nreads; i++) {
      fprintf (stdout, ">Read_%u\n", i);
      fprintf (stdout, "%s\n", reads[i]->seq);
    }
  }
  if (nreads < 10) {
    fprintf (stderr, "Final number of reads too low (%u)\n", nreads);
    return;
  } else if (nreads > MAX_READS) {
    fprintf (stderr, "Final number of reads too big (%u)\n", nreads);
    return;
  }

  /* Align all reads to reference */
  if (debug) fprintf (stderr, "Aligning reads to reference");
  GASMRead *a_reads[MAX_READS];
  NSeq *ref_seq = n_seq_new (ref, WORDLEN);
  int *a[1024];
  unsigned int na = 0;
  for (i = 0; i < nreads; i++) {
    unsigned int ref_p[2048], read_p[2048], align_len;
    a_reads[na] = reads[i];
    align_len = smith_waterman_seq (ref_p, read_p, ref_seq, a_reads[na]->nseq);

    /* Ignore too short alignments */
    if (align_len > 25) {
      if (((ref_p[0] > MAX_UNALIGNED_SIZE) && (read_p[0] > MAX_UNALIGNED_SIZE)) || ((ref_p[align_len - 1] < (ref_seq->len - MAX_UNALIGNED_SIZE)) && (read_p[align_len - 1] < (a_reads[na]->nseq->len - MAX_UNALIGNED_SIZE)))) {
        /* Potential long indel */
        if (debug == 1) fprintf (stderr, "*");
      } else {
        int j, r_p;
        unsigned int last;

        if (debug > 1) {
          fprintf (stderr, ">%u/%u\n", i, na);
          print_alignment (stderr, ref_p, read_p, align_len, ref_seq, a_reads[na]->nseq);
        }

        a[na] = (int *) malloc (ref_seq->len * 4);
        /* fixme: Remove this */
        for (j = 0; j < ref_seq->len; j++) a[na][j] = -1000;
        /* Initial part */
        for (j = 0; j < ref_p[0]; j++) {
          int d = j - (int) ref_p[0];
          r_p = (int) read_p[0] + d;
          a[na][j] = (r_p < 0) ? BEFORE : UNKNOWN;
        }
        /* Aligned part */
        a[na][ref_p[0]] = read_p[0];
        last = ref_p[0];
        for (j = 1; j < align_len; j++) {
          unsigned int k;
          for (k = last + 1; k < ref_p[j]; k++) a[na][k] = a[na][last];
          if (ref_p[j] > ref_p[j - 1]) a[na][ref_p[j]] = read_p[j];
          last = ref_p[j];
        }
        /* Final part */
        for (j = ref_p[align_len - 1] + 1; j < ref_seq->len; j++) {
          int d = j - (int) ref_p[align_len - 1];
          r_p = (int) read_p[align_len - 1] + d;
          a[na][j] = (r_p >= a_reads[na]->nseq->len) ? AFTER : UNKNOWN;
        }
        for (j = 0; j < ref_seq->len; j++) {
          assert (a[na][j] >= -3);
          assert (a[na][j] < 1000);
        }
        na += 1;
      }
    }
    if (debug == 1) fprintf (stderr, ".");
  }
  if (debug == 1) fprintf (stderr, "\n");

#define REF_NUCL(pos) aligned_ref[(pos)]
#define READ_NUCL(read,pos) _p[(1 + read) * 2048 + (pos)]

  /* Generate gapped alignment */
  unsigned int p_len = 0;
  int *_p;
  _p = (int *) malloc ((1 + na) * 2048 * 4);
  unsigned int aligned_ref[2048];
  int ref_pos[1024];
  int ref_p = 0, last_ref_p = UNKNOWN;
  int read_p[1024], last_read_p[1024];
  for (i = 0; i < na; i++) {
    read_p[i] = a[i][0];
    last_read_p[i] = UNKNOWN;
  }
  while (ref_p < ref_seq->len) {
    /* Write alignment */
    /* Ref */
    if ((last_ref_p < 0) || (ref_p > last_ref_p)) {
      _p[p_len] = ref_seq->pos[ref_p].nucl;
      aligned_ref[p_len] = ref_seq->pos[ref_p].nucl;
      ref_pos[p_len] = ref_start + ref_p;
      last_ref_p = ref_p;
    } else {
      _p[p_len] = GAP;
      aligned_ref[p_len] = GAP;
      ref_pos[p_len] = ref_start + ref_p;
    }
    /* Reads */
    for (i = 0; i < na; i++) {
      if ((read_p[i] >= 0) && ((last_read_p[i] < 0) || (read_p[i] > last_read_p[i]))) {
        _p[(1 + i) * 2048 + p_len] = a_reads[i]->nseq->pos[read_p[i]].nucl;
        last_read_p[i] = read_p[i];
      } else if (read_p[i] >= 0) {
        _p[(1 + i) * 2048 + p_len] = GAP;
      } else {
        _p[(1 + i) * 2048 + p_len] = NONE;
      }
    }
    /* Advance */
    int rgap = 1;
    if (ref_p < (ref_seq->len - 1)) {
      int next_ref_p = ref_p + 1;
      for (i = 0; i < na; i++) {
        int next_read_p = a[i][next_ref_p];
        if ((read_p[i] >= 0) && (next_read_p >= 0)) {
          int gap = next_read_p - read_p[i];
          if (gap > rgap) rgap = gap;
        }
      }
    }
    if (ref_p < (ref_seq->len - 1)) {
      int next_ref_p = ref_p + 1;
      for (i = 0; i < na; i++) {
        int next_read_p = a[i][next_ref_p];
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

#define COUNT(p,n) counts[(p) * (GAP + 1) + (n)]

  /* Calculate totals */
  unsigned int *totals, *counts;
  totals = (unsigned int *) malloc (p_len * 4);
  memset (totals, 0, p_len * 4);
  counts = (unsigned int *) malloc (p_len * (GAP + 1) * 4);
  memset (counts, 0, p_len * (GAP + 1) * 4);
  for (i = 0; i < p_len; i++) {
    unsigned int j;
    for (j = 0; j < na; j++) {
      if (READ_NUCL(j, i) <= GAP) {
        int nucl = READ_NUCL(j, i);
        counts[i * (GAP + 1) + nucl] += 1;
        totals[i] += 1;
      }
    }
  }

  /* Tag all reads by divergent positions */
  unsigned int n_divergent = 0;
  for (i = 0; i < p_len; i++) {
    unsigned int j;
    unsigned int diverges = 0;
    for (j = 0; j <= GAP; j++) {
      if (j == REF_NUCL(i)) continue;
      if (j == N) continue;
      if (counts[i * (GAP + 1) + j] >= 2) diverges = 1;
    }
    if (diverges) {
      unsigned int snv, known = 0;
      if (debug) fprintf (stderr, "Divergent position: %u\n", ref_pos[i]);
      snv = lookup_snv (snvs, n_snvs, ref_chr, ref_start + i);
      if ((snv < n_snvs) && (snvs[snv].chr == ref_chr) && (snvs[snv].pos == (ref_start + i))) {
        fprintf (stderr, "Known SNV %s (%c/%c)\n", snvs[snv].id, n2c[snvs[snv].ref_allele], n2c[snvs[snv].alt_allele]);
        known = 1;
      } else {
        fprintf (stderr, "Potential DeNovo\n");
      }
      for (j = 0; j < na; j++) {
        unsigned int ref = REF_NUCL(i);
        unsigned int nucl = READ_NUCL(j,i);
        unsigned int mask = 7;
        /* Do not count single nucleotides */
        if ((nucl <= GAP) && (COUNT(i, nucl) < 2)) mask = 0;
        /* N-s are counted same as reference */
        if (nucl == N) nucl = ref;
        /* Uncovered positions are counted as reference but masked */
        if (nucl > GAP) {
          nucl = ref;
          mask = 0;
        }
        a_reads[j]->unknown = a_reads[j]->unknown << 3;
        if (!known || ((nucl != snvs[snv].ref_allele) && (nucl != snvs[snv].alt_allele))) {
          a_reads[j]->unknown |= 7;
        }
        /* Make 000 reference variant */
        nucl = nucl ^ ref;
        a_reads[j]->tag = (a_reads[j]->tag << 3) | nucl;
        a_reads[j]->mask = (a_reads[j]->mask << 3) | mask;
      }
      n_divergent += 1;
    }
  }

  unsigned int read_group[1024];
  unsigned long long group_tags[1024];
  unsigned long long group_masks[1024];
  unsigned int group_sizes[1024];
  unsigned int n_groups = na;
  for (i = 0; i < na; i++) {
    read_group[i] = i;
    group_tags[i] = a_reads[i]->tag & a_reads[i]->mask;
    group_masks[i] = a_reads[i]->mask;
    group_sizes[i] = 1;
  }
  if (debug > 1) {
    for (i = 0; i < n_groups; i++) fprintf (stderr, "%llu\t", group_tags[i]);
    fprintf (stderr, "\n");
    for (i = 0; i < n_groups; i++) fprintf (stderr, "%llu\t", group_masks[i]);
    fprintf (stderr, "\n");
  }
  unsigned int *is_compat;
  unsigned int *n_common;
  is_compat = (unsigned int *) malloc (n_groups * n_groups * 4);
  memset (is_compat, 0, n_groups * n_groups * 4);
  n_common = (unsigned int *) malloc (n_groups * n_groups * 4);
  memset (n_common, 0, n_groups * n_groups * 4);
  while (n_groups > 1) {
    for (i = 0; i < n_groups; i++) {
      for (j = 0; j < n_groups; j++) {
        if (j == i) {
          is_compat[i * na + j] = 0;
          n_common[i * na + j] = 0;
          continue;
        }
        unsigned long long common = group_masks[i] & group_masks[j];
        if ((group_tags[i] & common) == (group_tags[j] & common)) {
          /* Groups are compatible */
          is_compat[i * na + j] = 1;
        } else {
          is_compat[i * na + j] = 0;
        }
        n_common[i * na + j] = 0;
        while (common != 0) {
          if (common & 7) n_common[i * na + j] += 1;
          common = common >> 3;
        }
      }
    }
    unsigned int max_i = 0;
    unsigned int max_j = 0;
    unsigned int found = 0;
    for (i = 0; i < n_groups; i++) {
      for (j = i + 1; j < n_groups; j++) {
        if (is_compat[i * na + j]) {
          if (!found) {
            max_i = i;
            max_j = j;
            found = 1;
          } else {
            if (n_common[i * na + j] > n_common[max_i * na + max_j]) {
              max_i = i;
              max_j = j;
            } else if (n_common[i * na + j] == n_common[max_i * na + max_j]) {
              if ((group_sizes[i] + group_sizes[j]) > (group_sizes[max_i] + group_sizes[max_j])) {
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
    if (debug > 0) fprintf (stderr, "Merging groups %u (size %u) and %u (size %u) (common %u): %llu %llu %llu %llu -> ", max_i, group_sizes[max_i], max_j, group_sizes[max_j], n_common[max_i * na + max_j], group_tags[max_i], group_masks[max_i], group_tags[max_j], group_masks[max_j]);
    group_tags[max_i] = (group_tags[max_i] & group_masks[max_i]) | (group_tags[max_j] & group_masks[max_j]);
    group_masks[max_i] = group_masks[max_i] | group_masks[max_j];
    group_sizes[max_i] += group_sizes[max_j];
    if (debug > 0) fprintf (stderr, "%llu %llu\n", group_tags[max_i], group_masks[max_i]);
    for (i = 0; i < na; i++) if (read_group[i] == max_j) read_group[i] = max_i;
    n_groups -= 1;
    group_tags[max_j] = group_tags[n_groups];
    group_masks[max_j] = group_masks[n_groups];
    group_sizes[max_j] = group_sizes[n_groups];
    for (i = 0; i < na; i++) if (read_group[i] == n_groups) read_group[i] = max_j;
  }
  fprintf (stderr, "Num remaining groups: %u\n", n_groups);
  if (debug > 0) {
    for (i = 0; i < n_groups; i++) fprintf (stderr, "%llu\t", group_tags[i]);
    fprintf (stderr, "\n");
    for (i = 0; i < n_groups; i++) fprintf (stderr, "%llu\t", group_masks[i]);
    fprintf (stderr, "\n");
  }
  free (is_compat);
  free (n_common);

  if (debug > 1) {
    fprintf (stderr, "Read groups:");
    for (i = 0; i < na; i++) fprintf (stderr, " %u:%u", i, read_group[i]);
    fprintf (stderr, "\n");
  }

  /* Create group consensus */
  unsigned int *g_cons = (unsigned int *) malloc (n_groups * p_len * 4);
  unsigned int g_divergent[32] = { 0 };
  for (i = 0; i < p_len; i++) {
    for (j = 0; j < n_groups; j++) {
      unsigned int c[10] = { 0 };
      for (k = 0; k < na; k++) {
        if (read_group[k] == j) c[READ_NUCL(k,i)] += 1;
      }
      unsigned int best = REF_NUCL(i);
      for (k = 0; k <= GAP; k++) {
        if (k == N) continue;
        if ((COUNT(i,k) > 1) && (c[k] > c[best])) best = k;
      }
      g_cons[j * p_len + i] = best;
      if (best != REF_NUCL(i)) {
        unsigned int snv;
        if (debug) fprintf (stderr, "Divergent position in group %u %u:%u\n", j, ref_chr, ref_pos[i]);
        snv = lookup_snv (snvs, n_snvs, ref_chr, ref_start + i);
        if ((snv < n_snvs) && (snvs[snv].chr == ref_chr) && (snvs[snv].pos == (ref_start + i))) {
          fprintf (stderr, "Known SNV (%c/%c)\n", n2c[snvs[snv].ref_allele], n2c[snvs[snv].alt_allele]);
        } else {
          fprintf (stderr, "Potential DeNovo\n");
          g_divergent[j] += 1;
        }
      }
    }
  }

  if (debug > 0) {
    for (i = 0; i < n_groups; i++) {
      unsigned int j;
      fprintf (stderr, "\nGroup %u size %u divergent %u\n", i, group_sizes[i], g_divergent[i]);
      for (j = 0; j < p_len; j++) fprintf (stderr, "%c", n2c[g_cons[i * p_len + j]]);
      fprintf (stderr, "\n");
      for (j = 0; j < na; j++) {
        if (read_group[j] == i) {
          fprintf (stderr, "%s\n", a_reads[j]->name);
        }
      }
    }
  }

#define MIN_GROUP_SIZE 8
#define MIN_GROUP_RSIZE 0.25
#define MAX_DIVERGENT 2
#define MAX_DIVERGENT_DELTA 1
  
  /* Discard groups with too many divergent positions */
  unsigned int min_div = g_divergent[0];
  for (i = 0; i < n_groups; i++) if (g_divergent[i] < min_div) min_div = g_divergent[i];
  unsigned int incl_group[32];
  for (i = 0; i < n_groups; i++) {
    incl_group[i] = 1;
    if (group_sizes[i] < MIN_GROUP_SIZE) {
      incl_group[i] = 0;
      fprintf (stderr, "Discarded group %u: size too small (%u < %u)\n", i, group_sizes[i], MIN_GROUP_SIZE);
      continue;
    }
    if (g_divergent[i] > MAX_DIVERGENT) {
      incl_group[i] = 0;
      fprintf (stderr, "Discarded group %u: too big divergence (%u > %u)\n", i, g_divergent[i], MAX_DIVERGENT);
      continue;
    }
    if (g_divergent[i] > (min_div + MAX_DIVERGENT_DELTA)) {
      incl_group[i] = 0;
      fprintf (stderr, "Discarded group %u: too big relative divergence (%u > %u)\n", i, g_divergent[i], min_div + MAX_DIVERGENT_DELTA);
      continue;
    }
    if (group_sizes[i] < (group_sizes[0] * MIN_GROUP_RSIZE)) {
      incl_group[i] = 0;
      fprintf (stderr, "Discarded group %u: relative size too small (%.2f < %.2f)\n", i, (double) group_sizes[i] / group_sizes[0], MIN_GROUP_RSIZE);
      continue;
    }
  }
  
  /* Recalculate totals */
  memset (totals, 0, p_len * 4);
  memset (counts, 0, p_len * (GAP + 1) * 4);
  for (i = 0; i < p_len; i++) {
    unsigned int j;
    for (j = 0; j < na; j++) {
      if (!incl_group[read_group[j]]) continue;
      if (READ_NUCL(j, i) <= GAP) {
        int nucl = READ_NUCL(j, i);
        counts[i * (GAP + 1) + nucl] += 1;
        totals[i] += 1;
      }
    }
  }
  
  /* Output alignment */
  fprintf (stdout, "CHR\tPOS      \tREF\tTOTAL\tA\tC\tG\tT\tN\tGAP\tCALL\tPROB\tRPROB");
  if (debug) {
    fprintf (stdout, "\t ");
    for (i = 0; i < n_groups; i++) {
      unsigned int j;
      if (!incl_group[i]) continue;
      fprintf (stdout, "       ");
      for (j = 0; j < na; j++) {
        if (read_group[j] == i) fprintf (stdout, "%c", 'A' + i);
      }
    }
  }
  fprintf (stdout, "\n");
  for (i = 0; i < p_len; i++) {
    if (REF_NUCL(i) != GAP) {
      fprintf (stdout, "%s\t%9u\t%c", chr_names[ref_chr], ref_pos[i], n2c[REF_NUCL(i)]);
    } else {
      fprintf (stdout, "\t         \t-");
    }
    /* Totals */
    fprintf (stdout, "\t%u", totals[i]);
    for (j = A; j <= GAP; j++) {
      fprintf (stdout, "\t%u", counts[i * (GAP + 1) + j]);
    }
    /* Call */
    if (totals[i] < 10) {
      fprintf (stdout, "\tNC\t\t");
    } else {
      unsigned int n1, n2, best_n1 = A, best_n2 = A;
      float best_prob = 0;
      float sum_probs = 0;
      for (n1 = A; n1 <= GAP; n1++) {
        if (n1 == N) continue;
        unsigned int c1 = counts[i * (GAP + 1) + n1];
        if (c1 < 2) continue;
        for (n2 = n1; n2 <= GAP; n2++) {
          float prob;
          if (n2 == N) continue;
          unsigned int c2 = counts[i * (GAP + 1) + n2];
          if (c2 < 2) continue;
          if (n2 == n1) {
            prob = gt1_prob (c1, totals[i] - counts[i * (GAP + 1) + N]);
          } else {
            prob = gt2_prob (c1, c2, totals[i] - counts[i * (GAP + 1) + N]);
          }
          if (prob > best_prob) {
            best_n1 = n1;
            best_n2 = n2;
            best_prob = prob;
          }
          sum_probs += prob;
        }
      }
      if (!sum_probs) {
        best_prob = 0;
        sum_probs = 1;
      }
      fprintf (stdout, "\t%c%c%c", n2c[best_n1], n2c[best_n2], ((best_n1 == REF_NUCL(i)) && (best_n2 == REF_NUCL(i))) ? ' ' : '*');
      fprintf (stdout, "\t%.3f\t%.3f", best_prob, best_prob / sum_probs);
    }
    if (debug) {
      /* Aligned reads */
      fprintf (stdout, "\t%c", n2c[REF_NUCL(i)]);
      for (j = 0; j < n_groups; j++) {
        if (!incl_group[j]) continue;
        fprintf (stdout, "  [%c%c] ", n2c[g_cons[j * p_len + i]], (g_cons[j * p_len + i] == REF_NUCL(i)) ? ' ' : '*');
        for (k = 0; k < na; k++) {
          if (read_group[k] == j) fprintf (stdout, "%c", n2c[READ_NUCL(k,i)]);
        }
      }
    }
    fprintf (stdout, "\n");
  }
  
  free (totals);
  free (counts);
  free (g_cons);
  free (_p);
  for (i = 0; i < na; i++) free (a[i]);
  n_seq_delete (ref_seq);
  for (i = 0; i < nreads; i++) gasm_read_delete (reads[i]);
}

static void
assemble_by_matrix (KMerDB *db, KMerDB *gdb, const char *seqs[], unsigned int nreads, const char *ref_chr, unsigned int ref_start, unsigned int ref_end, const char *ref)
{
  NMatrix *mat;
  unsigned int i;

  /* Create matrix */
  mat = n_matrix_new (nreads, seqs, WORDLEN);
  if (debug) {
    fprintf (stderr, "Kmers %u unique %u\n", mat->n_kmers, mat->n_unique_kmers);
  }
  if (print_kmers) {
    char c[32] = { 0 };
    for (i = 0; i < mat->n_unique_kmers; i++) {
      word2string (c, mat->unique_kmers[i].value, WORDLEN);
      fprintf (stdout, "%u %s %u\n", i, c, mat->unique_kmers[i].count);
    }
  }
  if (analyze_kmers) {
    print_kmer_statistics (mat);
  }

  NCell *alignment[2048];
  unsigned int alen = alignment_from_matrix (alignment, 2048, db, gdb, mat, min_coverage);
  if (alen < 10) {
    fprintf (stderr, "Error: cannot create alignment (len = %d) - too low coverage?\n", alen);
    return;
  }
  unsigned int a[2048];
  unsigned int b[2048];
  unsigned int a_pos[2048];
  unsigned int b_pos[2048];
  unsigned int r[2048] = { 0 };
  for (i = 0; i < alen; i++) a[i] = alignment[i]->consensus;
  NSeq *b_seq = n_seq_new (ref, WORDLEN);
  for (i = 0; i <= b_seq->len; i++) b[i] = b_seq->pos[i].nucl;
  unsigned int n_ref_align = smith_waterman (a_pos, b_pos, a, alen, b, b_seq->len);
  fprintf (stderr, "Length of reference alignment %u\n", n_ref_align);  
  for (i = 0; i < n_ref_align; i++) {
    if (a_pos[i] > 10000) continue;
    if (b_pos[i] > 10000) continue;
    r[a_pos[i]] = b_seq->pos[b_pos[i]].nucl;
    alignment[a_pos[i]]->ref_pos = ref_start + b_pos[i];
  }

  fprintf (stderr, "Length before final alignment %u\n", alen);  
  /* Compose alignment */
  int n_p[1024];
  for (i = 0; i < 1024; i++) n_p[i] = BEFORE;
  int file_idx = -1;
  unsigned int pos[1024] = { 0 };
  unsigned int ac[1024];
  unsigned int rc[1024];
  unsigned int nc[6 * 1024] = { 0 };
  unsigned int ac_len = 0;
  i = 0;
  unsigned int next_rp = 0;
  while (i < alen) {
    /* Check whether any sequence needs additional nucleotides */
    NCell *cell = alignment[i];
    /* fprintf (stdout, "%u\t%u\n", cell->file_idx, cell->ref_pos); */
    unsigned int has_gap = 0, has_rgap = 0;
    unsigned int j;
    for (j = 0; j < mat->n_seqs; j++) {
      int s_p = cell->links[j].pos;
      if ((s_p >= 0) && (n_p[j] >= 0) && (n_p[j] < s_p)) {
        has_gap = 1;
        break;
      }
    }
    if (next_rp && (cell->ref_pos > 0) && (cell->ref_pos > next_rp)) {
      fprintf (stderr, "Has gap %u %u\n", cell->ref_pos, next_rp);
      has_rgap = 1;
    }
    if (has_gap) {
      ac[ac_len] = GAP;
      rc[ac_len] = GAP;
      for (j = 0; j < mat->n_seqs; j++) {
        int s_p = cell->links[j].pos;
        if (n_p[j] < 0) {
        } else if ((n_p[j] >= 0) && (n_p[j] < s_p)) {
          nc[6 * ac_len + mat->seqs[j]->pos[n_p[j]].nucl] += 1;
          n_p[j] += 1;
          if (n_p[j] >= mat->seqs[j]->len) n_p[j] = -1;
        } else {
          nc[6 * ac_len + GAP] += 1;
        }
      }
      ac_len += 1;
      if (ac_len >= 1024) break;
    } else if (has_rgap) {
      int b_pos = (int) next_rp - (int) ref_start;
      pos[ac_len] = next_rp;
      ac[ac_len] = GAP;
      rc[ac_len] = ((b_pos >= 0) && (b_pos < b_seq->len)) ? b_seq->pos[b_pos].nucl : NONE;
      next_rp += 1;
      ac_len += 1;
      if (ac_len >= 1024) break;
    } else {
      if ((int) cell->file_idx > file_idx) file_idx = (int) cell->file_idx;
      if (cell->ref_pos > 0) pos[ac_len] = cell->ref_pos;
      ac[ac_len] = cell->consensus;
      rc[ac_len] = r[i];
      if (cell->ref_pos > 0) next_rp = cell->ref_pos + 1;
      for (j = 0; j < mat->n_seqs; j++) {
        int s_p = cell->links[j].pos;
        if (s_p >= 0) {
          if ((n_p[j] < 0) || (s_p == n_p[j])) {
            nc[6 * ac_len + mat->seqs[j]->pos[s_p].nucl] += 1;
            n_p[j] = s_p + 1;
            if (n_p[j] >= mat->seqs[j]->len) n_p[j] = -1;
          } else {
            nc[6 * ac_len + GAP] += 1;
          }
        } else {
        }
      }
      ac_len += 1;
      if (ac_len >= 1024) break;
      i += 1;
    }
  }
  fprintf (stderr, "Final alignment length %u\n", ac_len);

  fprintf (stdout, "CHR\tPOS\tREF\tCONS\tTOTAL\tA\tC\tG\tT\tN\tGAP\tALLELES\n");
  for (i = 0; i < ac_len; i++) {
    unsigned int j, sum;
    sum = 0;
    for (j = 0; j < 6; j++) sum += nc[6 * i + j];
    /*if (sum < 6) continue;*/
    if ((sum > 0) && (nc[6 * i + GAP] >= (2 * sum / 3))) continue;
    if (pos[i] > 0) {
      if (pos[i] < ref_start) continue;
      if (pos[i] >= ref_end) break;
      fprintf (stdout, "%s\t%u\t", ref_chr, pos[i] + 1);
    } else {
      fprintf (stdout, "\t\t");
    }
    fprintf (stdout, "%c\t%c\t%u", n2c[rc[i]], n2c[ac[i]], sum);
    unsigned int s[6] = {0};
    for (j = 0; j < 6; j++) {
      if ((sum >= 6) && (nc[6 * i + j] >= 3)) s[j] = 1;
      fprintf (stdout, "\t%d", nc[6 * i + j]);
    }
    fprintf (stdout, "\t");
    for (j = 0; j < 6; j++) {
      if (s[j]) fprintf (stdout, "%c", n2c[j]);
    }
    fprintf (stdout, "\n");
  }
  
  /* Free seqs */
  n_matrix_delete (mat);
}

/* Create all matrix cells linked by this kmer */

static void
link_seqs_by_kmer (NMatrix *mat, unsigned long long word)
{
  unsigned int start, end, nkmers, i;
  for (start = 0; start < mat->n_kmers; start++) if (mat->kmers[start].value == word) break;
  for (end = start + 1; end < mat->n_kmers; end++) if (mat->kmers[end].value != word) break;
  nkmers = end - start;
  for (i = 0; i < WORDLEN; i++) {
    unsigned int seqs[MAX_READS];
    unsigned int positions[MAX_READS];
    unsigned int j;
    memset (positions, 0, sizeof (positions));
    for (j = 0; j < nkmers; j++) {
      seqs[j] = mat->kmers[start + j].seq;
      positions[j] = mat->kmers[start + j].pos + i;
    }
    n_matrix_link_sequences (mat, seqs, positions, nkmers);
  }
}

/* Ensure that all sequences are linked at next positions, forking alignment if needed */

static void
link_to_next_recursive (NMatrix *mat, NCell *cell) {
  unsigned int i;
  NCell *next_cell[2048] = { 0 };
  unsigned int n_unique_cells = 0;
  unsigned int nucl_count[6] = { 0 };
  unsigned int total_count = 0;
  for (i = 0; i < mat->n_seqs; i++) {
    int next_pos;
    NCell *nc;
    /* Skip if this sequence is not linked */
    if (cell->links[i].pos < 0) continue;
    next_pos = cell->links[i].pos + 1;
    /* Skip if this sequence ends here */
    if (next_pos >= mat->seqs[i]->len) continue;
    /* Add next cell if present */
    for (nc = mat->seqs[i]->pos[next_pos].cells; nc; nc = nc->links[i].next) {
      next_cell[n_unique_cells++] = nc;
    }
    /* Count nucleotides at next position */
    nucl_count[i] += mat->seqs[i]->pos[next_pos].nucl;
    total_count += 1;
  }
  /* Find unique next cells */
  for (i = 0; i < n_unique_cells; i++) {
    unsigned int j;
    j = i + 1;
    while (j < n_unique_cells) {
      if (next_cell[j] == next_cell[i]) {
        n_unique_cells -= 1;
        next_cell[j] = next_cell[n_unique_cells];
      } else {
        j += 1;
      }
    }
  }
  /* Now, for all sequences not linked, find if any next cell fits */
  for (i = 0; i < mat->n_seqs; i++) {
    int next_pos;
    if (cell->links[i].pos < 0) continue;
    /* Skip if this sequence is not linked */
    if (cell->links[i].pos < 0) continue;
    next_pos = cell->links[i].pos + 1;
    /* Skip if this sequence ends here */
    if (next_pos >= mat->seqs[i]->len) continue;
    /* Skip if this sequence is already linked */
    if (mat->seqs[i]->pos[next_pos].cells) continue;
    
  }
}

static unsigned int
alignment_from_matrix (NCell *alignment[], unsigned int alignment_size, KMerDB *db, KMerDB *gdb, NMatrix *mat, unsigned int min_cov)
{
  unsigned int i;
  NCell *cell, *best_cell;
  /* Fill matrix */
  fprintf (stderr, "Aligning matrix\n");
  /* Iterate over all unique kmers */
  for (i = 0; i < mat->n_unique_kmers; i++) {
    link_seqs_by_kmer (mat, mat->unique_kmers[i].value);
  }
  best_cell = n_matrix_calculate_scores (mat);
  if (debug) fprintf (stderr, "Best cell score: %d\n", best_cell->score);

  /* Move forward from best cell linking sequences if needed */

  /* Now we have linked all reads by kmers */
  fprintf (stderr, "Linking remaining parts of reads\n");
  /* Link remaining parts of reads (N-s) */
  for (i = 0; i < mat->n_seqs; i++) {
    NSeq *seq = mat->seqs[i];
    unsigned int j;
    for (j = 0; j < seq->len; j++) {
      NCell *cell = seq->pos[j].cells;
      if (!cell) {
        cell = n_matrix_new_cell (mat);
        n_matrix_link_cell (mat, cell, i, j);
      }
    }
  }

  fprintf (stderr, "Chaining cells with biggest neighbours\n");
  /* Chain each cell with neighbour with largest number of per-seq links */
  for (cell = mat->cells; cell; cell = cell->next_allocated) {
    NCell *prev = NULL, *next = NULL;
    unsigned int max_prev, max_next, sum_prev, sum_next;
    /* Reset link counter */
    for (i = 0; i < mat->n_seqs; i++) {
      if (cell->links[i].pos < 0) continue;
      NCell *c = n_matrix_get_seq_cell (mat, i, (int) cell->links[i].pos - 1);
      if (c) c->count = 0;
      c = n_matrix_get_seq_cell (mat, i, (int) cell->links[i].pos + 1);
      if (c) c->count = 0;
    }
    /* Count sequences for links */
    sum_prev = sum_next = 0;
    for (i = 0; i < mat->n_seqs; i++) {
      if (cell->links[i].pos < 0) continue;
      /* fixme: If next and prev point to the same cell it messes up */
      NCell *c = n_matrix_get_seq_cell (mat, i, (int) cell->links[i].pos - 1);
      if (c) {
        c->count += 1;
        sum_prev += 1;
      }
      c = n_matrix_get_seq_cell (mat, i, (int) cell->links[i].pos + 1);
      if (c) {
        c->count += 1;
        sum_next += 1;
      }
    }
    max_prev = max_next = 0;
    /* Find prev and next links with max sequences */
    for (i = 0; i < mat->n_seqs; i++) {
      NCell *c = n_matrix_get_seq_cell (mat, i, (int) cell->links[i].pos - 1);
      if (c && (c->count > max_prev)) {
        prev = c;
        max_prev = c->count;
      }
      c = n_matrix_get_seq_cell (mat, i, (int) cell->links[i].pos + 1);
      if (c && (c->count > max_next)) {
        next = c;
        max_next = c->count;
      }
    }
    if (debug > 1) fprintf (stderr, "Max prev %u max next %u\n", max_prev, max_next);
    /* Analyze */
    for (i = 0; i < mat->n_seqs; i++) {
      NCell *cc = n_matrix_get_seq_cell (mat, i, (int) cell->links[i].pos - 1);
      if (cc && (cc != prev)) {
        unsigned int j;
        fprintf (stderr, "\nPrev inconsistency for sequence %u pos %d\n", i, cell->links[i].pos);
        fprintf (stderr, "Best link pos %d count %u\n", prev->links[i].pos, prev->count);
        fprintf (stderr, "Sequence pos  %d count %u\n", cc->links[i].pos, cc->count);
        fprintf (stderr, "Total links %u\n", sum_prev);
        for (j = 0; j < mat->n_seqs; j++) {
          cc = n_matrix_get_seq_cell (mat, j, (int) cell->links[j].pos - 1);
          if (cc) {
            const char *c = "ACGTN-*";
            const char *d = "?><";
            int p0 = cc->links[j].pos;
            char n0 = (p0 >= 0) ? c[mat->seqs[j]->pos[p0].nucl] : d[p0 + 3];
            int p1 = cell->links[j].pos;
            char n1 = (p1 >= 0) ? c[mat->seqs[j]->pos[p1].nucl] : d[p1 + 3];
            fprintf (stderr, " %u:%d/%d:%c/%c", j, p0, p1, n0, n1);
          }
        }
        fprintf (stderr, "\n");
        break;
      }
      cc = n_matrix_get_seq_cell (mat, i, (int) cell->links[i].pos + 1);
      if (cc && (cc != next)) {
        unsigned int j;
        fprintf (stderr, "\nNext inconsistency for sequence %u pos %d\n", i, cell->links[i].pos);
        fprintf (stderr, "Best link pos %d count %u\n", next->links[i].pos, next->count);
        fprintf (stderr, "Sequence pos  %d count %u\n", cc->links[i].pos, cc->count);
        fprintf (stderr, "Total links %u\n", sum_next);
        for (j = 0; j < mat->n_seqs; j++) {
          cc = n_matrix_get_seq_cell (mat, j, (int) cell->links[j].pos + 1);
          if (cc) {
            const char *c = "ACGTN-*";
            const char *d = "?><";
            int p0 = cc->links[j].pos;
            char n0 = (p0 >= 0) ? c[mat->seqs[j]->pos[p0].nucl] : d[p0 + 3];
            int p1 = cell->links[j].pos;
            char n1 = (p1 >= 0) ? c[mat->seqs[j]->pos[p1].nucl] : d[p1 + 3];
            fprintf (stderr, " %u:%d/%d:%c/%c", j, p0, p1, n0, n1);
          }
        }
        fprintf (stderr, "\n");
        break;
      }
    }
    if (debug > 2) fprintf (stderr, "Chained cells\n");
    if (prev) cell->prev = prev;
    if (next) cell->next = next;
  }

  /* Get longest alignment */
  fprintf (stderr, "Finding longest alignment\n");
  unsigned int max = 0;
  NCell *max_cell = NULL;
  for (cell = mat->cells; cell; cell = cell->next_allocated) {
    unsigned int count = 0;
    for (i = 0; i < mat->n_seqs; i++) {
      NCell *pc = n_matrix_get_seq_cell (mat, i, (int) cell->links[i].pos - 1);
      NCell *nc = n_matrix_get_seq_cell (mat, i, (int) cell->links[i].pos + 1);
      if (pc && nc) count += 1;
    }
    if (count > max) {
      max_cell = cell;
      max = count;
    }
  }
  fprintf (stderr, "Finding first cell\n");
  for (cell = mat->cells; cell; cell = cell->next_allocated) cell->count = 0;
  NCell *first_cell = max_cell;
  while (first_cell->prev) {
    first_cell->count = 1;
    if (first_cell->prev && first_cell->prev->count) {
      /* Remove all prev links */
      for (i = 0; i < mat->n_seqs; i++) {
        if (first_cell->links[i].prev) {
          n_matrix_unlink_cell (mat, first_cell->links[i].prev, i);
        }
      }
      first_cell->prev = NULL;
      break;
    }
    first_cell = first_cell->prev;
  }
  for (cell = mat->cells; cell; cell = cell->next_allocated) cell->count = 0;
  unsigned int alen = 0;
  cell = first_cell;
  while (cell) {
    cell->count = 1;
    alen += 1;
    if (cell->next && cell->next->count) {
      /* Remove all next links */
      for (i = 0; i < mat->n_seqs; i++) {
        if (cell->links[i].next) {
          n_matrix_unlink_cell (mat, cell->links[i].next, i);
        }
      }
      cell->next = NULL;
    }
    cell = cell->next;
  }
  if (debug) fprintf (stderr, "Initial alignment length %u\n", alen);
  
  /* Compact alignment */
  for (cell = mat->cells; cell; cell = cell->next_allocated) cell->count = 0;
  cell = first_cell;
  while (cell->next) {
    cell->count = 1;
    if (cell->next && cell->next->count) break;
    if (n_matrix_can_merge_cells (mat, cell, cell->next)) {
      NCell *next = cell->next->next;
      n_matrix_merge_cells (mat, cell, cell->next);
      cell->next = next;
      alen -= 1;
    } else {
      cell = cell->next;
    }
  }
  if (debug) fprintf (stderr, "Trimmed alignment length %u\n", alen);
  
  /* Create alignment */
  for (cell = mat->cells; cell; cell = cell->next_allocated) cell->count = 0;
  unsigned int apos = 0;
  for (cell = first_cell; cell; cell = cell->next) {
    cell->count = 1;
    alignment[apos++] = cell;
    if (cell->next && cell->next->count) break;
  }
  alen = apos;

  /* Calculate consensus nucleotides */
  for (i = 0; i < alen; i++) {
    unsigned int n_counts[4] = { 0 };
    unsigned int j, max;
    NCell *cell = alignment[i];
    for (j = 0; j < mat->n_seqs; j++) {
      if (cell->links[j].pos >= 0) {
        unsigned int n = mat->seqs[j]->pos[cell->links[j].pos].nucl;
        if (n < 4) n_counts[n] += 1;
      }
    }
    max = A;
    for (j = 1; j < 4; j++) if (n_counts[j] > n_counts[max]) max = j;
    cell->consensus = max;
  }

  /* Align cells back to consensus */
  for (i = 0; i < mat->n_seqs; i++) {
    align_seq (mat, alignment, alen, i);
  }

  /* Try to merge cells in alignment */
  for (cell = mat->cells; cell; cell = cell->next_allocated) cell->count = 0;
  for (i = 0; i < alen; i++) {
    unsigned int j;
    NCell *cell = alignment[i];
    for (j = 0; j < mat->n_seqs; j++) {
      unsigned int k;
      if (cell->links[j].pos < 0) continue;
      for (k = 0; k < i; k++) {
        int s_p = (int) cell->links[j].pos - (int) (i - k);
        /* Check whether this cell has given sequence linked */
        if (alignment[k]->links[j].pos != UNKNOWN) continue;
        /* Need link for this sequence */
        if (s_p < 0) {
          alignment[k]->links[j].pos = BEFORE;
        } else {
          NCell *cc;
          for (cc = mat->seqs[j]->pos[s_p].cells; cc; cc = cc->links[j].next) {
            unsigned int tag = i + 1;
            if (cc->count == tag) {
              /* Cycle */
              break;
            }
            cc->count = tag;
            if (n_matrix_can_merge_cells (mat, alignment[k], cc)) {
              if (debug > 1) fprintf (stderr, "Merging cells: alignment pos %u seq %u pos %u\n", i, j, s_p);
              n_matrix_merge_cells (mat, alignment[k], cc);
              break;
            } else {
              if (debug > 1) fprintf (stderr, "Cannot merge: alignment pos %u seq %u pos %u\n", i, j, s_p);
            }
          }
        }
      }
      for (k = i + 1; k < alen; k++) {
        int s_p = (int) cell->links[j].pos + (int) (k - i);
        /* Check whether this cell has given sequence linked */
        if (alignment[k]->links[j].pos != UNKNOWN) continue;
        /* Need link for this sequence */
        if (s_p >= mat->seqs[j]->len) {
          alignment[k]->links[j].pos = AFTER;
        } else {
          NCell *cc;
          for (cc = mat->seqs[j]->pos[s_p].cells; cc; cc = cc->links[j].next) {
            unsigned int tag = i + 1000000;
            if (cc->count == tag) {
              /* Cycle */
              break;
            }
            cc->count = tag;
            if (n_matrix_can_merge_cells (mat, alignment[k], cc)) {
              if (debug > 1) fprintf (stderr, "Merging cells: alignment pos %u seq %u pos %u\n", i, j, s_p);
              n_matrix_merge_cells (mat, alignment[k], cc);
              break;
            } else {
              if (debug > 1) fprintf (stderr, "Cannot merge: alignment pos %u seq %u pos %u\n", i, j, s_p);
            }
          }
        }
      }
    }
    cell->count = 0;
    for (j = 0; j < mat->n_seqs; j++) {
      if (cell->links[j].pos >= 0) cell->count += 1;
    }
  }
  /* Trim alignment */
  int a_start = 0;
  while (a_start < alen) {
    if (alignment[a_start]->count >= min_cov) break;
    a_start += 1;
  }
  int a_end = alen - 1;
  while (a_end >= 0) {
    if (alignment[a_end]->count >= min_cov) break;
    a_end -= 1;
  }
  if (a_start >= a_end) return 0;
  if (a_start > 0) {
    int j;
    for (j = a_start; j < a_end; j++) alignment[j - a_start] = alignment[j];
    alen = a_end - a_start + 1;
  }
  return alen;
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
  last_a = -1;
  last_b = -1;
  for (i = 0; i < len; i++) {
    if ((int) a_pos[i] > last_a) {
      fprintf (ofs, "%c", n2c[a->pos[a_pos[i]].nucl]);
    } else {
      fprintf (ofs, "-");
    }
    last_a = a_pos[i];
    last_b = b_pos[i];
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
  last_a = -1;
  last_b = -1;
  for (i = 0; i < len; i++) {
    if (((int) a_pos[i] > last_a) && ((int) b_pos[i] > last_b) && (a->pos[a_pos[i]].nucl != N) && (b->pos[b_pos[i]].nucl != N) && (a->pos[a_pos[i]].nucl == b->pos[b_pos[i]].nucl)) {
      fprintf (ofs, "|");
    } else {
      fprintf (ofs, " ");
    }
    last_a = a_pos[i];
    last_b = b_pos[i];
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
  last_a = -1;
  last_b = -1;
  for (i = 0; i < len; i++) {
    if ((int) b_pos[i] > last_b) {
      fprintf (ofs, "%c", n2c[b->pos[b_pos[i]].nucl]);
    } else {
      fprintf (ofs, "-");
    }
    last_a = a_pos[i];
    last_b = b_pos[i];
  }
  /* Final unaligned part */
  for (i = b_pos[len - 1] + 1; i < b->len; i++) {
    fprintf (ofs, "%c", n2c[b->pos[i].nucl]);
  }
  fprintf (ofs, "\n");
}

/* Smith-Waterman algorithm */

#define M_SCORE 1
#define N_SCORE 0
#define MM_SCORE -2
#define GAP_SCORE -4

#define NROWS (n + 1)
#define NCOLS (m + 1)

typedef struct _SWCell SWCell;
struct _SWCell {
  int score : 24;
  int sx : 4;
  int sy : 4;
};

static unsigned int
smith_waterman_seq (unsigned int a_pos[], unsigned int b_pos[], const NSeq *a, const NSeq *b)
{
  unsigned int i, j;
  unsigned int n = a->len, m = b->len;
  /* Matrix of (n + 1) rows and (m + 1) columns */
  static SWCell *t = NULL;
  static unsigned int t_size = 0;
  if ((NROWS * NCOLS) > t_size) {
    unsigned int h = (NROWS + 15) & 0xfffffff0;
    unsigned int w = (NCOLS + 15) & 0xfffffff0;
    t = (SWCell *) realloc (t, h * w * sizeof (SWCell));
    t_size = h * w;
  }
  /* Fill first column with zeroes */
  /* Fill first row with zeroes */
  memset (t, 0, NROWS * NCOLS * sizeof (SWCell));
  /* Fill table starting from 1-st row and pick maximum score */
  int max_i = 0, max_j = 0;
  for (i = 1; i <= n; i++) {
    for (j = 1; j <= m; j++) {
      int score = ((a->pos[i - 1].nucl >= N) || (b->pos[j - 1].nucl >= N)) ? N_SCORE : (a->pos[i - 1].nucl == b->pos[j - 1].nucl) ? M_SCORE : MM_SCORE;
      t[i * NCOLS + j].score = 0;
      /* Diagonal */
      if ((t[(i - 1) * NCOLS + (j - 1)].score + score) > t[i * NCOLS + j].score) {
        t[i * NCOLS + j].score = t[(i - 1) * NCOLS + (j - 1)].score + score;
        t[i * NCOLS + j].sx = -1;
        t[i * NCOLS + j].sy = -1;
      }
      /* Left */
      if ((t[i * NCOLS + (j - 1)].score + GAP_SCORE) > t[i * NCOLS + j].score) {
        t[i * NCOLS + j].score = t[i * NCOLS + (j - 1)].score + GAP_SCORE;
        t[i * NCOLS + j].sx = -1;
        t[i * NCOLS + j].sy = 0;
      }
      /* Top */
      if ((t[(i - 1) * NCOLS + j].score + GAP_SCORE) > t[i * NCOLS + j].score) {
        t[i * NCOLS + j].score = t[(i - 1) * NCOLS + j].score + GAP_SCORE;
        t[i * NCOLS + j].sx = 0;
        t[i * NCOLS + j].sy = -1;
      }
      if (t[i * NCOLS + j].score > t[max_i * NCOLS + max_j].score) {
        max_i = i;
        max_j = j;
      }
    }
  }
  unsigned int len = 0;
  while ((max_i > 0) && (max_j > 0)) {
    a_pos[len] = max_i - 1;
    b_pos[len] = max_j - 1;
    len += 1;
    int sx = t[max_i * NCOLS + max_j].sx;
    int sy = t[max_i * NCOLS + max_j].sy;
    if (!sx && !sy) break;
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
  return len;
}

static unsigned int
smith_waterman (unsigned int a_pos[], unsigned int b_pos[], const unsigned int a[], unsigned int n, const unsigned int b[], unsigned int m)
{
  unsigned int i, j;
  /* Matrix of (n + 1) rows and (m + 1) columns */
  int *t = (int *) malloc (NROWS * NCOLS * 4);
  /* Fill first column with zeroes */
  for (i = 0; i <= n; i++) t[i * NCOLS + 0] = 0;
  /* Fill first row with zeroes */
  for (j = 0; j <= m; j++) t[0 * NCOLS + j] = 0;
  /* Fill table starting from 1-st row and pick maximum score */
  int max_i = 0, max_j = 0;
  for (i = 1; i <= n; i++) {
    for (j = 1; j <= m; j++) {
      int score = ((a[i - 1] == N) || (b[j - 1] == N)) ? N_SCORE : (a[i - 1] == b[j - 1]) ? M_SCORE : MM_SCORE;
      /* Diagonal */
      t[i * NCOLS + j] = t[(i - 1) * NCOLS + (j - 1)] + score;
      /* Left */
      if ((t[i * NCOLS + (j - 1)] + GAP_SCORE) > t[i * NCOLS + j]) t[i * NCOLS + j] = t[i * NCOLS + (j - 1)] + GAP_SCORE;
      /* Top */
      if ((t[(i - 1) * NCOLS + j] + GAP_SCORE) > t[i * NCOLS + j]) t[i * NCOLS + j] = t[(i - 1) * NCOLS + j] + GAP_SCORE;
      if (t[i * NCOLS + j] > t[max_i * NCOLS + max_j]) {
        max_i = i;
        max_j = j;
      }
    }
  }
  unsigned int len = 0;
  while ((max_i > 0) && (max_j > 0)) {
    a_pos[len] = max_i - 1;
    b_pos[len] = max_j - 1;
    len += 1;
    /* Diagonal */
    int best_i = max_i - 1;
    int best_j = max_j - 1;
    /* Left */
    if (t[max_i * NCOLS + (max_j - 1)] > t[best_i * NCOLS + best_j]) {
      best_i = max_i;
      best_j = max_j - 1;
    }
    /* Top */
    if (t[(max_i - 1) * NCOLS + max_j] > t[best_i * NCOLS + best_j]) {
      best_i = max_i - 1;
      best_j = max_j;
    }
    max_i = best_i;
    max_j = best_j;
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
  free (t);
  return len;
}

static void
align_seq (NMatrix *mat, NCell *a[], unsigned int alen, unsigned int seq_idx)
{
  unsigned int an[4096];
  unsigned int bn[4096];
  unsigned int a_pos[4096], b_pos[4096];
  unsigned int i, j;
  for (j = 0; j < alen; j++) {
    an[j] = a[j]->consensus;
  }
  for (i = 0; i < mat->seqs[seq_idx]->len; i++) {
    bn[i] = mat->seqs[seq_idx]->pos[i].nucl;
  }
  unsigned int len = smith_waterman (a_pos, b_pos, an, alen, bn, mat->seqs[seq_idx]->len);
  for (i = 0; i < alen; i++) {
    if (a[i]->links[seq_idx].pos >= 0) n_matrix_unlink_cell (mat, a[i], seq_idx);
  }
  for (i = 0; i < len; i++) {
    if (a_pos[i] > 10000) continue;
    if (b_pos[i] > 10000) continue;
    if (a[a_pos[i]]->links[seq_idx].pos < 0) n_matrix_link_cell (mat, a[a_pos[i]], seq_idx, b_pos[i]);
  }
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
        unsigned int nstok;
        char chr[32];
        nstok = split_line_chr (tokenz[0], lengths[0], stok, slen, 5, ':');
        if (slen[0] > 31) slen[0] = 31;
        memcpy (chr, stok[0], slen[0]);
        chr[slen[0]] = 0;
        snvs[*n_snvs].chr = chr_from_text (chr);
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
    } else if ((snvs[mid].chr < chr) || (snvs[mid].pos < pos)) {
      min = mid;
    } else if ((snvs[mid].chr > chr) || (snvs[mid].pos > pos)) {
      max = mid;
    } else {
      break;
    }
    mid = (min + max) / 2;
  }
  return mid;
}

static void
load_db_or_die (KMerDB *db, const char *db_name, const char *id)
{
  const unsigned char *cdata;
  unsigned long long csize;
  /* Read database */
  if (debug) fprintf (stderr, "Loading %s database %s... ", id, db_name);
  cdata = gt4_mmap (db_name, &csize);
  if (!cdata) {
    fprintf (stderr, "cannot mmap (no such file?)\n");
    exit (1);
  }
  /* scout_mmap (cdata, csize); */
  if (!read_database_from_binary (db, cdata, csize)) {
    fprintf (stderr, "cannot read (wrong file format?)\n");
    exit (1);
  }
  if (!db->index.read_blocks) {
    fprintf (stderr, "no index\n");
    exit (1);
  }
  if (debug) fprintf (stderr, "done\n");
}

static SeqFile *
map_sequences (KMerDB *db)
{
  unsigned int i;
  SeqFile *files = (SeqFile *) malloc (db->index.n_files * sizeof (SeqFile));
  for (i = 0; i < db->index.n_files; i++) {
    files[i].name = db->index.files[i];
    if (!files[i].cdata) {
      files[i].cdata = gt4_mmap (files[i].name, &files[i].csize);
      if (!files[i].cdata) {
        fprintf (stderr, "Cannot memory map %s\n", files[i].name);
        free (files);
        return NULL;
      }
    }
  }
  return files;
}

static unsigned int
get_unique_reads (ReadInfo reads[], unsigned int max_reads, KMerDB *db, SeqFile files[], const char *kmers[], unsigned int nkmers)
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
    if (n_reads > MAX_READS_PER_KMER) {
      fprintf (stderr, "Kmer %u has too many reads: %u\n", i, n_reads);
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
        if (nreads >= max_reads) break;
      } else {
        if (debug > 2) fprintf (stderr, "  Already registered as %u\n", k);
      }
    }
    if (debug) fprintf (stderr, "Kmer %u %s reads %u new %u\n", i, kmers[i], n_reads, n_new_reads);
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
    if (len > 2047) len = 2047;
    memcpy (seq, p, len);
    seq[len] = 0;
    if (reads[i].dir) gt4_string_revcomp_inplace (seq, len);
    seqs[i] = gasm_read_new (name, seq, WORDLEN);
    if (debug) fprintf (stderr, "Read %2u: %s\n", i, seqs[i]->seq);
  }
  return 1;
}

unsigned int
remove_bad_reads (GASMRead *seqs[], unsigned int nseqs, KMerDB *gdb)
{
  unsigned int i;
  /* Sanitize */
  unsigned int idx = 0;
  i = 0;
  while (i < nseqs) {
    unsigned int j, invalid = 0;
    NSeq *seq = seqs[i]->nseq;
    unsigned long long ref_pos = 0;
    for (j = 0; j < seq->len; j++) {
      if (seq->pos[j].has_kmer) {
        unsigned int num_seqs = 0, file_idx = 0, dir = 0;
#if 0
        if (seq->pos[j].non_unique_kmer) {
          fprintf (stderr, "Invalid read %u: pos %u has non-unique k-mer %s\n", idx, j, word_to_string (seq->pos[j].kmer, WORDLEN));
          fprintf (stderr, "Seq: %s\n", seqs[i]->seq);
          invalid = 1;
          break;
        }
#endif
        unsigned long long kmer_pos = get_kmer_location (gdb, seq->pos[j].kmer, &num_seqs, &file_idx, &dir);
        if (num_seqs == 1) {
#if 0
          if (ref_pos) {
            long long delta = (long long) kmer_pos - (long long) ref_pos;
            if (delta < 0) delta = -delta;
            if (delta > 16) {
              fprintf (stderr, "Invalid read %u: pos %u reference gap > 16 bp at(%llu / %llu)\n", idx, j, ref_pos, kmer_pos);
              invalid = 1;
              break;
            }
          }
#endif
          ref_pos = kmer_pos;
        } else if (num_seqs > 10) {
          fprintf (stderr, "Invalid read %u: pos %u k-mer has %u reference locations\n", idx, j, num_seqs);
          invalid = 1;
          break;
        }
      }
      if (ref_pos) ref_pos += 1;
    }
    if (invalid) {
      gasm_read_delete (seqs[i]);
      nseqs -= 1;
      seqs[i] = seqs[nseqs];
    } else {
      i += 1;
    }
    idx += 1;
  }
  return nseqs;
}

static unsigned long long
get_kmer_location (KMerDB *gdb, unsigned long long word, unsigned int *num_seqs, unsigned int *file_idx, unsigned int *dir)
{
  unsigned long long kmer_pos = 0;
  unsigned int code = trie_lookup (&gdb->trie, word);
  if (code) {
    *dir = ((code & 0x8000000) != 0);
    code &= 0x7fffffff;
    unsigned int node = (code >> gdb->kmer_bits) - 1;
    unsigned int kmer = code & ((1 << gdb->kmer_bits) - 1);
    unsigned int kmer_idx = gdb->nodes[node].kmers + kmer;
    unsigned long long first_read;
    first_read = gt4_index_get_kmer_info (&gdb->index, kmer_idx, num_seqs);
    if (*num_seqs == 1) {
      unsigned long long name_pos;
      kmer_pos = gt4_index_get_read_info (&gdb->index, first_read, file_idx, &name_pos, dir);
    }
  }
  return kmer_pos;
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

static void
print_kmer_statistics (NMatrix *mat)
{
  unsigned int i;
  for (i = 0; i < mat->n_unique_kmers; i++) {
    char c[32] = { 0 };
    unsigned int j, k, l;
    unsigned int t[1024] = { 0 };
    for (j = 0; j < mat->n_kmers; j++) {
      if (mat->kmers[j].value != mat->unique_kmers[i].value) continue;
      for (k = 0; k < mat->n_kmers; k++) {
        if (mat->kmers[k].value == mat->kmers[j].value) continue;
        if (mat->kmers[k].seq != mat->kmers[j].seq) continue;
        for (l = 0; l < mat->n_unique_kmers; l++) {
          if (mat->unique_kmers[l].value == mat->kmers[k].value) {
            /* Kmers i & l are in the same read */
            t[l] += 1;
            break;
          }
        }
      }
    }
    word2string (c, mat->unique_kmers[i].value, WORDLEN);
    fprintf (stdout, "%u\t%s\t%u", i, c, mat->unique_kmers[i].count);
    for (j = 0; j < mat->n_unique_kmers; j++) {
      fprintf (stdout, "\t%u", t[j]);
    }
    fprintf (stdout, "\n");
  }
}

