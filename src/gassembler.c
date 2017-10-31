#define __GASSEMBLER_C__

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <database.h>
#include <matrix.h>
#include <queue.h>
#include <sequence.h>
#include <utils.h>

int debug = 1;

#define WORDLEN 25
#define MAX_KMERS 1024
#define MAX_READS 1024

typedef struct _SeqFile SeqFile;

struct _SeqFile {
  const char *name;
  const unsigned char *cdata;
  unsigned long long csize;
};

static void load_db_or_die (KMerDB *db, const char *db_name, const char *id);
static SeqFile *map_sequences (KMerDB *db);

static void assemble (KMerDB *db, SeqFile *files, KMerDB *gdb, SeqFile *g_files, const char *kmers[], unsigned int nkmers, const char *ref_chr, unsigned int ref_start, unsigned int ref_end, const char *ref);

/* Get list of unique reads containing at lest one kmer from list */
static unsigned int get_unique_reads (KMerDB *db, Read reads[], SeqFile files[], const char *kmers[], unsigned int nkmers);
/* Get proper read sequences (forward direction) */
static unsigned int get_read_sequences (char *seqs[], const Read reads[], unsigned int nreads, SeqFile files[]);
/* Remove reads that have multiple instances of same k-mer, too long gaps or too comon k-mers */
unsigned int remove_bad_reads (const char *seqs[], unsigned int nseqs, KMerDB *gdb);
static unsigned long long get_kmer_location (KMerDB *gdb, unsigned long long word, unsigned int *num_seqs, unsigned int *file_idx, unsigned int *dir);
static void print_db_reads (GT4Index *index, SeqFile *files, unsigned long long kmer_idx, unsigned int kmer_dir, FILE *ofs);
static void print_kmer_statistics (NMatrix *mat);

static unsigned int create_alignment (unsigned int a_pos[], unsigned int b_pos[], const unsigned int a[], unsigned int na, const unsigned int b[], unsigned int nb);
static unsigned int build_alignment (NCell *alignment[], unsigned int alignment_size, KMerDB *db, KMerDB *gdb, NMatrix *mat, unsigned int min_cov);
static void align_seq (NMatrix *mat, NCell *a[], unsigned int alen, unsigned int seq_idx);

const char *db_name = NULL;
const char *gdb_name = NULL;
static unsigned int print_chains = 0;
static unsigned int print_reads = 0, print_kmers = 0, analyze_kmers = 0;

int
main (int argc, const char *argv[])
{
  unsigned int i;
  const char *input_name = NULL;
  const char *kmers[MAX_KMERS];
  unsigned int nkmers = 0;
  const char *ref_chr = NULL;
  unsigned int ref_start = 0, ref_end = 0;
  const char *ref = NULL;
  unsigned int max_regions = 1000000000;

  KMerDB db, gdb;
  SeqFile *files, *g_files;
    
  for (i = 1; i < argc; i++) {
    if (!strcmp (argv[i], "-dbb") || !strcmp (argv[i], "-db")) {
      i += 1;
      if (i >= argc) exit (1);
      db_name = argv[i];
    } else if (!strcmp (argv[i], "-gdb")) {
      i += 1;
      if (i >= argc) exit (1);
      gdb_name = argv[i];
    } else if (!strcmp (argv[i], "-reference")) {
      if ((i + 4) >= argc) exit (1);
      ref_chr = argv[i + 1];
      ref_start = atoi (argv[i + 2]);
      ref_end = atoi (argv[i + 3]);
      ref = (const char *) argv[i + 4];
      i += 4;
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
    fprintf (stderr, "No DB/GDB specified\n");
    exit (1);
  }

  /* Read databases */
  load_db_or_die (&db, db_name, "reads");
  load_db_or_die (&gdb, gdb_name, "genome");

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
        ref_chr = strndup ((const char *) tokenz[0], lengths[0]);
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

static void
assemble (KMerDB *db, SeqFile *files, KMerDB *gdb, SeqFile *g_files, const char *kmers[], unsigned int nkmers, const char *ref_chr, unsigned int ref_start, unsigned int ref_end, const char *ref)
{
  NMatrix *mat;
  char **seqs;
  Read reads[8192];
  unsigned int nreads, i;

  /* Print virtual command line to simplify debugging */
  if (debug) {
    fprintf (stderr, "Arguments: -db %s -gdb %s -reference %s %u %u %s", db_name, gdb_name, ref_chr, ref_start, ref_end, ref);
    for (i = 0; i < nkmers; i++) {
      fprintf (stderr, " %s", kmers[i]);
    }
    fprintf (stderr, "\n");
  }

  /* Get all unique reads */
  nreads = get_unique_reads (db, reads, files, kmers, nkmers);
  if (debug) fprintf (stderr, "Got %u unique reads\n", nreads);

  /* Create actual sequences in correct direction */
  seqs = (char **) malloc (nreads * sizeof (char *));
  get_read_sequences (seqs, reads, nreads, files);

  /* Now we have bunch of reads */
  if (print_reads) {
    for (i = 0; i < nreads; i++) {
      fprintf (stdout, ">Read_%u\n", i);
      fprintf (stdout, "%s\n", seqs[i]);
    }
  }
  if (debug) fprintf (stderr, "Number of unique reads: %u\n", nreads);

  /* Sanitize */
  nreads = remove_bad_reads ((const char **) seqs, nreads, gdb);

  if (debug) fprintf (stderr, "Number of usable reads: %u\n", nreads);
  /* Now we have bunch of reads */
  if (print_reads) {
    for (i = 0; i < nreads; i++) {
      fprintf (stdout, ">Read_%u\n", i);
      fprintf (stdout, "%s\n", seqs[i]);
    }
  }
  if (nreads < 10) {
    fprintf (stderr, "Final number of reads too low (%u)\n", nreads);
    return;
  } else if (nreads > MAX_READS) {
    fprintf (stderr, "Final number of reads too big (%u)\n", nreads);
    return;
  }

  /* Create matrix */
  mat = n_matrix_new (nreads, (const char **) seqs, WORDLEN);
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
  unsigned int alen = build_alignment (alignment, 2048, db, gdb, mat, 15);
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
  unsigned int n_ref_align = create_alignment (a_pos, b_pos, a, alen, b, b_seq->len);
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

static unsigned int
build_alignment (NCell *alignment[], unsigned int alignment_size, KMerDB *db, KMerDB *gdb, NMatrix *mat, unsigned int min_cov)
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

#define MM_PENALTY -1
#define GAP_PENALTY -4

static unsigned int
create_alignment (unsigned int a_pos[], unsigned int b_pos[], const unsigned int a[], unsigned int na, const unsigned int b[], unsigned int nb)
{
  unsigned int i, j;
  int *t = (int *) malloc (na * nb * 4);
  for (j = 0; j < na; j++) {
    int score = 0;
    if ((a[j] != N) && (a[j] == b[0])) score = 1;
    t[0 * na + j] = score;
  }
  for (i = 0; i < nb; i++) {
    int score = 0;
    if ((b[i] != N) && (b[i] == a[0])) score = 1;
    t[i * na + 0] = score;
  }
  for (i = 1; i < nb; i++) {
    for (j = 1; j < na; j++) {
      int score = t[(i - 1) * na + (j - 1)];
      if ((t[i * na + (j - 1)] + GAP_PENALTY) > score) score = t[i * na + (j - 1)] + GAP_PENALTY;
      if ((t[(i - 1) * na + j] + GAP_PENALTY) > score) score = t[(i - 1) * na + j] + GAP_PENALTY;
      if ((a[j] != N) && (b[i] != N)) {
        if (a[j] == b[i]) {
          score += 1;
        } else {
          score += MM_PENALTY;
        }
      }
      /* if (score < 0) score = 0; */
      t[i * na + j] = score;
    }
  }

  int max_i = nb - 1, max_j = 0;
  for (j = 0; j < na; j++) {
    if (t[(nb - 1) * na + j] > t[(nb - 1) * na + max_j]) {
      max_i = nb - 1;
      max_j = j;
    }
  }
  for (i = 0; i < nb; i++) {
    if (t[i * na + (na - 1)] > t[max_i * na + max_j]) {
      max_i = i;
      max_j = na - 1;
    }
  }
  /* fprintf (stderr, "%d %d\n", max_i, max_j); */
  unsigned int len = 0;
  a_pos[0] = max_j;
  b_pos[0] = max_i;
  len += 1;
  while ((max_i >= 0) && (max_j >= 0)) {
    int score = -1;
    int best_i = i - 1;
    int best_j = j - 1;
    if ((max_i > 0) && (max_j > 0)) {
      best_i = max_i - 1;
      best_j = max_j - 1;
      score = t[best_i * na + best_j];
    }
    if ((max_i > 0) && (t[(max_i - 1) * na + max_j] > score)) {
      best_i = max_i - 1;
      best_j = max_j;
      score = t[best_i * na + best_j];
    }
    if ((max_j > 0) && (t[max_i * na + (max_j - 1)] > score)) {
      best_i = max_i;
      best_j = max_j - 1;
      score = t[best_i * na + best_j];
    }
    if (score <= 0) break;
    max_i = best_i;
    max_j = best_j;
    /* fprintf (stderr, "%d %d\n", max_i, max_j); */
    a_pos[len] = max_j;
    b_pos[len] = max_i;
    len += 1;
  }
  for (i = 0; i < len; i++) {
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
  unsigned int len = create_alignment (a_pos, b_pos, an, alen, bn, mat->seqs[seq_idx]->len);
  for (i = 0; i < alen; i++) {
    if (a[i]->links[seq_idx].pos >= 0) n_matrix_unlink_cell (mat, a[i], seq_idx);
  }
  for (i = 0; i < len; i++) {
    if (a_pos[i] > 10000) continue;
    if (b_pos[i] > 10000) continue;
    if (a[a_pos[i]]->links[seq_idx].pos < 0) n_matrix_link_cell (mat, a[a_pos[i]], seq_idx, b_pos[i]);
  }
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
  scout_mmap (cdata, csize);
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

#define MAX_READS_PER_KMER 100

static unsigned int
get_unique_reads (KMerDB *db, Read reads[], SeqFile files[], const char *kmers[], unsigned int nkmers)
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
        if (nreads >= MAX_READS) break;
      } else {
        if (debug > 2) fprintf (stderr, "  Already registered as %u\n", k);
      }
    }
    if (debug) fprintf (stderr, "Kmer %u %s reads %u new %u\n", i, kmers[i], n_reads, n_new_reads);
  }
  return nreads;
}

static unsigned int
get_read_sequences (char *seqs[], const Read reads[], unsigned int nreads, SeqFile files[])
{
  unsigned int i;
  /* Create actual sequences */
  /* Change directionality of reads if needed */
  for (i = 0; i < nreads; i++) {
    if (!files[reads[i].file_idx].cdata) {
      files[reads[i].file_idx].cdata = gt4_mmap (files[reads[i].file_idx].name, &files[reads[i].file_idx].csize);
      if (!files[reads[i].file_idx].cdata) {
        fprintf (stderr, "Cannot mmap %s\n", files[reads[i].file_idx].name);
        return 0;
      }
    }
    const unsigned char *p = files[reads[i].file_idx].cdata + reads[i].name_pos;
    unsigned int len;
    while (*p != '\n') p += 1;
    p += 1;
    len = 0;
    while (p[len] >= 'A') len += 1;
    seqs[i] = malloc (len + 1);
    memcpy (seqs[i], p, len);
    seqs[i][len] = 0;
    if (reads[i].dir) gt4_string_revcomp_inplace (seqs[i], len);
    fprintf (stderr, "%u: %s\n", i, seqs[i]);
  }
  return 1;
}

unsigned int
remove_bad_reads (const char *seqs[], unsigned int nseqs, KMerDB *gdb)
{
  unsigned int i;
  /* Sanitize */
  unsigned int idx = 0;
  i = 0;
  while (i < nseqs) {
    unsigned int j, invalid = 0;
    NSeq *seq = n_seq_new (seqs[i], WORDLEN);
    unsigned long long ref_pos = 0;
    for (j = 0; j < seq->len; j++) {
      if (seq->pos[j].has_kmer) {
        unsigned int num_seqs = 0, file_idx = 0, dir = 0;
        if (seq->pos[j].non_unique_kmer) {
          fprintf (stderr, "Invalid read %u: pos %u has non-unique k-mer %s\n", idx, j, word_to_string (seq->pos[j].kmer, WORDLEN));
          fprintf (stderr, "Seq: %s\n", seqs[i]);
          invalid = 1;
          break;
        }
        unsigned long long kmer_pos = get_kmer_location (gdb, seq->pos[j].kmer, &num_seqs, &file_idx, &dir);
        if (num_seqs == 1) {
          if (ref_pos) {
            long long delta = (long long) kmer_pos - (long long) ref_pos;
            if (delta < 0) delta = -delta;
            if (delta > 16) {
              fprintf (stderr, "Invalid read %u: pos %u reference gap > 16 bp at(%llu / %llu)\n", idx, j, ref_pos, kmer_pos);
              invalid = 1;
              break;
            }
          }
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

