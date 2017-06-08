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

typedef struct _SeqFile SeqFile;

struct _SeqFile {
  const char *name;
  const unsigned char *cdata;
  unsigned long long csize;
};

static void load_db_or_die (KMerDB *db, const char *db_name, const char *id);
static unsigned int get_reads (KMerDB *db, KMerDB *gdb, Read reads[], SeqFile files[], const char *kmers[], unsigned int nkmers);
static void print_db_reads (GT4Index *index, SeqFile *files, unsigned long long kmer_idx, unsigned int kmer_dir, FILE *ofs);
static void print_kmer_statistics (NMatrix *mat);

static unsigned int create_alignment (unsigned int a_pos[], unsigned int b_pos[], const unsigned int a[], unsigned int na, const unsigned int b[], unsigned int nb);
static unsigned int build_alignment (NCell *alignment[], unsigned int alignment_size, KMerDB *db, KMerDB *gdb, NMatrix *mat, unsigned int min_cov);
static void align_seq (NMatrix *mat, NCell *a[], unsigned int alen, unsigned int seq_idx);

static unsigned int print_chains = 0;

int
main (int argc, const char *argv[])
{
  unsigned int i;
  const char *db_name = NULL;
  const char *gdb_name = NULL;
  const char *kmers[1024];
  unsigned int nkmers = 0;
  const unsigned char *ref_chr = NULL;
  unsigned int ref_start = 0, ref_end = 0;
  const char *ref = NULL;
  unsigned int print_reads = 0, print_kmers = 0, analyze_kmers = 0;
  unsigned int die = 0;

  KMerDB db, gdb;
  SeqFile *files, *g_files;
  unsigned int nreads;
  Read reads[1024];
  char *seqs[1024];
  NMatrix *mat;
    
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
      ref_chr = (const unsigned char *) argv[i + 1];
      ref_start = atoi (argv[i + 2]);
      ref_end = atoi (argv[i + 3]);
      ref = (const unsigned char *) argv[i + 4];
      i += 4;
    } else if (!strcmp (argv[i], "--print_reads")) {
      print_reads = 1;
    } else if (!strcmp (argv[i], "--print_kmers")) {
      print_kmers = 1;
    } else if (!strcmp (argv[i], "--print_chains")) {
      print_chains = 1;
    } else if (!strcmp (argv[i], "--analyze_kmers")) {
      analyze_kmers = 1;
    } else  {
      kmers[nkmers++] = argv[i];
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
  files = (SeqFile *) malloc (db.index.n_files * sizeof (SeqFile));
  for (i = 0; i < db.index.n_files; i++) {
    files[i].name = db.index.files[i];
    if (!files[i].cdata) {
      files[i].cdata = gt4_mmap (files[i].name, &files[i].csize);
      if (!files[i].cdata) {
        fprintf (stderr, "Cannot load read file %s\n", files[i].name);
        die = 1;
      }
    }
  }
  g_files = (SeqFile *) malloc (gdb.index.n_files * sizeof (SeqFile));
  for (i = 0; i < gdb.index.n_files; i++) {
    g_files[i].name = gdb.index.files[i];
    if (!g_files[i].cdata) {
      g_files[i].cdata = gt4_mmap (g_files[i].name, &g_files[i].csize);
      if (!g_files[i].cdata) {
        fprintf (stderr, "Cannot load genome file %s\n", g_files[i].name);
        die = 1;
      }
    }
  }
  if (die) {
    fprintf (stderr, "Terminating\n");
    exit (1);
  }

  /* Get all unique reads */
  nreads = get_reads (&db, &gdb, reads, files, kmers, nkmers);
  if (debug) fprintf (stderr, "Total %u unique reads\n", nreads);

  /* Create actual sequences */
  /* Change directionality of reads if needed */
  for (i = 0; i < nreads; i++) {
    if (!files[reads[i].file_idx].cdata) {
      files[reads[i].file_idx].cdata = gt4_mmap (files[reads[i].file_idx].name, &files[reads[i].file_idx].csize);
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
    fprintf (stderr, "%u %s\n", i, seqs[i]);
    if (reads[i].dir) gt4_string_revcomp_inplace (seqs[i], len);
  }

  /* Now we have bunch of reads */
  if (print_reads) {
    for (i = 0; i < nreads; i++) {
      fprintf (stdout, ">Read_%u\n", i);
      fprintf (stdout, "%s\n", seqs[i]);
    }
  }

  /* Create matrix */
  mat = n_matrix_new (nreads, (const char **) seqs, WORDLEN);
  if (debug > 2) {
    for (i = 0; i < mat->n_seqs; i++) {
      unsigned int j;
      fprintf (stderr, "Seq %u\n%s\n", i, mat->seqs[i]->str);
      for (j = 0; j < mat->seqs[i]->len; j++) {
        if (mat->seqs[i]->pos[j].has_kmer) {
          char c[32];
          unsigned int k;
          word2string (c, mat->seqs[i]->pos[j].kmer, WORDLEN);
          c[WORDLEN] = 0;
          for (k = 0; k < j; k++) fprintf (stderr, " ");
          fprintf (stderr, "%s\n", c);
        }
      }
    }
  }
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
  unsigned int alen = build_alignment (alignment, 2048, &db, &gdb, mat, 15);
  unsigned int a[2048];
  unsigned int b[2048];
  unsigned int a_pos[2048];
  unsigned int b_pos[2048];
  unsigned int r[2048] = { 0 };
  for (i = 0; i < alen; i++) a[i] = alignment[i]->consensus;
  NSeq *b_seq = n_seq_new (ref, WORDLEN);
  for (i = 0; i <= b_seq->len; i++) b[i] = b_seq->pos[i].nucl;
  unsigned int n_ref_align = create_alignment (a_pos, b_pos, a, alen, b, b_seq->len);
  for (i = 0; i < n_ref_align; i++) {
    r[a_pos[i]] = b_seq->pos[b_pos[i]].nucl;
    alignment[a_pos[i]]->ref_pos = ref_start + b_pos[i];
  }
  
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
    } else if (has_rgap) {
      pos[ac_len] = next_rp;
      ac[ac_len] = GAP;
      rc[ac_len] = b_seq->pos[next_rp - ref_start].nucl;
      next_rp += 1;
      ac_len += 1;
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
      i += 1;
    }
  }

  fprintf (stdout, "CHR\tREF\tPOS\tCONS\tTOTAL\tA\tC\tG\tT\tN\tGAP\tALLELES\n");
  for (i = 0; i < ac_len; i++) {
    unsigned int j, sum;
    sum = 0;
    for (j = 0; j < 6; j++) sum += nc[6 * i + j];
    /* if (sum < 10) continue;*/
    if ((sum > 0) && (nc[6 * i + GAP] >= (2 * sum / 3))) continue;
    if (pos[i] > 0) {
      fprintf (stdout, "%s\t%u\t", ref_chr, pos[i]);
    } else {
      fprintf (stdout, "\t\t");
    }
    fprintf (stdout, "%c\t%c\t%u", n2c[rc[i]], n2c[ac[i]], sum);
    unsigned int s[6] = {0};
    for (j = 0; j < 6; j++) {
      if ((sum >= 10) && (nc[6 * i + j] > (sum / 3))) s[j] = 1;
      fprintf (stdout, "\t%d", nc[6 * i + j]);
    }
    fprintf (stdout, "\t");
    for (j = 0; j < 6; j++) {
      if (s[j]) fprintf (stdout, "%c", n2c[j]);
    }
    fprintf (stdout, "\n");
  }

  return 0;
}

static unsigned int
build_alignment (NCell *alignment[], unsigned int alignment_size, KMerDB *db, KMerDB *gdb, NMatrix *mat, unsigned int min_cov)
{
  unsigned int i;
  NCell *cell;
  /* Fill matrix */
  fprintf (stderr, "Building matrix\n");
  /* Iterate over all unique kmers */
  for (i = 0; i < mat->n_unique_kmers; i++) {
    NCell *cells[WORDLEN] = { 0 };
    unsigned int j;

    unsigned int file_idx = 0;
    unsigned long long kmer_pos = 0;
    unsigned int dir = 0;
    unsigned int code = trie_lookup (&gdb->trie, mat->unique_kmers[i].value);
    if (code) {
      dir = ((code & 0x8000000) != 0);
      code &= 0x7fffffff;
      unsigned int node = (code >> gdb->kmer_bits) - 1;
      unsigned int kmer = code & ((1 << gdb->kmer_bits) - 1);
      unsigned int kmer_idx = gdb->nodes[node].kmers + kmer;
      unsigned long long first_read;
      unsigned int num_reads;
      first_read = gt4_index_get_kmer_info (&gdb->index, kmer_idx, &num_reads);
      if (num_reads == 1) {
        unsigned long long name_pos;
        kmer_pos = gt4_index_get_read_info (&gdb->index, first_read, &file_idx, &name_pos, &dir);
      }
    }
    
    if (print_chains) {
      fprintf (stdout, "Kmer %u, reads", i);
    }
    for (j = 0; j < mat->n_kmers; j++) {
      unsigned int seq_idx, pos, k;
      if (mat->kmers[j].value != mat->unique_kmers[i].value) continue;
      seq_idx = mat->kmers[j].seq;
      pos = mat->kmers[j].pos;
      for (k = 0; k < WORDLEN; k++) {
        if (!cells[k]) cells[k] = mat->seqs[seq_idx]->pos[pos + k].cells;
      }
    }
    /* Add all instances of this kmer */
    for (j = 0; j < mat->n_kmers; j++) {
      unsigned int seq_idx, pos, k;
      if (mat->kmers[j].value != mat->unique_kmers[i].value) continue;
      seq_idx = mat->kmers[j].seq;
      pos = mat->kmers[j].pos;
      if (print_chains) fprintf (stdout, " %u", seq_idx);
      for (k = 0; k < WORDLEN; k++) {
        if (!cells[k]) {
          if (debug > 2) fprintf (stderr, "Create cell: seq %u start %u k %u pos %u\n", seq_idx, pos, k, pos + k);
          assert (!mat->seqs[seq_idx]->pos[pos + k].cells);
          cells[k] = n_matrix_new_cell (mat);
        }
        if (kmer_pos > 0) {
          cells[k]->file_idx = file_idx;
          cells[k]->ref_pos = (dir) ? kmer_pos + WORDLEN - 1 - k : kmer_pos + k;
        }
        if (!mat->seqs[seq_idx]->pos[pos + k].cells) {
          /* Link cell from current kmer to seq/pos */
          n_matrix_link_cell (mat, cells[k], seq_idx, pos + k);
        }
      }
    }
    if (print_chains) fprintf (stdout, "\n");
  }

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
        cell->ref_pos = 666;
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
    fprintf (stderr, "Max prev %u max next %u\n", max_prev, max_next);
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
  fprintf (stderr, "Finding fist cell\n");
  NCell *first_cell = max_cell;
  while (first_cell->prev) first_cell = first_cell->prev;
  unsigned int alen = 0;
  cell = first_cell;
  while (cell) {
    alen += 1;
    cell = cell->next;
  }
  if (debug) fprintf (stderr, "Initial alignment length %u\n", alen);
  
  /* Compact alignment */
  cell = first_cell;
  while (cell->next) {
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
  unsigned int apos = 0;
  for (cell = first_cell; cell; cell = cell->next) {
    alignment[apos++] = cell;
  }

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
            if (n_matrix_can_merge_cells (mat, alignment[k], cc)) {
              fprintf (stderr, "Merging cells: alignment pos %u seq %u pos %u\n", i, j, s_p);
              n_matrix_merge_cells (mat, alignment[k], cc);
              break;
            } else {
              fprintf (stderr, "Cannot merge: alignment pos %u seq %u pos %u\n", i, j, s_p);
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
            if (n_matrix_can_merge_cells (mat, alignment[k], cc)) {
              fprintf (stderr, "Merging cells: alignment pos %u seq %u pos %u\n", i, j, s_p);
              n_matrix_merge_cells (mat, alignment[k], cc);
              break;
            } else {
              fprintf (stderr, "Cannot merge: alignment pos %u seq %u pos %u\n", i, j, s_p);
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

#if 0
  char *n = "ACGTN";
  fprintf (stderr, " ");
  for (j = 0; j < na; j++) {
    fprintf (stderr, "   %c", n[a[j]]);
  }
  fprintf (stderr, "\n");
  for (i = 0; i < nb; i++) {
    fprintf (stderr, "%c", n[b[i]]);
    for (j = 0; j < na; j++) {
      fprintf (stderr, " %3d", t[i * na + j]);
    }
    fprintf (stderr, "\n");
  }
  fprintf (stderr, "\n");
#endif

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
    if ((max_i > 0) && (max_j >= 0)) {
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
#if 0
  for (i = 0; i < len; i++) {
    fprintf (stderr, "%c", n[a[a_pos[i]]]);
  }
  fprintf (stderr, "\n");
  for (i = 0; i < len; i++) {
    fprintf (stderr, "%c", n[b[b_pos[i]]]);
  }
  fprintf (stderr, "\n");
#endif
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
    if (a[a_pos[i]]->links[seq_idx].pos < 0) n_matrix_link_cell (mat, a[a_pos[i]], seq_idx, b_pos[i]);
  }
}

static void
load_db_or_die (KMerDB *db, const char *db_name, const char *id)
{
  const unsigned char *cdata;
  unsigned long long csize;
  /* Read database */
  if (debug) fprintf (stderr, "Loading [%s] database %s\n", id, db_name);
  cdata = gt4_mmap (db_name, &csize);
  if (!cdata) {
    fprintf (stderr, "Cannot mmap %s\n", db_name);
    exit (1);
  }
  scout_mmap (cdata, csize);
  if (!read_database_from_binary (db, cdata, csize)) {
    fprintf (stderr, "Cannot read [%s] database %s\n", id, db_name);
    exit (1);
  }
  if (debug) fprintf (stderr, "Finished loading [%s] database (index = %u)\n", id, db->index.read_blocks != NULL);
}

static unsigned int
get_reads (KMerDB *db, KMerDB *gdb, Read reads[], SeqFile files[], const char *kmers[], unsigned int nkmers)
{
  unsigned int nreads = 0, i;

  for (i = 0; i < nkmers; i++) {
    unsigned long long first_read;
    unsigned int num_reads, j, kmer_dir;
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
    if (debug > 0) fprintf (stderr, "Kmer %s word %llu code %u\n", kmers[i], word, code);
    code &= 0x7fffffff;
    node_idx = (code >> db->kmer_bits) - 1;
    node_kmer = code & ((1 << db->kmer_bits) - 1);
    kmer_idx = db->nodes[node_idx].kmers + node_kmer;
    if (debug > 0) fprintf (stderr, "Node %u kmer %u idx %u dir %u\n", node_idx, node_kmer, kmer_idx, kmer_dir);
  
    if (debug > 1) print_db_reads (&db->index, files, kmer_idx, kmer_dir, stderr);

    first_read = gt4_index_get_kmer_info (&db->index, kmer_idx, &num_reads);
    if (debug > 0) fprintf (stderr, "Num reads %u\n", num_reads);
    for (j = 0; j < num_reads; j++) {
      unsigned long long kmer_pos, name_pos;
      unsigned int file_idx, dir;
      unsigned int k;
      kmer_pos = gt4_index_get_read_info (&db->index, first_read + j, &file_idx, &name_pos, &dir);
      for (k = 0; k < nreads; k++) {
        if ((reads[k].file_idx == file_idx) && (reads[k].name_pos == name_pos)) break;
      }
      if (k >= nreads) {
        if (debug) fprintf (stderr, "Adding read %u dir %u\n", nreads, dir);
        reads[nreads].name_pos = name_pos;
        reads[nreads].kmer_pos = kmer_pos;
        reads[nreads].file_idx = file_idx;
        /* fixme: What to do if two kmers have conflicting read directions? */
        reads[nreads].dir = (dir != kmer_dir);
        nreads += 1;
      } else {
        if (debug > 2) fprintf (stderr, "  Already registered as %u\n", k);
      }
    }
    
    if (debug > 0) {
      first_read = gt4_index_get_kmer_info (&gdb->index, kmer_idx, &num_reads);
      for (j = 0; j < num_reads; j++) {
        unsigned long long kmer_pos, name_pos;
        unsigned int file_idx, dir;
        kmer_pos = gt4_index_get_read_info (&gdb->index, first_read + j, &file_idx, &name_pos, &dir);
        fprintf (stderr, "Genome: %s:%llu\n", gdb->index.files[file_idx], kmer_pos);
      }
    }
  }
  return nreads;
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

