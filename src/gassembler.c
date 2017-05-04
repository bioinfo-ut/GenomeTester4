#define __GASSEMBLER_C__

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <database.h>
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

static void print_db_reads (GT4Index *index, SeqFile *files, unsigned long long kmer_idx, unsigned int kmer_dir, FILE *ofs);

/* Nucleotide values */
#define A 0
#define C 1
#define G 2
#define T 3
#define N 4
#define GAP 5
#define NONE 6

static const char *n2c = "ACGTN- ";

/* Position codes */
#define BEFORE -1
#define AFTER -2
#define UNKNOWN -3

typedef struct _NKMer NKMer;
typedef struct _NPos NPos;
typedef struct _NSeq NSeq;
typedef struct _NLink NLink;
typedef struct _NCell NCell;
typedef struct _NMatrix NMatrix;

struct _NKMer {
  unsigned long value;
  unsigned int seq;
  unsigned int pos;
  unsigned int count;
};

struct _NPos {
  unsigned int nucl;
  unsigned int has_kmer;
  unsigned long long kmer;
  /* Sorted list of cells aligned with this position */
  NCell *cells;
};

struct _NSeq {
  unsigned int len;
  const char *str;
  NPos pos[1];
};

struct _NLink {
  int pos;
  /* Chain of the same seq position */
  NCell *prev;
  NCell *next;
};

struct _NCell {
  /* Neighbours in alignment */
  NCell *prev;
  NCell *next;
  NCell *next_allocated;
  /* Help for path solution */
  unsigned int count;
  NLink links[1];
};

struct _NMatrix {
  unsigned int n_seqs;
  /* KMers, sorted by value, seq, pos */
  unsigned int n_kmers;
  NKMer *kmers;
  /* Unique kmers, sorted by count */
  unsigned int n_unique_kmers;
  NKMer *unique_kmers;
  /* Cells, in no particular order */
  NCell *cells;
  NCell *free_cells;
  NSeq *seqs[1];
};

static NSeq *n_seq_new (const char *str);
static unsigned int n_seq_get_kmer (NSeq *seq, unsigned int pos, unsigned long long *kmer);
static NMatrix *n_matrix_new (unsigned int n_seqs, const char *seqs[]);
static NCell *n_matrix_new_cell (NMatrix *mat);
static void n_matrix_free_cell (NMatrix *mat, NCell *cell);
static int n_matrix_compare_cells (NMatrix *mat, NCell *lhs, NCell *rhs);
static void n_matrix_link_cell (NMatrix *mat, NCell *cell, unsigned int seq, unsigned int pos);
static void n_matrix_unlink_cell (NMatrix *mat, NCell *cell, unsigned int seq);
static unsigned int n_matrix_can_merge_cells (NMatrix *mat, NCell *a, NCell *b);
static void n_matrix_merge_cells (NMatrix *mat, NCell *cell, NCell *other);
static NCell *n_matrix_get_seq_cell (NMatrix *mat, unsigned int seq, int pos);

static int
n_matrix_compare_cells (NMatrix *mat, NCell *lhs, NCell *rhs)
{
  unsigned int i;
  for (i = 0; i < mat->n_seqs; i++) {
    if (lhs->links[i].pos == UNKNOWN) continue;
    if (rhs->links[i].pos == UNKNOWN) continue;
    if (lhs->links[i].pos == BEFORE) {
      if (rhs->links[i].pos == BEFORE) continue;
      if (rhs->links[i].pos == AFTER) return -1;
      return -1;
    }
    if (lhs->links[i].pos == AFTER) {
      if (rhs->links[i].pos == AFTER) continue;
      if (rhs->links[i].pos == BEFORE) return 1;
      return 1;
    }
    if (rhs->links[i].pos == BEFORE) return 1;
    if (rhs->links[i].pos == AFTER) return -1;
    if (lhs->links[i].pos < rhs->links[i].pos) return -1;
    if (lhs->links[i].pos > rhs->links[i].pos) return 1;
  }
  return 0;
}

static void
n_matrix_link_cell (NMatrix *mat, NCell *cell, unsigned int seq, unsigned int pos)
{
  NCell *prev, *cur;
  if (cell->links[seq].pos != UNKNOWN) {
    fprintf (stderr, "n_matrix_link_cell: Cell already linked at seq %u pos %u\n", seq, pos);
  }
  cell->links[seq].pos = pos;
  prev = NULL;
  for (cur = mat->seqs[seq]->pos[pos].cells; cur; cur = cur->links[seq].next) {
    if (n_matrix_compare_cells (mat, cell, cur) < 0) break;
    prev = cur;
  }
  cell->links[seq].next = cur;
  if (cur) cur->links[seq].prev = cell;
  cell->links[seq].prev = prev;
  if (prev) {
    prev->links[seq].next = cell;
  } else {
    mat->seqs[seq]->pos[pos].cells = cell;
  }
}

static void
n_matrix_unlink_cell (NMatrix *mat, NCell *cell, unsigned int seq)
{
  NCell *prev, *cur, *next;
  unsigned int pos = cell->links[seq].pos;
  if (pos == UNKNOWN) {
    fprintf (stderr, "n_matrix_unlink_cell; Cell not linked at seq %u\n", seq);
  }
  prev = NULL;
  for (cur = mat->seqs[seq]->pos[pos].cells; cur; cur = cur->links[seq].next) {
    if (cur == cell) break;
    prev = cur;
  }
  next = cur->links[seq].next;
  if (prev) {
    prev->links[seq].next = next;
  } else {
    mat->seqs[seq]->pos[pos].cells = next;
  }
  if (next) next->links[seq].prev = prev;
  cur->links[seq].prev = NULL;
  cur->links[seq].next = NULL;
  cell->links[seq].pos = UNKNOWN;
}

static unsigned int
n_matrix_can_merge_cells (NMatrix *mat, NCell *a, NCell *b)
{
  unsigned int i;
  for (i = 0; i < mat->n_seqs; i++) {
    if ((a->links[i].pos != UNKNOWN) && (b->links[i].pos != UNKNOWN)) return 0;
  }
  return 1;
}

static void
n_matrix_merge_cells (NMatrix *mat, NCell *cell, NCell *other)
{
  unsigned int i;
  for (i = 0; i < mat->n_seqs; i++) {
    unsigned int pos = other->links[i].pos;
    if (pos != UNKNOWN) {
      n_matrix_unlink_cell (mat, other, i);
      n_matrix_link_cell (mat, cell, i, pos);
    }
  }
  n_matrix_free_cell (mat, other);
}

static NCell *
n_matrix_get_seq_cell (NMatrix *mat, unsigned int seq, int pos)
{
  if (pos < 0) return NULL;
  if (pos >= mat->seqs[seq]->len) return NULL;
  return mat->seqs[seq]->pos[pos].cells;
}

static int
compare_kmers (const NKMer *lhs, const NKMer *rhs)
{
  if (lhs->value < rhs->value) return -1;
  if (lhs->value > rhs->value) return 1;
  if (lhs->seq < rhs->seq) return -1;
  if (lhs->seq > rhs->seq) return 1;
  if (lhs->pos < rhs->pos) return -1;
  if (lhs->pos > rhs->pos) return 1;
  return 0;
}

static int
compare_kmer_counts (const NKMer *lhs, const NKMer *rhs)
{
  if (lhs->count < rhs->count) return 1;
  if (lhs->count > rhs->count) return -1;
  return 0;
}

int
main (int argc, const char *argv[])
{
  unsigned int i;
  const char *dbb = NULL;
  const char *kmers[1024];
  unsigned int nkmers = 0;
  unsigned int print_reads = 0, print_kmers = 0, print_chains = 0;
  unsigned long long word, rword;
  unsigned int code, node_idx, node_kmer, kmer_idx;

  KMerDB db;
  SeqFile *files;
  unsigned int nreads;
  Read reads[1024];
  char *seqs[1024];
  NMatrix *mat;
  NCell *cell;
    
  for (i = 1; i < argc; i++) {
    if (!strcmp (argv[i], "-dbb")) {
      i += 1;
      if (i >= argc) exit (1);
      dbb = argv[i];
    } else if (!strcmp (argv[i], "--print_reads")) {
      print_reads = 1;
    } else if (!strcmp (argv[i], "--print_kmers")) {
      print_kmers = 1;
    } else if (!strcmp (argv[i], "--print_chains")) {
      print_chains = 1;
    } else  {
      kmers[nkmers++] = argv[i];
    }
  }

  /* Read binary database */
  const unsigned char *cdata;
  unsigned long long csize;

  if (debug) fprintf (stderr, "Loading binary database %s\n", dbb);
  cdata = gt4_mmap (dbb, &csize);
  if (!cdata) {
    fprintf (stderr, "Cannot mmap %s\n", dbb);
    exit (1);
  }
  scout_mmap (cdata, csize);
  if (!read_database_from_binary (&db, cdata, csize)) {
    fprintf (stderr, "Cannot read binary database %s\n", dbb);
    exit (1);
  }
  if (debug) fprintf (stderr, "Finished loading binary database (index = %u)\n", db.index.read_blocks != NULL);

  files = (SeqFile *) malloc (db.index.n_files * sizeof (SeqFile));
  for (i = 0; i < db.index.n_files; i++) {
    files[i].name = db.index.files[i];
  }

  /* Build a list of unique reads */
  nreads = 0;
  for (i = 0; i < nkmers; i++) {
    unsigned long long first_read;
    unsigned int num_reads, j, kmer_dir;
    word = string_to_word (kmers[i], strlen (kmers[i]));
    rword = get_reverse_complement (word, strlen (kmers[i]));
    if (rword < word) word = rword;
    code = trie_lookup (&db.trie, word);
    if (!code) {
      fprintf (stderr, "No such kmer\n");
      exit (0);
    }
    kmer_dir = ((code & 0x80000000) != 0);
    if (debug > 0) fprintf (stderr, "Kmer %s word %llu code %u\n", kmers[i], word, code);
    code &= 0x7fffffff;
    node_idx = (code >> db.kmer_bits) - 1;
    node_kmer = code & ((1 << db.kmer_bits) - 1);
    kmer_idx = db.nodes[node_idx].kmers + node_kmer;
    if (debug > 0) fprintf (stderr, "Node %u kmer %u idx %u dir %u\n", node_idx, node_kmer, kmer_idx, kmer_dir);
  
    if (debug > 1) print_db_reads (&db.index, files, kmer_idx, kmer_dir, stderr);

    first_read = gt4_index_get_kmer_info (&db.index, kmer_idx, &num_reads);
    if (debug > 0) fprintf (stderr, "Num reads %u\n", num_reads);
    for (j = 0; j < num_reads; j++) {
      unsigned long long kmer_pos, name_pos;
      unsigned int file_idx, dir;
      unsigned int k;
      kmer_pos = gt4_index_get_read_info (&db.index, first_read + j, &file_idx, &name_pos, &dir);
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
  }
  if (debug) fprintf (stderr, "Total %u unique reads\n", nreads);
  /* Chenge directionality of reads if needed */
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
  mat = n_matrix_new (nreads, (const char **) seqs);
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
  /* Fill matrix */
  fprintf (stderr, "Building matrix\n");
  /* Iterate over all unique kmers */
  for (i = 0; i < mat->n_unique_kmers; i++) {
    NCell *cells[WORDLEN] = { 0 };
    NCell *prev = NULL;
    unsigned int j;
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
        if (!mat->seqs[seq_idx]->pos[pos + k].cells) {
          /* Link cell from current kmer to seq/pos */
          n_matrix_link_cell (mat, cells[k], seq_idx, pos + k);
        }
#if 0
        if (k > 0) {
          /* Chain cells on this sequence */
          if (prev->links[seq_idx].next && (prev->links[seq_idx].next != cells[k])) {
            fprintf (stderr, "Inconsistent next link for seq %u\n", seq_idx);
            fprintf (stderr, "Pos by current kmer: %d\n", cells[k]->links[seq_idx].pos);
            fprintf (stderr, "Pos by current chain: %d\n", prev->links[seq_idx].next->links[seq_idx].pos);
          } else if (cells[k]->links[seq_idx].prev && (cells[k]->links[seq_idx].prev != prev)) {
            fprintf (stderr, "Inconsistent prev link for seq %u\n", seq_idx);
            fprintf (stderr, "Pos by current kmer: %d\n", prev->links[seq_idx].pos);
            fprintf (stderr, "Pos by current chain: %d\n", cells[k]->links[seq_idx].prev->links[seq_idx].pos);
          } else {
            prev->links[seq_idx].next = cells[k];
            cells[k]->links[seq_idx].prev = prev;
          }
        }
#endif
        prev = cells[k];
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

  fprintf (stderr, "Finding longest alignment\n");
  /* Get longest alignment */
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
  if (debug) fprintf (stderr, "Alignment length %u\n", alen);
  
  /* Create alignment */
  NCell **alignment = (NCell **) malloc (alen * sizeof (NCell *));
  unsigned int apos = 0;
  for (cell = first_cell; cell; cell = cell->next) {
    alignment[apos++] = cell;
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
              unsigned int l;
              fprintf (stderr, "Cannot merge: alignment pos %u seq %u pos %u\n", i, j, s_p);
              for (l = 0; l < mat->n_seqs; l++) {
                fprintf (stderr, "%c", 'A' + alignment[k]->links[l].pos);
              }
              fprintf (stderr, "\n");
              for (l = 0; l < mat->n_seqs; l++) {
                fprintf (stderr, "%c", 'A' + cc->links[l].pos);
              }
              fprintf (stderr, "\n");
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
  }

  for (i = 0; i < mat->n_seqs; i++) {
    unsigned int j;
    for (j = 0; j < alen; j++) {
      NCell *cell = alignment[j];
      char c;
      if (cell->links[i].pos == BEFORE) {
        c = '<';
      } else if (cell->links[i].pos == AFTER) {
        c = '>';
      } else if (cell->links[i].pos == UNKNOWN) {
        c = '?';
      } else {
        char *n = "ACGTN.,";
        c = n[mat->seqs[i]->pos[cell->links[i].pos].nucl];
      }
      fprintf (stderr, "%c", c);
    }
    fprintf (stderr, "\n");
  }

  int *align = (int *) malloc (mat->n_seqs * alen * sizeof (int));
  for (i = 0; i < mat->n_seqs * alen; i++) align[i] = -1;
  for (i = 0; i < alen; i++) {
    unsigned int j;
    for (j = 0; j < mat->n_seqs; j++) {
      if (alignment[i]->links[j].pos >= 0) align[j * alen + i] = alignment[i]->links[j].pos;
    }
  }

  /* Printout */
  unsigned int max_s = 0;
  unsigned int s[4096];
  float q[4096];
  for (i = 0; i < alen; i++) {
    unsigned int j;
    unsigned int s_n[4] = { 0 };
    for (j = 0; j < mat->n_seqs; j++) {
      if (align[j * alen + i] < 0) continue;
      if (align[j * alen + i] >= mat->seqs[j]->len) continue;
      unsigned int nucl = mat->seqs[j]->pos[align[j * alen + i]].nucl;
      if (nucl > T) continue;
      s_n[nucl] += 1;
    }
    s[i] = s_n[A] + s_n[C] + s_n[G] + s_n[T];
    if (s[i] > max_s) max_s = s[i];
    unsigned int max = A;
    for (j = A; j <= T; j++) {
      if (s_n[j] > s_n[max]) max = j;
    }
    if (s[i] != 0) {
      q[i] = (float) s_n[max] / s[i];
    } else {
      q[i] = 0;
    }
  }
  for (i = 0; i < alen; i++) {
    fprintf (stdout, "%c", '1' + (unsigned int) (8 * s[i] / max_s));
  }
  fprintf (stdout, "\n");
  for (i = 0; i < alen; i++) {
    fprintf (stdout, "%c", '1' + (unsigned int) (8 * q[i]));
  }
  fprintf (stdout, "\n");
  for (i = 0; i < mat->n_seqs; i++) {
    unsigned int j;
    for (j = 0; j < alen; j++) {
      if (align[i * alen + j] < 0) {
        fprintf (stdout, " ");
      } else if (align[i * alen + j] < mat->seqs[i]->len) {
        unsigned int nucl = mat->seqs[i]->pos[align[i * alen + j]].nucl;
        fprintf (stdout, "%c", n2c[nucl]);
      } else {
        break;
      }
    }
    fprintf (stdout, "\n");
  }
  
  return 0;
}

static NSeq *
n_seq_new (const char *str)
{
  NSeq *seq;
  unsigned long long len;
  unsigned int i;
  len = strlen (str);
  seq = (NSeq *) malloc (sizeof (NSeq) + (len - 1) * sizeof (NPos));
  seq->len = len;
  seq->str = str;
  for (i = 0; i < len; i++) {
    static unsigned int *c2n = NULL;
    if (!c2n) {
      unsigned int j;
      c2n = (unsigned int *) malloc (256 * 4);
      for (j = 0; j < 256; j++) c2n[j] = N;
      c2n['a'] = c2n['A'] = A;
      c2n['c'] = c2n['C'] = C;
      c2n['g'] = c2n['G'] = G;
      c2n['t'] = c2n['T'] = c2n['u'] = c2n['U'] = T;
    }
    seq->pos[i].nucl = c2n[(unsigned int) str[i]];
    seq->pos[i].has_kmer = 0;
    seq->pos[i].kmer = 0;
    seq->pos[i].cells = NULL;
  }
  for (i = 0; i <= (len - WORDLEN); i++) {
    seq->pos[i].has_kmer = n_seq_get_kmer (seq, i, &seq->pos[i].kmer);
  }
  return seq;
}

static unsigned int
n_seq_get_kmer (NSeq *seq, unsigned int pos, unsigned long long *kmer)
{
  unsigned long long val = 0;
  unsigned int i;
  for (i = 0; i < WORDLEN; i++) {
    if (seq->pos[pos + i].nucl >= N) return 0;
    val = val << 2;
    val |= seq->pos[pos + i].nucl;
  }
  *kmer = val;
  return 1;
}

static NMatrix *
n_matrix_new (unsigned int n_seqs, const char *seqs[])
{
  NMatrix *mat;
  unsigned int i, size_kmers = 256;
  mat = (NMatrix *) malloc (sizeof (NMatrix) + (n_seqs - 1) * sizeof (NSeq *));
  mat->n_seqs = n_seqs;
  mat->n_kmers = 0;
  mat->kmers = (NKMer *) malloc (size_kmers * sizeof (NKMer));
  mat->n_unique_kmers = 0;
  mat->unique_kmers = NULL;
  mat->cells = NULL;
  mat->free_cells = NULL;
  for (i = 0; i < n_seqs; i++) {
    unsigned int j;
    mat->seqs[i] = n_seq_new (seqs[i]);
    for (j = 0; j < mat->seqs[i]->len; j++) {
      if (mat->seqs[i]->pos[j].has_kmer) {
        if (mat->n_kmers >= size_kmers) {
          size_kmers = size_kmers << 1;
          mat->kmers = (NKMer *) realloc (mat->kmers, size_kmers * sizeof (NKMer));
        }
        mat->kmers[mat->n_kmers].value = mat->seqs[i]->pos[j].kmer;
        mat->kmers[mat->n_kmers].seq = i;
        mat->kmers[mat->n_kmers].pos = j;
        mat->kmers[mat->n_kmers].count = 1;
        mat->n_kmers += 1;
      }
    }
  }
  /* Sort kmers */
  qsort (mat->kmers, mat->n_kmers, sizeof (NKMer), (int (*) (const void *, const void *)) compare_kmers);
  /* Fill unique kmers */
  mat->unique_kmers = (NKMer *) malloc (mat->n_kmers * sizeof (NKMer));
  if (mat->n_kmers) {
    i = 0;
    while (i < mat->n_kmers) {
      unsigned int j;
      mat->unique_kmers[mat->n_unique_kmers] = mat->kmers[i];
      for (j = i + 1; (j < mat->n_kmers) && (mat->kmers[j].value == mat->unique_kmers[mat->n_unique_kmers].value); j++) {
        mat->unique_kmers[mat->n_unique_kmers].count += 1;
      }
      mat->n_unique_kmers += 1;
      i = j;
    }
    qsort (mat->unique_kmers, mat->n_unique_kmers, sizeof (NKMer), (int (*) (const void *, const void *)) compare_kmer_counts);
  }
  
  return mat;
}

static NCell *
n_matrix_new_cell (NMatrix *mat)
{
  NCell *cell;
  unsigned int i;
  if (mat->free_cells) {
    cell = mat->free_cells;
    mat->free_cells = cell->next_allocated;
  } else {
    cell = (NCell *) malloc (sizeof (NCell) + (mat->n_seqs - 1) * sizeof (NLink));
  }
  cell->prev = NULL;
  cell->next = NULL;
  cell->count = 0;
  for (i = 0; i < mat->n_seqs; i++) {
    cell->links[i].pos = UNKNOWN;
    cell->links[i].prev = NULL;
    cell->links[i].next = NULL;
  }
  cell->next_allocated = mat->cells;
  mat->cells = cell;
  return cell;
}

static void
n_matrix_free_cell (NMatrix *mat, NCell *cell)
{
  if (cell == mat->cells) {
    mat->cells = cell->next_allocated;
  } else {
    NCell *prev;
    prev = mat->cells;
    while (prev->next_allocated != cell) prev = prev->next_allocated;
    prev->next_allocated = cell->next_allocated;
  }
  cell->next_allocated = mat->free_cells;
  mat->free_cells = cell;
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
