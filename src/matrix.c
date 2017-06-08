#define __GT4_MATRIX_C__

#include <malloc.h>
#include <string.h>
#include <stdlib.h>

#include <matrix.h>

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

NSeq *
n_seq_new (const char *str, unsigned int wlen)
{
  NSeq *seq;
  unsigned long long len;
  unsigned int i;
  len = strlen (str);
  seq = (NSeq *) malloc (sizeof (NSeq) + (len - 1) * sizeof (NPos));
  seq->wlen = wlen;
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
  for (i = 0; i <= (len - wlen); i++) {
    seq->pos[i].has_kmer = n_seq_get_kmer (seq, i, &seq->pos[i].kmer);
  }
  return seq;
}

unsigned int
n_seq_get_kmer (NSeq *seq, unsigned int pos, unsigned long long *kmer)
{
  unsigned long long val = 0;
  unsigned int i;
  for (i = 0; i < seq->wlen; i++) {
    if (seq->pos[pos + i].nucl >= N) return 0;
    val = val << 2;
    val |= seq->pos[pos + i].nucl;
  }
  *kmer = val;
  return 1;
}

NMatrix *
n_matrix_new (unsigned int n_seqs, const char *seqs[], unsigned int wlen)
{
  NMatrix *mat;
  unsigned int i, size_kmers = 256;
  mat = (NMatrix *) malloc (sizeof (NMatrix) + (n_seqs - 1) * sizeof (NSeq *));
  mat->wlen = wlen;
  mat->n_seqs = n_seqs;
  mat->n_kmers = 0;
  mat->kmers = (NKMer *) malloc (size_kmers * sizeof (NKMer));
  mat->n_unique_kmers = 0;
  mat->unique_kmers = NULL;
  mat->cells = NULL;
  mat->free_cells = NULL;
  for (i = 0; i < n_seqs; i++) {
    unsigned int j;
    mat->seqs[i] = n_seq_new (seqs[i], wlen);
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
      unsigned int j, last_seq;
      unsigned int duplicate = 0;
      mat->unique_kmers[mat->n_unique_kmers] = mat->kmers[i];
      last_seq = mat->kmers[i].seq;
      for (j = i + 1; (j < mat->n_kmers) && (mat->kmers[j].value == mat->unique_kmers[mat->n_unique_kmers].value); j++) {
        if (mat->kmers[j].seq == last_seq) duplicate = 1;
        mat->unique_kmers[mat->n_unique_kmers].count += 1;
        last_seq = mat->kmers[j].seq;
      }
      if (!duplicate) mat->n_unique_kmers += 1;
      i = j;
    }
    qsort (mat->unique_kmers, mat->n_unique_kmers, sizeof (NKMer), (int (*) (const void *, const void *)) compare_kmer_counts);
  }
  
  return mat;
}

NCell *
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
  memset (cell, 0, sizeof (NCell));
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

void
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

int
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

void
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

void
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

unsigned int
n_matrix_can_merge_cells (NMatrix *mat, NCell *a, NCell *b)
{
  unsigned int i;
  for (i = 0; i < mat->n_seqs; i++) {
    if ((a->links[i].pos != UNKNOWN) && (b->links[i].pos != UNKNOWN)) return 0;
  }
  /* fixme */
  if (a->ref_pos && b->ref_pos && (a->ref_pos != b->ref_pos)) return 0;
  return 1;
}

void
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
  if (other->ref_pos && !cell->ref_pos) {
    cell->file_idx = other->file_idx;
    cell->ref_pos = other->ref_pos;
  }
  n_matrix_free_cell (mat, other);
}

NCell *
n_matrix_get_seq_cell (NMatrix *mat, unsigned int seq, int pos)
{
  if (pos < 0) return NULL;
  if (pos >= mat->seqs[seq]->len) return NULL;
  return mat->seqs[seq]->pos[pos].cells;
}

