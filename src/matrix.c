#define __GT4_MATRIX_C__

#include <assert.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>

#include <matrix.h>

unsigned int
c2n (unsigned char c)
{
  static unsigned int *_c2n = NULL;
  if (!_c2n) {
    unsigned int j;
    _c2n = (unsigned int *) malloc (256 * 4);
    for (j = 0; j < 256; j++) _c2n[j] = N;
    _c2n['a'] = _c2n['A'] = A;
    _c2n['c'] = _c2n['C'] = C;
    _c2n['g'] = _c2n['G'] = G;
    _c2n['t'] = _c2n['T'] = _c2n['u'] = _c2n['U'] = T;
    _c2n['-'] = GAP;
  }
  return _c2n[c];
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

NSeq *
n_seq_new_length (const char *str, unsigned int len, unsigned int wlen)
{
  NSeq *seq;
  unsigned int i, j;
  seq = (NSeq *) malloc (sizeof (NSeq));
  memset (seq, 0, sizeof (NSeq));
  seq->wlen = wlen;
  seq->len = len;
  seq->str = str;
  seq->pos = (NPos *) malloc (len * sizeof (NPos));
  memset (seq->pos, 0, len * sizeof (NPos));
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
  for (i = 0; (i + wlen) <= len; i++) {
    if (n_seq_get_kmer (seq, i, &seq->pos[i].kmer)) {
      seq->pos[i].has_kmer = 1;
      for (j = 0; j < i; j++) {
        if (seq->pos[j].has_kmer && (seq->pos[j].kmer == seq->pos[i].kmer)) {
          seq->pos[j].non_unique_kmer = 1;
          seq->pos[i].non_unique_kmer = 1;
        }
      }
    }
  }
  return seq;
}

NSeq *
n_seq_new (const char *str, unsigned int wlen)
{
  return n_seq_new_length (str, (unsigned int) strlen (str), wlen);
}

void
n_seq_delete (NSeq *seq)
{
  free (seq->pos);
  free (seq);
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

unsigned int
n_seq_get_kmer_unique_pos (NSeq *seq, unsigned long long word)
{
  unsigned int i;
  for (i = 0; i < seq->len; i++) {
    if (seq->pos[i].has_kmer && !seq->pos[i].non_unique_kmer && (seq->pos[i].kmer == word)) break;
  }
  return i;
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
      if (mat->seqs[i]->pos[j].has_kmer && !mat->seqs[i]->pos[j].non_unique_kmer) {
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
      unsigned long long value;
      mat->unique_kmers[mat->n_unique_kmers] = mat->kmers[i];
      value = mat->kmers[i].value;
      last_seq = mat->kmers[i].seq;
      for (j = i + 1; (j < mat->n_kmers) && (mat->kmers[j].value == value); j++) {
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

void
n_matrix_delete (NMatrix *mat)
{
  unsigned int i;
  while (mat->cells) n_matrix_free_cell (mat, mat->cells);
  while (mat->free_cells) {
    NCell *cell = mat->free_cells;
    mat->free_cells = cell->next_allocated;
    free (cell);
  }
  for (i = 0; i < mat->n_seqs; i++) n_seq_delete (mat->seqs[i]);
  free (mat->kmers);
  free (mat->unique_kmers);
  free (mat);
}

unsigned int
n_matrix_get_kmer_first_index (NMatrix *mat, unsigned long long value)
{
  unsigned int i;
  for (i = 0; i < mat->n_kmers; i++) {
    if (mat->kmers[i].value == value) break;
  }
  return i;
}

unsigned int
n_matrix_get_kmer_unique_index (NMatrix *mat, unsigned long long value)
{
  unsigned int i;
  for (i = 0; i < mat->n_unique_kmers; i++) {
    if (mat->unique_kmers[i].value == value) break;
  }
  return i;
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
  cell->score = 0;
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

NCell *
n_matrix_link_sequences (NMatrix *mat, unsigned int seqs[], unsigned int positions[], unsigned int nseqs)
{
  unsigned int best_count = 0, i;
  NCell *best_cell = NULL;
  for (i = 0; i < nseqs; i++) {
    NCell *cell;
    for (cell = mat->seqs[seqs[i]]->pos[positions[i]].cells; cell; cell = cell->links[seqs[i]].next) {
      unsigned int count = 0, j;
      assert (cell->links[seqs[i]].pos == positions[i]);
      for (j = 0; j < nseqs; j++) {
        if (cell->links[seqs[j]].pos < 0) continue;
        if (cell->links[seqs[j]].pos != positions[j]) break;
        count += 1;
      }
      if (j >= nseqs) {
        /* This cell can be linked */
        if (count > best_count) {
          best_count = count;
          best_cell = cell;
        }
      }
    }
  }
  if (!best_cell) best_cell = n_matrix_new_cell (mat);
  for (i = 0; i < nseqs; i++) {
    if (best_cell->links[seqs[i]].pos < 0) {
      n_matrix_link_cell (mat, best_cell, seqs[i], positions[i]);
    } else {
      assert (best_cell->links[seqs[i]].pos == positions[i]);
    }
  }
  return best_cell;
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

int
n_matrix_calculate_cell_score (NMatrix *mat, NCell *cell)
{
  unsigned int max, i;
  unsigned int count[7] = { 0 };
  int best;
  for (i = 0; i < mat->n_seqs; i++) {
    if (cell->links[i].pos < 0) continue;
    count[mat->seqs[i]->pos[cell->links[i].pos].nucl] += 1;
  }
  max = A;
  for (i = C; i <= T; i++) if (count[i] > count[max]) max = i;
  best = count[max];
  for (i = A; i <= T; i++) if (i != max) best -= count[i];
  best -= 2 * count[GAP];
  return best;
}

NCell *
n_matrix_calculate_scores (NMatrix *mat)
{
  NCell *cell, *best_cell = NULL;
  for (cell = mat->cells; cell; cell = cell->next_allocated) {
    cell->score = n_matrix_calculate_cell_score (mat, cell);
    if (!best_cell || (best_cell->score < cell->score)) best_cell = cell;
  }
  return best_cell;
}
