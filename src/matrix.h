#ifndef __GT4_MATRIX_H__
#define __GT4_MATRIX_H__

/*
 * A n-dimensional alignment matrix
 */

/* Nucleotide values */
#define A 0
#define C 1
#define G 2
#define T 3
#define N 4
#define GAP 5
#define NONE 6

#ifndef __GT4_MATRIX_C__
extern const char *n2c;
#else
const char *n2c = "ACGTN- ";
#endif

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
  unsigned int wlen;
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
  /* Consensus nucleotide */
  unsigned int consensus;
  /* Reference position */
  unsigned int file_idx;
  unsigned int ref_pos;
  NLink links[1];
};

struct _NMatrix {
  unsigned int wlen;
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

NSeq *n_seq_new (const char *str, unsigned int wlen);
unsigned int n_seq_get_kmer (NSeq *seq, unsigned int pos, unsigned long long *kmer);

NMatrix *n_matrix_new (unsigned int n_seqs, const char *seqs[], unsigned int wlen);

NCell *n_matrix_new_cell (NMatrix *mat);
void n_matrix_free_cell (NMatrix *mat, NCell *cell);
int n_matrix_compare_cells (NMatrix *mat, NCell *lhs, NCell *rhs);
void n_matrix_link_cell (NMatrix *mat, NCell *cell, unsigned int seq, unsigned int pos);
void n_matrix_unlink_cell (NMatrix *mat, NCell *cell, unsigned int seq);
unsigned int n_matrix_can_merge_cells (NMatrix *mat, NCell *a, NCell *b);
void n_matrix_merge_cells (NMatrix *mat, NCell *cell, NCell *other);
NCell *n_matrix_get_seq_cell (NMatrix *mat, unsigned int seq, int pos);

#endif
