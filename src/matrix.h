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

unsigned int c2n (unsigned char c);

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
  unsigned int has_kmer : 1;
  unsigned int non_unique_kmer : 1;
  unsigned long long kmer;
  /* Sorted list of cells aligned with this position */
  NCell *cells;
};

struct _NSeq {
  const char *str;
  /* K-mer length */
  unsigned int wlen;
  /* Sequence length */
  unsigned int len;
  NPos *pos;
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
  int score;
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
void n_seq_delete (NSeq *seq);
unsigned int n_seq_get_kmer (NSeq *seq, unsigned int pos, unsigned long long *kmer);
unsigned int n_seq_get_kmer_unique_pos (NSeq *seq, unsigned long long word);

NMatrix *n_matrix_new (unsigned int n_seqs, const char *seqs[], unsigned int wlen);
void n_matrix_delete (NMatrix *mat);
unsigned int n_matrix_get_kmer_first_index (NMatrix *mat, unsigned long long value);
unsigned int n_matrix_get_kmer_unique_index (NMatrix *mat, unsigned long long value);

NCell *n_matrix_new_cell (NMatrix *mat);
void n_matrix_free_cell (NMatrix *mat, NCell *cell);
int n_matrix_compare_cells (NMatrix *mat, NCell *lhs, NCell *rhs);

/* Link cell to sequences (where pos >= 0), create new cell if needed */
NCell *n_matrix_link_sequences (NMatrix *mat, unsigned int seqs[], unsigned int positions[], unsigned int nseqs);
/* Recalculate all scores, return cell with biggest score */
NCell *n_matrix_calculate_scores (NMatrix *mat);

/* Link/unlink specific cell to sequence */
void n_matrix_link_cell (NMatrix *mat, NCell *cell, unsigned int seq, unsigned int pos);
void n_matrix_unlink_cell (NMatrix *mat, NCell *cell, unsigned int seq);

unsigned int n_matrix_can_merge_cells (NMatrix *mat, NCell *a, NCell *b);
void n_matrix_merge_cells (NMatrix *mat, NCell *cell, NCell *other);
NCell *n_matrix_get_seq_cell (NMatrix *mat, unsigned int seq, int pos);

int n_matrix_calculate_cell_score (NMatrix *mat, NCell *cell);

#endif
