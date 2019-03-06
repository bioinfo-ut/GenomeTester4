#ifndef __GT4_INDEX_H__
#define __GT4_INDEX_H__

#include <stdio.h>

typedef struct _GT4Index GT4Index;

struct _GT4Index {
  unsigned int code;
  unsigned int version_major;
  unsigned int version_minor;
  unsigned int filler;
  /* Read layout */
  unsigned int nbits_file;
  unsigned int nbits_npos;
  unsigned int nbits_kmer;
  unsigned int n_files;
  unsigned long long n_kmers;
  unsigned long long n_reads;
  /* Read files */
  char **files;
  /* Read starts and counts per kmer */
  /* Bits 63...24 - start, 23...0 count */
  unsigned long long *read_blocks;
  /* Read data */
  unsigned long long *reads;
};

/* Return the first read index */
unsigned long long gt4_index_get_kmer_info (GT4Index *index, unsigned int kmer, unsigned int *num_reads);
/* Return kmer position in sequence */
unsigned long long gt4_index_get_read_info (GT4Index *index, unsigned long long read, unsigned int *file_idx, unsigned long long *name_pos, unsigned int *dir);

/*
 * Index layout
 *
 * 0  : file_idx_bits (4)
 * 4  : name_pos_bits (4)
 * 8  : kmer_pos_bits (4)
 * 12 : n_files (4)
 * 16 : n_kmers (8)
 * 24 : n_reads (8)
 * 32 : files_start (8)
 * 40 : blocks_start (8)
 * 48 : reads_start (8)
 *    : files
 *    : kmer_reads
 *    : reads
 */

unsigned int gt4_index_init_from_data (GT4Index *index, const unsigned char *cdata, unsigned long long csize, unsigned long long n_kmers, unsigned int compatibility_mode);
unsigned long long gt4_index_write (GT4Index *index, FILE *ofs, unsigned long long n_kmers);
unsigned long long gt4_index_write_with_reads_callback (GT4Index *index, FILE *ofs, unsigned long long n_kmers, unsigned long long (*write_reads) (GT4Index *index, FILE *ofs, void *data), void *data);

#endif
