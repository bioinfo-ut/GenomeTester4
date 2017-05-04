#define __GM4_INDEX_C__

#include <malloc.h>
#include <string.h>

#include "index.h"

unsigned long long
gt4_index_get_kmer_info (GT4Index *index, unsigned int kmer, unsigned int *num_reads)
{
  *num_reads = index->read_blocks[kmer] & 0xffffff;
  return index->read_blocks[kmer] >> 24;
}

unsigned long long
gt4_index_get_read_info (GT4Index *index, unsigned long long read, unsigned int *file_idx, unsigned long long *name_pos, unsigned int *dir)
{
  unsigned long long code;
  code = read = index->reads[read];
  *name_pos = (code >> index->nbits_kmer) & ((1ULL << index->nbits_npos) - 1);
  *file_idx = (code >> (index->nbits_npos + index->nbits_kmer)) & ((1ULL << index->nbits_file) - 1);
  *dir = (code >> (index->nbits_npos + index->nbits_kmer + index->nbits_file)) & 1;
  return code & ((1ULL << index->nbits_kmer) - 1);
}

unsigned int
gt4_index_init_from_data (GT4Index *index, const unsigned char *cdata, unsigned long long csize, unsigned long long n_kmers)
{
  unsigned long long cpos = 0, files_start, blocks_start, reads_start;
  unsigned int i;

  /* Bitsizes */
  memcpy (&index->nbits_file, cdata + cpos, 4);
  memcpy (&index->nbits_npos, cdata + cpos + 4, 4);
  memcpy (&index->nbits_kmer, cdata + cpos + 8, 4);
  cpos += 12;
  /* Sizes */
  memcpy (&index->n_files, cdata + cpos, 4);
  memcpy (&index->n_kmers, cdata + cpos + 4, 8);
  memcpy (&index->n_reads, cdata + cpos + 12, 8);
  cpos += 20;
  /* Starts */
  memcpy (&files_start, cdata + cpos, 8);
  memcpy (&blocks_start, cdata + cpos + 8, 8);
  memcpy (&reads_start, cdata + cpos + 16, 8);
  cpos += 24;
  /* Files */
  cpos = files_start;
  index->files = (char **) malloc (index->n_files * sizeof (char *));
  for (i = 0; i < index->n_files; i++) {
    index->files[i] = (char *) cdata + cpos;
    cpos += strlen ((const char *) cdata + cpos);
    cpos += 1;
  }
  /* Blocks */
  cpos = blocks_start;
  index->read_blocks = (unsigned long long *) (cdata + cpos);
  /* Reads */
  cpos = reads_start;
  index->reads = (unsigned long long *) (cdata + cpos);

  return 1;
}

unsigned int
gt4_index_write (GT4Index *index, FILE *ofs, unsigned long long n_kmers)
{
  unsigned long long fpos;
  unsigned long long written = 0, files_start, blocks_start, reads_start;
  unsigned int i;

  fpos = ftello (ofs);

  /* Bitsizes */
  fwrite (&index->nbits_file, 4, 1, ofs);
  fwrite (&index->nbits_npos, 4, 1, ofs);
  fwrite (&index->nbits_kmer, 4, 1, ofs);
  written += 12;
  /* Sizes */
  fwrite (&index->n_files, 4, 1, ofs);
  fwrite (&index->n_kmers, 8, 1, ofs);
  fwrite (&index->n_reads, 8, 1, ofs);
  written += 20;
  /* Starts */
  fseek (ofs, fpos + written + 24, SEEK_SET);
  written += 24;
  /* Files */
  files_start = written;
  fseek (ofs, fpos + files_start, SEEK_SET);
  for (i = 0; i < index->n_files; i++) {
    unsigned int len = (unsigned int) strlen (index->files[i]) + 1;
    fwrite (index->files[i], len, 1, ofs);
    written += len;
  }
  written = (written + 15) & 0xfffffffffffffff0;
  /* Blocks */
  blocks_start = written;
  fseek (ofs, fpos + blocks_start, SEEK_SET);
  if (index->read_blocks) {
    fwrite (index->read_blocks, 8, n_kmers, ofs);
    written += n_kmers * 8;
  }
  written = (written + 15) & 0xfffffffffffffff0;
  /* Reads */
  reads_start = written;
  fseek (ofs, fpos + reads_start, SEEK_SET);
  if (index->reads) {
    fwrite (index->reads, 8, index->n_reads, ofs);
    written += index->n_reads * 8;
  }
  written = (written + 15) & 0xfffffffffffffff0;
  /* Starts */
  fseek (ofs, fpos + 32, SEEK_SET);
  fwrite (&files_start, 8, 1, ofs);
  fwrite (&blocks_start, 8, 1, ofs);
  fwrite (&reads_start, 8, 1, ofs);

  fseek (ofs, fpos + written, SEEK_SET);

  return written;
}

