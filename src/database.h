#ifndef __GM4_DATABASE_H__
#define __GM4_DATABASE_H__

#include <stdio.h>

#include "index.h"
#include "trie.h"

#ifndef __DATABASE_C__
extern unsigned int db_debug;
#endif

typedef struct _Node Node;

struct _Node {
  /* Index in db name table */
  unsigned int name;
  /* Index in db kmer table */
  unsigned int kmers;
  unsigned int nkmers;
};

typedef struct _GT4GmerDB GT4GmerDB;

struct _GT4GmerDB {
  unsigned int major;
  unsigned int minor;
  unsigned int wordsize;
  unsigned int node_bits;
  unsigned int kmer_bits;
  unsigned int count_bits;
  unsigned long long n_nodes;
  unsigned long long n_kmers;
  unsigned long long names_size;
  /* Table of nodes */
  Node *nodes;
  /* Table of kmer counts */
  union {
    unsigned short *kmers_16;
    unsigned int *kmers_32;
  };
  /* Table of names */
  char *names;
  /* Trie mapping kmers to nodes/counts */
  Trie trie;
  /* Read index */
  GT4Index index;
};

/* Reads */

typedef struct _Read Read;

struct _Read {
  /* Subsequence index */
  unsigned int subseq;
  /* Start of kmer relative to the start of name */
  unsigned int kmer_pos : 18;
  /* File number */
  unsigned int source_idx : 12;
  /* Direction */
  unsigned int dir : 1;
};

/* Read list */

typedef struct _ReadList ReadList;

struct _ReadList {
  ReadList *next;
  Read read;
};

ReadList *gm4_read_list_new (void);

/* Return number of nodes successfully read */
GT4GmerDB *gt4_gmer_db_new_from_text (const unsigned char *cdata, unsigned long long csize, unsigned int max_kmers_per_node, unsigned int count_bits);

/*
 * Binary representation
 *
 * 0  : "GMDB" (4)
 * 4  : major (2)
 * 6  : minor (2)
 * 8  : wordsize (4)
 * 12 : node_bits (4)
 * 16 : kmer_bits (4)
 * 20 : count_bits (4)
 * 24 : n_nodes (8)
 * 32 : n_kmers (8)
 * 40 : names_size (8)
 * 48 : nodes_start (8)
 * 56 : kmers_start (8)
 * 64 : names_start (8)
 * 72 : trie_start (8)
 * 80 : index start (8)
 *    : nodes_blocksize (8)
 *    : Nodes
 *    : kmers_blocksize (8)
 *    : Kmers
 *    : names_blocksize (8)
 *    : Names
 *    : trie_blocksize (8)
 *    : Trie
 *    : index_blocksize (8)
 *    : Index
*/

/* Return number of bytes written */
unsigned int write_db_to_file (GT4GmerDB *db, FILE *ofs, unsigned int kmers);
unsigned int write_db_to_file_with_reads_callback (GT4GmerDB *db, FILE *ofs, unsigned int kmers, unsigned long long (*write_reads) (GT4Index *index, FILE *ofs, void *data), void *data);

GT4GmerDB *gt4_gmer_db_new_from_binary (const unsigned char *cdata, unsigned long long csize);

void gt4_db_clear_index (GT4GmerDB *db);

/* Debug */
void gt4_db_dump (GT4GmerDB *db, FILE *ofs);

#endif
