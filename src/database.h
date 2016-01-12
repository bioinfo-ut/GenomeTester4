#ifndef __DATABASE_H__
#define __DATABASE_H__

#include <stdio.h>

#include "trie.h"

typedef struct _Node Node;

struct _Node {
  /* Index in db name table */
  unsigned int name;
  /* Index in db kmer table */
  unsigned int kmers;
  unsigned int nkmers;
};

typedef struct _KMerDB KMerDB;

struct _KMerDB {
  unsigned int wordsize;
  unsigned int node_bits;
  unsigned int kmer_bits;
  unsigned long long n_nodes;
  unsigned long long n_kmers;
  unsigned long long names_size;
  /* Table of nodes */
  Node *nodes;
  /* Table of kmer counts */
  unsigned short *kmers;
  /* Table of names */
  char *names;
  /* Trie mapping kmers to nodes/counts */
  Trie trie;
};

/* Return number of nodes successfully read */
unsigned int read_db_from_text (KMerDB *db, const unsigned char *cdata, unsigned long long csize, unsigned int max_kmers_per_node);

/*
 * Binary representation
 *
 * 0  : wordsize (4)
 * 4  : node_bits (4)
 * 8  : kmer_bits (4)
 * 12 : dummy (4)
 * 16 : n_nodes (8)
 * 24 : n_kmers (8)
 * 32 : names_size (8)
 * 40 : nodes_blocksize (8)
 * 48 : Nodes
 *    : kmers_blocksize (8)
 *    : Kmers
 *    : names_blocksize (8)
 *    : Names
 *    : Trie
*/

/* Return number of bytes written */
unsigned int write_db_to_file (KMerDB *db, FILE *ofs);

unsigned int read_database_from_binary (KMerDB *db, const unsigned char *cdata, unsigned long long csize);

#endif
