#ifndef __GT4_TRIE_H__
#define __GT4_TRIE_H__

#include <stdio.h>
#include <pthread.h>

#ifndef __GT4_TRIE_C__
extern unsigned int gt4_trie_debug;
#endif

/* Block of Trie structures */

typedef struct _Trie Trie;
typedef struct _TrieNodeBranch TrieNodeBranch;
typedef unsigned long long TrieRef;

struct _TrieNodeBranch {
  unsigned long long _nbits_this : 5;
  unsigned long long nbits_children : 6;
  unsigned long long word : 26;
  TrieRef children[2];
};

/* Debug flags */
#define GT4_TRIE_COUNT_ALLOCATIONS 1

/* Allocation */
#define TRIE_BLOCK_BITS 30
#define TRIE_BLOCK_SIZE (1ULL << TRIE_BLOCK_BITS)

#define TRIE_BLOCK_FROM_REF(r) ((((r) >> 2) >> TRIE_BLOCK_BITS) & 0x3ff)
#define TRIE_INDEX_FROM_REF(r) (((r) >> 2) & (TRIE_BLOCK_SIZE - 1))
#define TRIE_ADDRESS_FROM_REF(t,r) ((t)->branches[TRIE_BLOCK_FROM_REF(r)] + TRIE_INDEX_FROM_REF(r))
#define TRIE_REF_FROM_ADDRESS(t,b,i) (((b) << TRIE_BLOCK_BITS) | (i)) << 2;

#define TYPE_BRANCH 0ULL
#define TYPE_KMER 1ULL

/* Ref */

#define REF_IS_EMPTY(r) (!(r))
#define REF_IS_KMER(r) (((r) & 1) == TYPE_KMER)
#define REF_IS_BRANCH(r) (((r) & 1) == TYPE_BRANCH)

#define BRANCH_FROM_REF(t,r) TRIE_ADDRESS_FROM_REF(t,r)

/* Kmer */

/* nbits:5 word:26 count:32 type:1 */
#define KMER_MAX_BITS 26

#define KMER_GET_NBITS(n) (((unsigned long long) (n) >> 59) & 0x1f)
#define KMER_GET_WORD(n) (((unsigned long long) (n) >> 33) & 0x3ffffff)
#define KMER_GET_COUNT(n) (((unsigned long long) (n) >> 1) & 0xffffffff)

#define MAKE_KMER(b,w,c) (((unsigned long long) (b) << 59) | ((unsigned long long) (w) << 33) | ((unsigned long long) (c) << 1) | TYPE_KMER)

/* Branch */
/* nbits_this:5 nbits_children:5 word:32 */
/* INVARIANT: branch max bits >= kmer max bits */

#define BRANCH_MAX_BITS_THIS 52
#define BRANCH_MAX_BITS_CHILDREN 64

#define BRANCH_GET_NBITS_THIS(n) (n)->_nbits_this
#define BRANCH_GET_NBITS_CHILDREN(n) (n)->nbits_children
#define BRANCH_GET_NUM_CHILDREN(t,n,l) (1ULL << (n)->nbits_children)

#define BRANCH_GET_CHILDREN(t,n) ((TrieNodeBranch *) (n))->children

typedef struct _TrieAllocator TrieAllocator;
struct _TrieAllocator {
  unsigned long long next;
  /* Span cache line */
  char dummy[120];
};

struct _Trie {
  unsigned int nbits;
  unsigned int nbits_root;
  TrieRef *roots;
  /* Allocators */
  pthread_mutex_t mutex;
  unsigned int nallocators;
  TrieAllocator *allocators;
  /* Blocks */
  unsigned long long nbranches;
  TrieNodeBranch *branches[1024];
  /* Debug */
  unsigned int num_allocations;
  unsigned long long total_memory;
};

Trie *trie_new (unsigned int nbits, unsigned int nbits_root, unsigned int nallocators);
void trie_setup (Trie *trie, unsigned int nbits, unsigned int nbits_root);
void trie_setup_full (Trie *trie, unsigned int nbits, unsigned int nbits_root, unsigned int nallocators);
void trie_release (Trie *trie);
void trie_add_word (Trie *trie, unsigned long long word, unsigned int count);
void trie_add_word_with_allocator (Trie *trie, unsigned long long word, unsigned int count, unsigned int aidx);
unsigned int trie_lookup (Trie *trie, unsigned long long word);

/* If callback return not 0 lookup stops */
int trie_foreach (Trie *trie, int (*callback) (Trie *trie, unsigned long long word, unsigned int code, unsigned long long word_idx, void *data), void *data);

unsigned int trie_setup_from_file (Trie *trie, FILE *ofs);
unsigned int trie_setup_from_data (Trie *trie, const unsigned char *cdata);
unsigned long long trie_write_to_file (Trie *trie, FILE *ofs);

#endif
