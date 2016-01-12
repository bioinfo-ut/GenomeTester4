#include <malloc.h>
#include <string.h>
#include <assert.h>
#include <strings.h>

#include "trie.h"

#define debug 1

static TrieRef trie_node_add_word (Trie *trie, TrieRef ref, unsigned int level, unsigned long long word, unsigned int nbits, unsigned int count, unsigned int aidx);
static unsigned int trie_node_lookup (Trie *trie, TrieRef ref, unsigned long long word, unsigned int nbits);

#define ALLOCATOR_BLOCK_SIZE 65536

Trie *
trie_new (unsigned int nbits, unsigned int nbits_root, unsigned int nallocators)
{
  Trie *trie = (Trie *) malloc (sizeof (Trie));
  trie_setup_full (trie, nbits, nbits_root, nallocators);
  return trie;
}

void
trie_setup_full (Trie *trie, unsigned int nbits, unsigned int nbits_root, unsigned int nallocators)
{
  memset (trie, 0, sizeof (Trie));
  trie->nbits = nbits;
  trie->nbits_root = (nbits > nbits_root) ? nbits_root : nbits;
  trie->roots = (TrieRef *) malloc ((1ULL << trie->nbits_root) * sizeof (TrieRef));
  memset (trie->roots, 0, (1ULL << trie->nbits_root) * sizeof (TrieRef));

  pthread_mutex_init (&trie->mutex, NULL);
  trie->nallocators = nallocators;
  trie->allocators = (TrieAllocator *) memalign (64, nallocators * sizeof (TrieAllocator));
  memset (trie->allocators, 0, nallocators * sizeof (TrieAllocator));
}

void
trie_setup (Trie *trie, unsigned int nbits, unsigned int nbits_root)
{
  trie_setup_full (trie, nbits, nbits_root, 1);
}

void
trie_release (Trie *trie)
{
  free (trie->allocators);
  pthread_mutex_destroy (&trie->mutex);
  
  /* fixme: Release nodes (Lauris) */
  free (trie->roots);
}

void
trie_add_word (Trie *trie, unsigned long long word, unsigned int count)
{
  trie_add_word_with_allocator (trie, word, count, 0);
}

void
trie_add_word_with_allocator (Trie *trie, unsigned long long word, unsigned int count, unsigned int aidx)
{
  unsigned int cbits = trie->nbits - trie->nbits_root;
  trie->roots[word >> cbits] = trie_node_add_word (trie, trie->roots[word >> cbits], 0, word % (1ULL << cbits), cbits, count, aidx);
}

unsigned int
trie_lookup (Trie *trie, unsigned long long word)
{
  unsigned int cbits = trie->nbits - trie->nbits_root;
  return trie_node_lookup (trie, trie->roots[word >> cbits],  word % (1ULL << cbits), cbits);
}

unsigned int
trie_setup_from_file (Trie *trie, FILE *ifs)
{
  unsigned long long i, nbranches;
  /* Trie */
  memset (trie, 0, sizeof (Trie));
  fread (&trie->nbits, 4, 1, ifs);
  fread (&trie->nbits_root, 4, 1, ifs);
  fread (&trie->nbranches, 8, 1, ifs);
  /* Roots */
  unsigned long long nroots = 1ULL << trie->nbits_root;
  trie->roots = (TrieRef *) malloc (nroots * sizeof (TrieRef));
  for (i = 0; i < nroots; i++) {
    fread (&trie->roots[i], sizeof (TrieRef), 1, ifs);
  }
  /* Blocks */
  nbranches = trie->nbranches;
  while (nbranches > 0) {
    unsigned long long size = nbranches;
    if (size > TRIE_BLOCK_SIZE) size = TRIE_BLOCK_SIZE;
    trie->branches[i] = (TrieNodeBranch *) malloc (size * sizeof (TrieNodeBranch));
    fread (trie->branches[i], sizeof (TrieNodeBranch), size, ifs);
    nbranches -= size;
  }
  return 0;
}

unsigned int
trie_setup_from_data (Trie *trie, const unsigned char *cdata)
{
  const unsigned char *p;
  int i;
  unsigned long long nbranches;
  p = cdata;
  /* Trie */
  memset (trie, 0, sizeof (Trie));
  memcpy (&trie->nbits, p, 4);
  p += 4;
  memcpy (&trie->nbits_root, p, 4);
  p += 4;
  memcpy (&trie->nbranches, p, 8);
  p += 8;
  /* Roots */
  unsigned long long nroots = 1ULL << trie->nbits_root;
  trie->roots = (TrieRef *) p;
  p += nroots * sizeof (TrieRef);
  /* Blocks */
  nbranches = trie->nbranches;
  i = 0;
  while (nbranches > 0) {
    unsigned long long size = nbranches;
    if (size > TRIE_BLOCK_SIZE) size = TRIE_BLOCK_SIZE;
    trie->branches[i] = (TrieNodeBranch *) p;
    p += size * sizeof (TrieNodeBranch);
    nbranches -= size;
    i += 1;
  }
  return 0;
}

unsigned long long
trie_write_to_file (Trie *trie, FILE *ofs)
{
  unsigned long long len = 0, i;
  /* Trie */
  fwrite (&trie->nbits, 4, 1, ofs);
  len += 4;
  fwrite (&trie->nbits_root, 4, 1, ofs);
  len += 4;
  fwrite (&trie->nbranches, 8, 1, ofs);
  len += 8;
  /* Roots */
  unsigned long long nroots = 1ULL << trie->nbits_root;
  for (i = 0; i < nroots; i++) {
    fwrite (&trie->roots[i], sizeof (TrieRef), 1, ofs);
    len += sizeof (TrieRef);
  }
  /* Blocks */
  for (i = 0; i < 1024; i++) {
    if (trie->branches[i] != NULL) {
      unsigned long long size = trie->nbranches - i * TRIE_BLOCK_SIZE;
      if (size > TRIE_BLOCK_SIZE) size = TRIE_BLOCK_SIZE;
      fwrite (trie->branches[i], sizeof (TrieNodeBranch), size, ofs);
      len += size * sizeof (TrieNodeBranch);
    }
  }
  return len;
}

static TrieRef
trie_allocate_branch (Trie *trie, unsigned int aidx)
{
  unsigned long long block, idx;
  /* Try allocator */
  if ((trie->allocators[aidx].next & (ALLOCATOR_BLOCK_SIZE - 1)) == 0) {
    /* Grab new block */
    pthread_mutex_lock (&trie->mutex);
    block = trie->nbranches / TRIE_BLOCK_SIZE;
    idx = trie->nbranches % TRIE_BLOCK_SIZE;
    idx = ((idx + ALLOCATOR_BLOCK_SIZE - 1) / ALLOCATOR_BLOCK_SIZE) * ALLOCATOR_BLOCK_SIZE;
    if ((idx + ALLOCATOR_BLOCK_SIZE) > TRIE_BLOCK_SIZE) {
      block += 1;
      idx = 0;
    }
    if (!trie->branches[block]) {
      /* fprintf (stderr, "trie_allocate_branch: new block %llu\n", block); */
      trie->branches[block] = (TrieNodeBranch *) malloc (TRIE_BLOCK_SIZE * sizeof (TrieNodeBranch));
    }
    trie->allocators[aidx].next = block * TRIE_BLOCK_SIZE + idx;
    trie->nbranches += ALLOCATOR_BLOCK_SIZE;
    pthread_mutex_unlock (&trie->mutex);
  }
  block = trie->allocators[aidx].next / TRIE_BLOCK_SIZE;
  idx = trie->allocators[aidx].next % TRIE_BLOCK_SIZE;
  trie->allocators[aidx].next += 1;
  memset (trie->branches[block] + idx, 0, sizeof (TrieNodeBranch));
  return TRIE_REF_FROM_ADDRESS(trie, block, idx)
}

static TrieRef
trie_node_kmer_new (Trie *trie, unsigned long long word, unsigned int nbits, unsigned int count)
{
  return (TrieRef) MAKE_KMER (nbits, word, count);
}

static TrieRef
trie_node_branch_new (Trie *trie, unsigned long long word, unsigned int nbits_this, unsigned int nbits_children, unsigned int aidx)
{
  TrieRef ref;
  TrieNodeBranch *branch;

  assert (nbits_this <= BRANCH_MAX_BITS_THIS);
  /* assert (nbits_children <= BRANCH_MAX_BITS_CHILDREN); */
  assert (nbits_children == 1);
  /* ref = trie_allocate_branch (trie, 1ULL << nbits_children, aidx); */
  ref = trie_allocate_branch (trie, aidx);
  branch = BRANCH_FROM_REF (trie, ref);
  branch->_nbits_this = nbits_this;
  branch->nbits_children = nbits_children;
  branch->word = word;
  return ref;
}

static TrieRef
trie_node_kmer_add_word (Trie *trie, TrieRef ref, unsigned int level, unsigned long long word, unsigned int nbits, unsigned int count, unsigned int aidx)
{
  unsigned long long kmer = ref;
  assert (REF_IS_KMER (ref));
  assert (nbits <= KMER_MAX_BITS);
  assert (nbits == KMER_GET_NBITS (ref));
  if (KMER_GET_WORD (ref) == word) {
    /* Same kmer, increase count */
    /* kmer->count = kmer->count + count; */
    kmer = MAKE_KMER (KMER_GET_NBITS (kmer), KMER_GET_WORD (kmer), KMER_GET_COUNT (kmer) + count);
    
    return (TrieRef) kmer;
  } else {
    unsigned int bit, new_this_bits, new_child_bits, child_kmer_bits;
    unsigned int old_idx;
    TrieRef new_ref;
    TrieNodeBranch *new_branch;
    TrieRef new_node_ref;

    /* Different kmer, split */
    bit = 63 - __builtin_clzll (KMER_GET_WORD (kmer) ^ word);

    new_this_bits = KMER_GET_NBITS (kmer) - bit - 1;
    new_child_bits = 1;
    child_kmer_bits = bit;

    old_idx = (KMER_GET_WORD (kmer) >> bit) & 1;

    /* fixme: Add children directly */
    new_ref = trie_node_branch_new (trie, word >> (bit + 1), new_this_bits, new_child_bits, aidx);
    new_branch = BRANCH_FROM_REF (trie, new_ref);

    /* new_node = trie_node_add_word (trie, (TrieNode *) new_branch, level + 1, kmer->word, nbits, kmer->count); */
    /* Update this */
    kmer = MAKE_KMER (child_kmer_bits, KMER_GET_WORD (kmer) % (1ULL << child_kmer_bits), KMER_GET_COUNT (kmer));
    /* Assign to new node */
    new_branch->children[old_idx] = (TrieRef) kmer;
    
    new_node_ref = trie_node_add_word (trie, new_ref, level + 1, word, nbits, count, aidx);

    return new_node_ref;
  }
}

static TrieRef
trie_node_branch_split (Trie *trie, TrieRef ref, unsigned int level, unsigned int bit, unsigned int aidx)
{
    unsigned int new_this_bits, new_child_bits, child_kmer_bits;
    unsigned int old_idx;
    TrieRef new_ref;
    TrieNodeBranch *new_branch;
    TrieNodeBranch *branch = BRANCH_FROM_REF(trie,ref);

    /* Different kmer, split */
    new_this_bits = BRANCH_GET_NBITS_THIS (branch) - bit - 1;
    new_child_bits = 1;
    child_kmer_bits = bit;

    old_idx = (branch->word >> bit) & 1;

    /* Create new parent */
    new_ref = trie_node_branch_new (trie, branch->word >> (bit + 1), new_this_bits, new_child_bits, aidx);
    new_branch = BRANCH_FROM_REF(trie,new_ref);
    /* Update this */
    branch->_nbits_this = child_kmer_bits;
    branch->word = branch->word % (1ULL << child_kmer_bits);
    /* Assign to new node */
    new_branch->children[old_idx] = ref;

    return new_ref;
}

static TrieRef
trie_node_branch_add_word (Trie *trie, TrieRef ref, unsigned int level, unsigned long long word, unsigned int nbits, unsigned int count, unsigned int aidx)
{
  unsigned long long lword;
  TrieNodeBranch *branch = BRANCH_FROM_REF(trie, ref);
  assert (branch != NULL);
  assert (REF_IS_BRANCH (ref));
  lword = word >> (nbits - BRANCH_GET_NBITS_THIS (branch));
  if (branch->word == lword) {
    /* Same kmer */
    unsigned long long nbits_this, nbits_children, cword, dword;
    nbits_this = BRANCH_GET_NBITS_THIS (branch);
    nbits_children = BRANCH_GET_NBITS_CHILDREN (branch);
    cword = (word >> (nbits - nbits_this - nbits_children)) % (1ULL << nbits_children);
    dword = word % (1ULL << (nbits - nbits_this - nbits_children));
    branch->children[cword] = trie_node_add_word (trie, branch->children[cword], level + 1, dword, nbits - nbits_this - nbits_children, count, aidx);
    return ref;
  } else {
    unsigned int bit;
    bit = 63 - __builtin_clzll (branch->word ^ lword);
    ref = trie_node_branch_split (trie, ref, level, bit, aidx);
    return trie_node_branch_add_word (trie, ref, level + 1, word, nbits, count, aidx);
  }
}

static TrieRef
trie_node_add_word (Trie *trie, TrieRef ref, unsigned int level, unsigned long long word, unsigned int nbits, unsigned int count, unsigned int aidx)
{
  if (REF_IS_EMPTY (ref)) {
    assert (nbits <= 64);
    if (nbits <= KMER_MAX_BITS) {
      return trie_node_kmer_new (trie, word, nbits, count);
    } else {
      unsigned int nbits_remaining = nbits;
      /* kmer */
      nbits_remaining -= KMER_MAX_BITS;
      /* children */
      nbits_remaining -= 1;
      /* clip */
      if (nbits_remaining > BRANCH_MAX_BITS_THIS) nbits_remaining = BRANCH_MAX_BITS_THIS;
      TrieRef branch = trie_node_branch_new (trie, word >> (nbits - nbits_remaining), nbits_remaining, 1, aidx);
      return trie_node_branch_add_word (trie, branch, level + 1, word, nbits, count, aidx);
    }
  }
  if (REF_IS_KMER (ref)) {
    return trie_node_kmer_add_word (trie, ref, level, word, nbits, count, aidx);
  } else {
    return trie_node_branch_add_word (trie, ref, level, word, nbits, count, aidx);
  }
}

static unsigned int
trie_node_kmer_lookup (TrieRef kmer, unsigned long long word, unsigned int nbits)
{
  assert (!REF_IS_EMPTY(kmer));
  assert (REF_IS_KMER (kmer));
  assert (nbits <= KMER_MAX_BITS);
  assert (nbits == KMER_GET_NBITS (kmer));
  if (KMER_GET_WORD (kmer) == word) {
    return KMER_GET_COUNT (kmer);
  } else {
    return 0;
  }
}

static unsigned int
trie_node_branch_lookup (Trie *trie, TrieRef ref, unsigned long long word, unsigned int nbits)
{
  unsigned long long lword;
  TrieNodeBranch *branch = BRANCH_FROM_REF(trie, ref);
  assert (REF_IS_BRANCH(ref));
  assert (branch != NULL);
  lword = word >> (nbits - BRANCH_GET_NBITS_THIS (branch));
  if (branch->word == lword) {
    /* Same kmer */
    unsigned long long nbits_this, nbits_children, cword, dword;
    nbits_this = BRANCH_GET_NBITS_THIS (branch);
    nbits_children = BRANCH_GET_NBITS_CHILDREN (branch);
    cword = (word >> (nbits - nbits_this - nbits_children)) % (1ULL << nbits_children);
    dword = word % (1ULL << (nbits - nbits_this - nbits_children));
    return trie_node_lookup (trie, branch->children[cword], dword, nbits - nbits_this - nbits_children);
  } else {
    return 0;
  }
}

static unsigned int
trie_node_lookup (Trie *trie, TrieRef ref, unsigned long long word, unsigned int nbits)
{
  if (REF_IS_EMPTY(ref)) return 0;
  if (REF_IS_KMER (ref)) {
    return trie_node_kmer_lookup (ref, word, nbits);
  } else {
    return trie_node_branch_lookup (trie, ref, word, nbits);
  }
}
