#define __DATABASE_C__

#include <string.h>

#include "sequence.h"
#include "utils.h"

#include "database.h"

extern unsigned int debug;

static unsigned int
count_lines_from_text (const unsigned char *cdata, size_t csize, unsigned int *wordsize, unsigned int *n_kmers, unsigned int *max_kmers, unsigned int *names_size)
{
  unsigned long long cpos, last, div;
  unsigned int n_lines, node_n_kmers;
  double t_s;
  
  div = csize / 100;

  /* Count lines */
  t_s = get_time ();
  if (debug) fprintf (stderr, "Counting lines ");

  *n_kmers = 0;
  *max_kmers = 0;
  *names_size = 0;

  n_lines = 0;
  cpos = 0;
  last = 0;
  while (cpos < csize) {
    const unsigned char *tokenz[4];
    unsigned int lengths[4];
    unsigned int ntokenz;

    ntokenz = split_line (cdata + cpos, csize - cpos, tokenz, lengths, 3);
    if (ntokenz < 2) {
      fprintf (stderr, "Line %u has <2 (%u) tokens\n", n_lines, ntokenz);
      break;
    }
    if (n_lines == 0) {
      *wordsize = lengths[2];
      if (debug) fprintf (stderr, "Wordsize %u\n", *wordsize);
    }
    node_n_kmers = strtol ((const char *) tokenz[1], NULL, 10);
    *n_kmers += node_n_kmers;
    if (node_n_kmers > *max_kmers) {
      *max_kmers = node_n_kmers;
    }
    *names_size += (lengths[0] + 1);

    while ((cdata[cpos] < csize) && (cdata[cpos] != '\n')) cpos++;
    if (cpos < csize) cpos += 1;
    n_lines += 1;
    if ((debug > 0) && ((cpos / div) > last)) {
      fprintf (stderr, ".");
      last += 1;
    }
  }
  double t_e = get_time ();
  if (debug) {
    fprintf (stderr, " %u\n", n_lines);
    fprintf (stderr, "Time %.1f (%.0f lines/s)\n", t_e - t_s, n_lines / (t_e - t_s));
  }
  return n_lines;
}

static unsigned int
get_bits (unsigned int value) {
  unsigned int bits = 0;
  while (value > 0) {
    bits += 1;
    value /= 2;
  }
  return bits;
}

unsigned int
read_db_from_text (KMerDB *db, const unsigned char *cdata, unsigned long long csize, unsigned int max_kmers_per_node)
{
  unsigned long long cpos, last, div;
  unsigned int nlines, wordsize, n_kmers, max_kmers, names_size;
  unsigned int names_pos, kmers_pos, idx;

  double t_s, t_e;

  div = csize / 100;
  
  wordsize = 0;  
  nlines = count_lines_from_text (cdata, csize, &wordsize, &n_kmers, &max_kmers, &names_size);
  if (debug) {
    fprintf (stderr, "Lines %u wordisze %u kmers %u (max per node %u) name size %u\n", nlines, wordsize, n_kmers, max_kmers, names_size);
  }
  if (max_kmers > max_kmers_per_node) {
    max_kmers = max_kmers_per_node;
  }
  /* Set up DB */
  db->wordsize = wordsize;
  db->node_bits = get_bits (nlines + 1);
  db->kmer_bits = get_bits (max_kmers);
  db->n_nodes = 0;
  db->nodes = (Node *) malloc (nlines * sizeof (Node));
  memset (db->nodes, 0, nlines * sizeof (Node));
  db->n_kmers = n_kmers;
  db->kmers = (unsigned short *) malloc (n_kmers * 2);
  memset (db->kmers, 0, n_kmers * 2);
  db->names_size = names_size;
  db->names = (char *) malloc (names_size);
  memset (db->names, 0, names_size);
  
  if (debug) {
    fprintf (stderr, "Node bits %u kmer bits %u\n", db->node_bits, db->kmer_bits);
  }

  trie_setup (&db->trie, db->wordsize * 2, 28);

  /* Fill table */
  t_s = get_time ();
  if (debug) fprintf (stderr, "Building database ");
  names_pos = 0;
  kmers_pos = 0;
  idx = 0;
  
  cpos = 0;
  last = 0;
  while (cpos < csize) {
    const unsigned char *tokenz[256];
    unsigned int lengths[256];
    unsigned int ntokenz, n_kmers;
    unsigned int i;
    /* Initialize */
    memset (db->nodes + idx, 0, sizeof (Node));
    /* Parse ID + number of kmers */
    ntokenz = split_line (cdata + cpos, csize - cpos, tokenz, lengths, 3);
    if (ntokenz < 2) {
      fprintf (stderr, "Line %u has <2 (%u) tokens\n", idx, ntokenz);
      break;
    }
    /* Name */
    db->nodes[idx].name = names_pos;
    memcpy (db->names + names_pos, tokenz[0], lengths[0]);
    if (debug > 2) fprintf (stderr, "Names %u %u %s\n", names_pos, lengths[0], db->names + names_pos);
    names_pos += lengths[0];
    names_pos += 1;
    /* Kmers */
    db->nodes[idx].kmers = kmers_pos;
    /* Number of kmers */
    n_kmers = strtol ((const char *) tokenz[1], NULL, 10);
    if (n_kmers > max_kmers_per_node) {
      n_kmers = max_kmers_per_node;
    }
    db->nodes[idx].nkmers = n_kmers;
    kmers_pos += n_kmers;
    /* Add kmers to trie */
    if (ntokenz > 2) cpos = tokenz[2] - cdata;
    for (i = 0; i < n_kmers; i++) {
      unsigned long long word, rword;
      unsigned int code, code2;

      while ((cdata[cpos] < ' ') && (cpos < csize)) cpos += 1;
      if ((csize - cpos) < db->wordsize) break;

      word = string_to_word ((const char *) cdata + cpos, db->wordsize);
      rword = get_reverse_complement (word, db->wordsize);
      if (rword < word) word = rword;

      /* Calculate code */
      code = ((idx + 1) << db->kmer_bits) | i;

      if (debug) {
        code2 = trie_lookup (&db->trie, word);
        if (code2 != 0) {
          unsigned int idx2, kmer2;
          idx2 = (code2 >> db->kmer_bits) - 1;
          kmer2 = code2 & ((1 << db->kmer_bits) - 1);
          fprintf (stderr, "KMer already present (current node %u (%s) kmer %u (%s) code %u) previous %u (%s) kmer %u code %u\n",
            idx, db->names + db->nodes[idx].name, i, word_to_string (word, db->wordsize), code,
            idx2, db->names + db->nodes[idx2].name, kmer2, code2);
          break;
        }
      }
      trie_add_word (&db->trie, word, code);
      if (debug) {
        code2 = trie_lookup (&db->trie, word);
        if (code2 != code) {
          fprintf (stderr, "Trie inconsistency (node %u (%s) kmer %u (%s) code %u) retrieved code %u\n",
          idx, db->names + db->nodes[idx].name, i, word_to_string (word, db->wordsize), code, code2);
          fprintf (stderr, "The program may crash in lookup phase\n");
          break;
        }
      }
      while ((cdata[cpos] >= ' ') && (cpos < csize)) cpos += 1;
    }
    if (i == db->nodes[idx].nkmers) {
      idx += 1;
    } else {
      fprintf (stderr, "Inconsisten number of kmers at node %u: %u (should be %u)\n", idx, i, db->nodes[idx].nkmers);
    }
    if ((debug > 0) && ((cpos / div) > last)) {
      fprintf (stderr, ".");
      last += 1;
    }
    /* Go to line end */
    while ((cdata[cpos] != '\n') && (cpos < csize)) {
      cpos += 1;
    }
    if (cpos < csize) cpos += 1;
    if (idx > nlines) break;
  }
  t_e = get_time ();
  if (debug) {
    fprintf (stderr, " %u\n", idx);
    fprintf (stderr, "Time %.1f (%.0f SNV/s)\n", t_e - t_s, idx / (t_e - t_s));
  }
  
  db->n_nodes = idx;
  db->n_kmers = kmers_pos;
  db->names_size = names_pos;

  return idx;
}

unsigned int
write_db_to_file (KMerDB *db, FILE *ofs)
{
  unsigned long long written = 0, blocksize;
  fwrite (&db->wordsize, 4, 1, ofs);
  fwrite (&db->node_bits, 4, 1, ofs);
  fwrite (&db->kmer_bits, 4, 1, ofs);
  fwrite (&db, 4, 1, ofs);
  fwrite (&db->n_nodes, 8, 1, ofs);
  fwrite (&db->n_kmers, 8, 1, ofs);
  fwrite (&db->names_size, 8, 1, ofs);
  written = 40;
  /* Nodes */
  blocksize = (db->n_nodes * sizeof (Node) + 15) & 0xfffffffffffffff0;
  fwrite (&blocksize, 8, 1, ofs);
  written += 8;
  fwrite (db->nodes, sizeof (Node), db->n_nodes, ofs);
  written += blocksize;
  fseek (ofs, written, SEEK_SET);
  /* KMers */
  blocksize = (db->n_kmers * 2 + 15) & 0xfffffffffffffff0;
  fwrite (&blocksize, 8, 1, ofs);
  written += 8;
  fwrite (db->kmers, 2, db->n_kmers, ofs);
  written += blocksize;
  fseek (ofs, written, SEEK_SET);
  /* Names */
  blocksize = (db->names_size + 15) & 0xfffffffffffffff0;
  fwrite (&blocksize, 8, 1, ofs);
  written += 8;
  fwrite (db->names, 1, db->names_size, ofs);
  written += blocksize;
  fseek (ofs, written, SEEK_SET);
  /* Trie */
  written += trie_write_to_file (&db->trie, ofs);
  return written;
}

unsigned int
read_database_from_binary (KMerDB *db, const unsigned char *cdata, unsigned long long csize)
{
  unsigned long long cpos = 0, blocksize;
  memcpy (&db->wordsize, cdata + cpos, 4);
  cpos += 4;
  memcpy (&db->node_bits, cdata + cpos, 4);
  cpos += 4;
  memcpy (&db->kmer_bits, cdata + cpos, 4);
  cpos += 4;
  cpos += 4;
  memcpy (&db->n_nodes, cdata + cpos, 8);
  cpos += 8;
  memcpy (&db->n_kmers, cdata + cpos, 8);
  cpos += 8;
  memcpy (&db->names_size, cdata + cpos, 8);
  cpos += 8;
  /* Nodes */
  memcpy (&blocksize, cdata + cpos, 8);
  cpos += 8;
  db->nodes = (Node *) (cdata + cpos);
  cpos += blocksize;
  /* Kmers */
  db->kmers = (unsigned short *) malloc (db->n_kmers * 2);
  memcpy (&blocksize, cdata + cpos, 8);
  cpos += 8;
  memcpy (db->kmers, cdata + cpos, db->n_kmers * 2);
  cpos += blocksize;
  /* Names */
  memcpy (&blocksize, cdata + cpos, 8);
  cpos += 8;
  db->names = (char *) (cdata + cpos);
  cpos += blocksize;
  trie_setup_from_data (&db->trie, cdata + cpos);
  return 1;
}
