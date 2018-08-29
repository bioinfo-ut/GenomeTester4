#define __DATABASE_C__

#include <string.h>

#include "sequence.h"
#include "utils.h"

#include "database.h"

unsigned int db_debug = 0;
#define debug db_debug

void
gt4_db_clear_index (KMerDB *db)
{
  memset (&db->index, 0, sizeof (db->index));
}

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
      n_lines = 0;
      break;
    }
    if (!*wordsize) {
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
read_db_from_text (KMerDB *db, const unsigned char *cdata, unsigned long long csize, unsigned int max_kmers_per_node, unsigned int count_bits)
{
  unsigned long long cpos, last, div;
  unsigned int nlines, wordsize, n_kmers, max_kmers, names_size;
  unsigned int node_bits, kmer_bits;
  unsigned int names_pos, kmers_pos, idx;
  double t_s, t_e;

  if (csize < 8) return 0;
  if (!cdata[5] || !cdata[7]) return 0;

  div = csize / 100;
  
  wordsize = 0;  
  nlines = count_lines_from_text (cdata, csize, &wordsize, &n_kmers, &max_kmers, &names_size);
  if (debug) {
    fprintf (stderr, "Lines %u wordisze %u kmers %u (max per node %u) name size %u\n", nlines, wordsize, n_kmers, max_kmers, names_size);
  }
  if (!nlines) {
    fprintf (stderr, "File is not text-format kmer database (maybe binary?)\n");
    return 0;
  }
  if (max_kmers > max_kmers_per_node) {
    max_kmers = max_kmers_per_node;
  }
  node_bits = get_bits (nlines + 1);
  kmer_bits = get_bits (max_kmers);
  if ((node_bits + kmer_bits) > 31) {
    fprintf (stderr, "Too many nodes and kmers (%u (%u bits), %u (%u bits)\n", nlines + 1, max_kmers, node_bits, kmer_bits);
    return 0;
  }
  /* Set up DB */
  db->wordsize = wordsize;
  db->node_bits = node_bits;
  db->kmer_bits = kmer_bits;
  db->count_bits = count_bits;
  db->n_nodes = 0;
  db->nodes = (Node *) malloc (nlines * sizeof (Node));
  memset (db->nodes, 0, nlines * sizeof (Node));
  db->n_kmers = n_kmers;
  if (count_bits == 16) {
    db->kmers_16 = (unsigned short *) malloc (n_kmers * 2);
    memset (db->kmers_16, 0, n_kmers * 2);
  } else {
    db->kmers_32 = (unsigned int *) malloc (n_kmers * 4);
    memset (db->kmers_32, 0, n_kmers * 4);
  }
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
    const unsigned char *tokenz[65536];
    unsigned int lengths[65536];
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
      unsigned int dir = 0;

      while ((cdata[cpos] < ' ') && (cpos < csize)) cpos += 1;
      if ((csize - cpos) < db->wordsize) break;

      word = string_to_word ((const char *) cdata + cpos, db->wordsize);
      rword = get_reverse_complement (word, db->wordsize);
      if (rword < word) {
        word = rword;
        dir = 0x80000000;
      }

      /* Calculate code */
      code = dir | ((idx + 1) << db->kmer_bits) | i;

      if (debug) {
        code2 = trie_lookup (&db->trie, word);
        if (code2 != 0) {
          unsigned int idx2, kmer2;
          idx2 = ((code2 & 0x7fffffff) >> db->kmer_bits) - 1;
          kmer2 = code2 & ((1 << db->kmer_bits) - 1);
          fprintf (stderr, "KMer already present (current node %u (%s) kmer %u/%u (%s) code %u) previous %u (%s) kmer %u/%u code %u\n",
            idx, db->names + db->nodes[idx].name, i, (dir != 0), word_to_string (word, db->wordsize), code,
            idx2, db->names + db->nodes[idx2].name, kmer2, ((code2 & 0x7fffffff) != 0), code2);
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

  if (debug) {
    fprintf (stderr, "Database layout\n");
    fprintf (stderr, "  Wordsize: %u\n", db->wordsize);
    fprintf (stderr, "  Node bits: %u\n", db->node_bits);
    fprintf (stderr, "  KMer bits: %u\n", db->kmer_bits);
    fprintf (stderr, "  Count bits: %u\n", db->count_bits);
    fprintf (stderr, "  Nodes: %llu\n", db->n_nodes);
    fprintf (stderr, "  Kmers: %llu\n", db->n_kmers);
    fprintf (stderr, "  Names size: %llu\n", db->names_size);
  }
  
  return idx;
}

static const char *DBKEY = "GMDB";

unsigned int
write_db_to_file (KMerDB *db, FILE *ofs, unsigned int kmers)
{
  unsigned long long written = 0, blocksize, nodes_start, kmers_start, names_start, trie_start, index_start;
  static const unsigned short major = 0, minor = 3;
  fwrite (DBKEY, 4, 1, ofs);
  fwrite (&major, 2, 1, ofs);
  fwrite (&minor, 2, 1, ofs);

  fwrite (&db->wordsize, 4, 1, ofs);
  fwrite (&db->node_bits, 4, 1, ofs);
  fwrite (&db->kmer_bits, 4, 1, ofs);
  fwrite (&db->count_bits, 4, 1, ofs);
  fwrite (&db->n_nodes, 8, 1, ofs);
  fwrite (&db->n_kmers, 8, 1, ofs);
  fprintf (stderr, "Names size %llu\n", db->names_size);
  fwrite (&db->names_size, 8, 1, ofs);
  written = 48;
  /* Starts (5 * 8 bytes) */
  fseek (ofs, written + 40, SEEK_SET);
  written += 40;
  /* Nodes */
  nodes_start = written;
  blocksize = (db->n_nodes * sizeof (Node) + 15) & 0xfffffffffffffff0;
  fwrite (&blocksize, 8, 1, ofs);
  written += 8;
  fwrite (db->nodes, sizeof (Node), db->n_nodes, ofs);
  written += blocksize;
  fseek (ofs, written, SEEK_SET);
  /* KMers */
  kmers_start = written;
  if (kmers) {
    if (db->count_bits == 16) {
      blocksize = (db->n_kmers * 2 + 15) & 0xfffffffffffffff0;
      fwrite (&blocksize, 8, 1, ofs);
      written += 8;
      fwrite (db->kmers_16, 2, db->n_kmers, ofs);
      written += blocksize;
    } else if (db->count_bits == 32) {
      blocksize = (db->n_kmers * 4 + 15) & 0xfffffffffffffff0;
      fwrite (&blocksize, 8, 1, ofs);
      written += 8;
      fwrite (db->kmers_32, 4, db->n_kmers, ofs);
      written += blocksize;
    }
    fseek (ofs, written, SEEK_SET);
  } else {
    blocksize = 0;
    fwrite (&blocksize, 8, 1, ofs);
    written += 8;
    fseek (ofs, written, SEEK_SET);
  }
  /* Names */
  names_start = written;
  blocksize = (db->names_size + 15) & 0xfffffffffffffff0;
  fwrite (&blocksize, 8, 1, ofs);
  written += 8;
  fwrite (db->names, 1, db->names_size, ofs);
  written += blocksize;
  fseek (ofs, written, SEEK_SET);
  /* Trie */
  trie_start = written;
  fwrite (&blocksize, 8, 1, ofs);
  written += 8;
  blocksize = trie_write_to_file (&db->trie, ofs);
  blocksize = (blocksize + 15) & 0xfffffffffffffff0;
  written += blocksize;
  fseek (ofs, trie_start, SEEK_SET);
  fwrite (&blocksize, 8, 1, ofs);
  fseek (ofs, written, SEEK_SET);
  /* Index */
  index_start = written;
  fwrite (&blocksize, 8, 1, ofs);
  written += 8;
  blocksize = gt4_index_write (&db->index, ofs, db->n_kmers);
  blocksize = (blocksize + 15) & 0xfffffffffffffff0;
  written += blocksize;
  fseek (ofs, index_start, SEEK_SET);
  fwrite (&blocksize, 8, 1, ofs);
  fseek (ofs, written, SEEK_SET);
  /* Rewrite start locations */
  fseek (ofs, 48, SEEK_SET);
  fwrite (&nodes_start, 8, 1, ofs);
  fwrite (&kmers_start, 8, 1, ofs);
  fwrite (&names_start, 8, 1, ofs);
  fwrite (&trie_start, 8, 1, ofs);
  fwrite (&index_start, 8, 1, ofs);

  if (debug) {
    fprintf (stderr, "Database layout\n");
    fprintf (stderr, "  Wordsize: %u\n", db->wordsize);
    fprintf (stderr, "  Node bits: %u\n", db->node_bits);
    fprintf (stderr, "  KMer bits: %u\n", db->kmer_bits);
    fprintf (stderr, "  Count bits: %u\n", db->count_bits);
    fprintf (stderr, "  Nodes: %llu\n", db->n_nodes);
    fprintf (stderr, "  Kmers: %llu\n", db->n_kmers);
    fprintf (stderr, "  Names size: %llu\n", db->names_size);
    fprintf (stderr, "  Nodes start: %llu\n", nodes_start);
    fprintf (stderr, "  KMers start: %llu\n", kmers_start);
    fprintf (stderr, "  Names start: %llu\n", names_start);
    fprintf (stderr, "  Trie start: %llu\n", trie_start);
    fprintf (stderr, "  Index start: %llu\n", index_start);
  }
  
  return written;
}

unsigned int
read_database_from_binary (KMerDB *db, const unsigned char *cdata, unsigned long long csize)
{
  unsigned long long cpos = 0, blocksize, nodes_start, kmers_start, names_start, trie_start, index_start;
  unsigned short major, minor;
  unsigned int version;
  unsigned int has_index = 0;

  if (memcmp (cdata + cpos, DBKEY, 4)) return 0;
  cpos += 4;
  memcpy (&major, cdata + cpos, 2);
  cpos += 2;
  memcpy (&minor, cdata + cpos, 2);
  cpos += 2;
  if (debug) fprintf (stderr, "Database version %u.%u\n", major, minor);
  version = (major << 16) | minor;
  if (version < 1) return 0;
  if (version >= 3) has_index = 1;

  memcpy (&db->wordsize, cdata + cpos, 4);
  cpos += 4;
  memcpy (&db->node_bits, cdata + cpos, 4);
  cpos += 4;
  memcpy (&db->kmer_bits, cdata + cpos, 4);
  cpos += 4;
  if (version == 0) {
    db->count_bits = 16;
  } else {
    memcpy (&db->count_bits, cdata + cpos, 4);
  }
  cpos += 4;
  memcpy (&db->n_nodes, cdata + cpos, 8);
  cpos += 8;
  memcpy (&db->n_kmers, cdata + cpos, 8);
  cpos += 8;
  memcpy (&db->names_size, cdata + cpos, 8);
  cpos += 8;
  if (version > 1) {
    /* Start positions */
    memcpy (&nodes_start, cdata + cpos, 8);
    cpos += 8;
    memcpy (&kmers_start, cdata + cpos, 8);
    cpos += 8;
    memcpy (&names_start, cdata + cpos, 8);
    cpos += 8;
    memcpy (&trie_start, cdata + cpos, 8);
    cpos += 8;
    memcpy (&index_start, cdata + cpos, 8);
    cpos += 8;
  } else {
    nodes_start = kmers_start = names_start = trie_start = index_start = 0;
  }

  if (debug) {
    fprintf (stderr, "Database layout\n");
    fprintf (stderr, "  Wordsize: %u\n", db->wordsize);
    fprintf (stderr, "  Node bits: %u\n", db->node_bits);
    fprintf (stderr, "  KMer bits: %u\n", db->kmer_bits);
    fprintf (stderr, "  Count bits: %u\n", db->count_bits);
    fprintf (stderr, "  Nodes: %llu\n", db->n_nodes);
    fprintf (stderr, "  Kmers: %llu\n", db->n_kmers);
    fprintf (stderr, "  Names size: %llu\n", db->names_size);
    fprintf (stderr, "  Nodes start: %llu\n", nodes_start);
    fprintf (stderr, "  KMers start: %llu\n", kmers_start);
    fprintf (stderr, "  Names start: %llu\n", names_start);
    fprintf (stderr, "  Trie start: %llu\n", trie_start);
    fprintf (stderr, "  Index start: %llu\n", index_start);
  }

  /* Nodes */
  if (version > 1) cpos = nodes_start;
  memcpy (&blocksize, cdata + cpos, 8);
  cpos += 8;
  db->nodes = (Node *) (cdata + cpos);
  cpos += blocksize;
  /* Kmers */
  if (db->count_bits == 16) {
    db->kmers_16 = (unsigned short *) malloc (db->n_kmers * 2);
    if (version > 1) cpos = kmers_start;
    memcpy (&blocksize, cdata + cpos, 8);
    cpos += 8;
    if (blocksize >= db->n_kmers * 2) {
      memcpy (db->kmers_16, cdata + cpos, db->n_kmers * 2);
    } else {
      memset (db->kmers_16, 0, db->n_kmers * 2);
    }
  } else if (db->count_bits == 32) {
    db->kmers_32 = (unsigned int *) malloc (db->n_kmers * 4);
    if (version > 1) cpos = kmers_start;
    memcpy (&blocksize, cdata + cpos, 8);
    cpos += 8;
    if (blocksize >= db->n_kmers * 4) {
      memcpy (db->kmers_32, cdata + cpos, db->n_kmers * 4);
    } else {
      memset (db->kmers_32, 0, db->n_kmers * 4);
    }
  }
  cpos += blocksize;
  /* Names */
  if (version > 1) cpos = names_start;
  memcpy (&blocksize, cdata + cpos, 8);
  cpos += 8;
  db->names = (char *) (cdata + cpos);
  cpos += blocksize;
  /* Trie */
  if (version > 1) {
    cpos = trie_start;
    memcpy (&blocksize, cdata + cpos, 8);
    cpos += 8;
  }
  trie_setup_from_data (&db->trie, cdata + cpos);
  /* Index */
  if (has_index) {
    cpos = index_start;
    memcpy (&blocksize, cdata + cpos, 8);
    cpos += 8;
    gt4_index_init_from_data (&db->index, cdata + cpos, blocksize, db->n_kmers);
  }
  return 1;
}

ReadList
*gm4_read_list_new (void)
{
  static ReadList *lists = NULL;
  static unsigned int n_lists = 0;
  static unsigned int block_size = 0;
  if (n_lists >= block_size) {
    block_size = 8 * 1024 * 1024;
    lists = (ReadList *) malloc (block_size * sizeof (ReadList));
    n_lists = 0;
  }
  return &lists[n_lists++];
}
