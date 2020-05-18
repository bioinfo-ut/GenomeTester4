#define __GT4_WORD_TABLE_C__

/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 * 
 * Cpyright (C) 2014 University of Tartu
 *
 * Authors: Maarja Lepamets and Lauris Kaplinski
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "sequence.h"
#include "common.h"
#include "utils.h"
#include "version.h"
#include "word-list.h"

#include "word-table.h"

unsigned int debug_tables = 0;
unsigned long long total_memory = 0;

void
gt4_word_table_setup (GT4WordTable *table, unsigned int wordlen, unsigned long long n_slots, unsigned int data_size)
{
  memset (table, 0, sizeof (GT4WordTable));
  table->wordlength = wordlen;
  table->data_size = data_size;
  gt4_word_table_ensure_size (table, n_slots);
}

void
gt4_word_table_release (GT4WordTable *table)
{
  if (debug_tables) {
    unsigned long long size =  table->n_word_slots * 8 + table->n_word_slots * table->data_size;
    fprintf (stderr, "wordtable_delete: Releasing %llu total %.2fG\n", (unsigned long long) size, (double) total_memory / 1073741824.0);
    total_memory -= size;
  }
  free (table->words);
  if (table->data) free (table->data);
}

GT4WordTable *
gt4_word_table_new (unsigned int wordlength, unsigned long long n_slots, unsigned int data_size)
{
  GT4WordTable *table = (GT4WordTable *) malloc (sizeof (GT4WordTable));
  if (!table) return NULL;
  memset (table, 0, sizeof (GT4WordTable));
  table->wordlength = wordlength;
  table->data_size = data_size;
  if (gt4_word_table_ensure_size (table, n_slots)) {
    gt4_word_table_delete (table);
    return NULL;
  }
  return table;
}

void 
gt4_word_table_delete (GT4WordTable *table)
{
  if (debug_tables) {
    unsigned long long size =  table->n_word_slots * 8 + table->n_word_slots * table->data_size;
    fprintf (stderr, "wordtable_delete: Releasing %llu total %.2fG\n", (unsigned long long) size, (double) total_memory / 1073741824.0);
    total_memory -= size;
  }
  free (table->words);
  if (table->data) free (table->data);
  free (table);
}

void 
gt4_word_table_clear (GT4WordTable *table)
{
  table->n_words = 0L;
}

int 
gt4_word_table_ensure_size (GT4WordTable *table, unsigned long long size)
{
	if (table->n_word_slots < size) {
		table->n_word_slots = size;
		table->words = (unsigned long long *) realloc (table->words, table->n_word_slots * sizeof (unsigned long long));
		if (!table->words) return GT_OUT_OF_MEMORY_ERROR;
		if (table->data_size) {
			table->data = (unsigned char *) realloc (table->data, table->n_word_slots * table->data_size);
			if (!table->data) return GT_OUT_OF_MEMORY_ERROR;
		}
	}
	return 0;
}

#define WORDTABLE_MIN_SIZE 20000000

static int 
wordtable_enlarge (GT4WordTable *table)
{
	unsigned long long nslots;
	nslots = table->n_word_slots * 2;
	if (nslots < WORDTABLE_MIN_SIZE) nslots = WORDTABLE_MIN_SIZE;
	return gt4_word_table_ensure_size (table, nslots);
}

static int 
wordtable_enlarge_nofreq (GT4WordTable *table)
{
	return wordtable_enlarge (table);
}

int 
gt4_word_table_add_word (GT4WordTable *table, unsigned long long word, void *data)
{
	if (table->n_words >= table->n_word_slots) {
		int v;
		v = wordtable_enlarge (table);
		if (v > 0) return v;
	}
	table->words[table->n_words] = word;
	memcpy (table->data + table->n_words * table->data_size, data, table->data_size);
	table->n_words += 1;
	return 0;
}

int 
gt4_word_table_add_word_nofreq (GT4WordTable *table, unsigned long long word)
{
	if (table->n_words >= table->n_word_slots) {
		int v;
		v = wordtable_enlarge_nofreq (table);
		if (v > 0) return v;
	}
	table->words[table->n_words] = word;
	table->n_words += 1;
	return 0;
}

int 
wordtable_merge (GT4WordTable *table, GT4WordTable *other)
{
	long long i = 0L, j = 0L, k;
	unsigned long long nnew, incr, nequals;
	int v;
	unsigned int *t_freqs = (unsigned int *) table->data;
	unsigned int *o_freqs = (unsigned int *) other->data;

	if (table->wordlength != other->wordlength) return GT_INCOMPATIBLE_WORDLENGTH_ERROR;

	nequals = 0L;
	nnew = 0L;
	while ((unsigned long long) j < other->n_words && (unsigned long long) i < table->n_words) {

		if (table->words[i] == other->words[j]) {
			t_freqs[i] += o_freqs[j];
			i += 1;
			j += 1;
			nequals += 1;
		} else if (table->words[i] < other->words[j]) {
			i += 1;
		} else {
			j += 1;
		}
	}
	nnew = other->n_words - nequals;
	/* if (nnew == 0) return 0; */
	incr = nnew;
	if ((table->n_words + incr) > table->n_word_slots) {
		unsigned long long step = (table->n_word_slots + 7) >> 3;
		if (incr < step) incr = step;
	}

	v = gt4_word_table_ensure_size (table, table->n_words + incr);
	if (v > 0) return v;

	i = table->n_words - 1;
	j = other->n_words - 1;

	for (k = table->n_words + nnew - 1; k >= 0; k--) {

		if (j >= 0 && i >= 0 && table->words[i] == other->words[j]) {
			table->words[k] = table->words[i];
			t_freqs[k] = t_freqs[i];
			i -= 1;
			j -= 1;
		} else if ((j < 0) || (i >= 0 && table->words[i] > other->words[j])) {
			table->words[k] = table->words[i];
			t_freqs[k] = t_freqs[i];
			i -= 1;
		} else {
			table->words[k] = other->words[j];
			t_freqs[k] = o_freqs[j];
			j -= 1;
		}
	}

	table->n_words += nnew;
	return 0;
}

void 
wordtable_sort (GT4WordTable *table, int sort_data)
{
  unsigned int shift = 0;
  if (table->n_words == 0) return;
  /* Calculate the number of shifted positions for making radix sort faster (no need to sort digits that are all zeros)*/
  while (shift + 8 < table->wordlength * 2) {
    shift += 8;
  }
  if (sort_data) {
    hybridInPlaceRadixSort256 (table->words, table->words + table->n_words, table->data, table->data_size, shift);
  } else {
    hybridInPlaceRadixSort256 (table->words, table->words + table->n_words, NULL, 0, shift);
  }
}

int 
wordtable_find_frequencies (GT4WordTable *table)
{
  unsigned long long ri, wi, count;
  unsigned int *freqs = (unsigned int *) table->data;
  assert (table->data_size == 4);
  if (table->n_words == 0) return 0;

  wi = 0;
  count = 1;
  for (ri = 1; ri < table->n_words; ri++) {
    if (table->words[ri] == table->words[ri - 1]) {
      count += 1;
    } else {
      table->words[wi] = table->words[ri - 1];
      freqs[wi] = count;
      count = 1;
      wi += 1;
    }
  }

  table->words[wi] = table->words[ri - 1];
  freqs[wi] = count;
  table->n_words = wi + 1;
  return 0;
}

void 
wordtable_merge_freqs (GT4WordTable *table)
{
  unsigned long long ri, wi, count;
  unsigned int *freqs = (unsigned int *) table->data;

  wi = 0;
  if (table->n_words == 0) return;

  count = freqs[0];
  for (ri = 1; ri < table->n_words; ri++) {
    if (table->words[ri] == table->words[ri - 1]) {
      count += freqs[ri];
    } else {
      table->words[wi] = table->words[ri - 1];
      freqs[wi] = count;
      count = freqs[ri];
      wi += 1;
    }
  }
  table->words[wi] = table->words[ri - 1];
  freqs[wi] = count;
  table->n_words = wi + 1;
}

unsigned long long 
wordtable_count_unique (GT4WordTable *table)
{
	unsigned long long i, count;
	count = 0;
	for (i = 0; i < table->n_words; i++) {
		if (i == 0 || table->words[i] != table->words[i - 1]) {
			count += 1;
		}
	}
	return count;
}

#define BSIZE 10000

unsigned int
wordtable_write_to_file (GT4WordTable *table, const char *outputname, unsigned int cutoff)
{
	unsigned long long i, count, totalfreq;
	char fname[256]; /* the length of the output name is limited and checked in main(..) method */
	FILE *f;
	GT4ListHeader h;
	char b[BSIZE + 12];
	unsigned int bp;
	unsigned int *freqs = (unsigned int *) table->data;
	if (table->n_words == 0) return 0;

	memset (&h, 0, sizeof (GT4ListHeader));

	sprintf (fname, "%s_%d.list", outputname, table->wordlength);
	f = fopen (fname, "w");
	if (!f) {
		fprintf (stderr, "Cannot open output file %s\n", fname);
		return 1;
	}
#if 0
	b = malloc (1024 * 1024);
	setvbuf (f, b, _IOFBF, 1024 * 1024);
#endif
	fwrite (&h, sizeof (GT4ListHeader), 1, f);

	count = 0;
	totalfreq = 0;
	bp = 0;
	for (i = 0; i < table->n_words; i++) {
		if (freqs[i] >= cutoff) {
			memcpy (b + bp, &table->words[i], 8);
			bp += 8;
			memcpy (b + bp, &freqs[i], 4);
			bp += 4;
			if (bp >= BSIZE) {
				fwrite (b, 1, bp, f);
				bp = 0;
			}
#if 0
			fwrite (&table->words[i], sizeof (table->words[i]), 1, f);
			fwrite (&freqs[i], sizeof (freqs[i]), 1, f);
#endif
			count += 1;
			totalfreq += freqs[i];
		}
	}
	if (bp) fwrite (b, 1, bp, f);
	gt4_list_header_init (&h, table->wordlength);
	h.n_words = count;
	h.total_count = totalfreq;
	fseek (f, 0, SEEK_SET);
	fwrite (&h, sizeof (GT4ListHeader), 1, f);
	fclose (f);
#if 0
	free (b);
#endif
	return 0;
}

unsigned long long
gt4_word_table_generate_mismatches (GT4WordTable *tbl, unsigned long long word, void *data, unsigned int n_mm, unsigned int start, unsigned int canonical, unsigned int count_only, unsigned int equal_mm_only)
{
  unsigned long long mask = 0L, count = 0L, mismatch = 0L;
  unsigned int i;
  /* The current word */
  if (!count_only && (!equal_mm_only || !n_mm)) {
    if (canonical) {
      gt4_word_table_add_word (tbl, get_canonical_word (word, tbl->wordlength), data);
    } else {
      gt4_word_table_add_word (tbl, word, data);
    }
  }
  if (!n_mm) return 1;
  /* Generate mismatches */
  for (i = start; i < tbl->wordlength; i++) {
    for (mismatch = 1; mismatch < 4; mismatch++) {
      if (!count_only) mask = mismatch << (2 * i);
      count += gt4_word_table_generate_mismatches (tbl, word ^ mask, data, n_mm - 1, i + 1, canonical, count_only, equal_mm_only);
    }
  }
  return count;
}
