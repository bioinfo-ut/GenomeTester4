#define WORDTABLE_C

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

#include <stdlib.h>
#include <string.h>

#include "sequence.h"
#include "wordtable.h"
#include "common.h"
#include "utils.h"
#include "version.h"

unsigned int debug_tables;
unsigned long long total_memory = 0;

wordtable *
wordtable_new (unsigned int wordlength, unsigned long long size)
{
	static int idx = 0;
	int v;
	wordtable *table;
	table = (wordtable *) malloc (sizeof (wordtable));
	if (!table) return NULL;
	memset (table, 0, sizeof (wordtable));
	/* Set ID */
	char *c = "AA";
	strcpy (table->id, c);
	table->id[0] += idx / 26;
	table->id[1] += idx % 26;
	idx += 1;
	table->wordlength = wordlength;
	v = wordtable_ensure_size (table, size, 0);
	if (v != 0) return NULL;
	return table;
}

void 
wordtable_delete (wordtable *table)
{
	if (debug_tables) {
		unsigned long long size =  table->nwordslots * 8 + table->nfreqslots * 4;
		fprintf (stderr, "Table %s releasing %llu total %.2fG\n", table->id, (unsigned long long) size, (double) total_memory / 1073741824.0);
		total_memory -= size;
	}
	free (table->words);
	free (table->frequencies);
	free (table);
	return;
}

void 
wordtable_empty (wordtable *table)
{
	table->wordlength = 0;
	table->nwords = 0L;
	return;
}

int 
wordtable_ensure_size (wordtable *table, unsigned long long size, unsigned long long freqsize)
{
	if (table->nwordslots < size) {
		if (debug_tables) {
			unsigned long long asize =  (size - table->nwordslots) * 8;
			total_memory += asize;
			fprintf (stderr, "Table %s allocating words %llu total %.2fG\n", table->id, (unsigned long long) asize, (double) total_memory / 1073741824.0);
		}
		table->nwordslots = size;
		table->words = (unsigned long long *) realloc (table->words, table->nwordslots * sizeof (unsigned long long));
		if (!table->words) return GT_OUT_OF_MEMORY_ERROR;
	}
	if (table->nfreqslots < freqsize) {
		if (debug_tables) {
			unsigned long long asize =  (freqsize - table->nfreqslots) * 4;
			total_memory += asize;
			fprintf (stderr, "Table %s allocating freqs %llu total %.2fG\n", table->id, (unsigned long long) asize, (double) total_memory / 1073741824.0);
		}
		table->nfreqslots = freqsize;
		table->frequencies = (unsigned int *) realloc (table->frequencies, table->nfreqslots * sizeof (unsigned int));
		if (!table->frequencies) return GT_OUT_OF_MEMORY_ERROR;
	}
	return 0;
}

#define WORDTABLE_MIN_SIZE 20000000

int 
wordtable_enlarge (wordtable *table)
{
	unsigned long long nslots;
	int v;
	if (table->nwordslots < WORDTABLE_MIN_SIZE && table->nfreqslots < WORDTABLE_MIN_SIZE) {
		nslots = WORDTABLE_MIN_SIZE;
	} else {
		nslots = (table->nwordslots > table->nfreqslots ? table->nwordslots : table->nfreqslots) * 2;
	}
	v = wordtable_ensure_size (table, nslots, nslots);
	if (v) return v;
	return 0;
}

int 
wordtable_enlarge_nofreq (wordtable *table)
{
	unsigned long long nslots;
	int v;
	if (table->nwordslots < 10000000 && table->nfreqslots < 10000000) {
		nslots = 10000000;
	} else {
		nslots = (table->nwordslots > table->nfreqslots ? table->nwordslots : table->nfreqslots) * 2;
	}
	v = wordtable_ensure_size (table, nslots, 0);
	if (v) return v;
	return 0;
}

int 
wordtable_add_word (wordtable *table, unsigned long long word, unsigned int freq, unsigned int wordlength)
{
	if (table->wordlength != wordlength) {
		return GT_INCOMPATIBLE_WORDLENGTH_ERROR;
	}
	if (table->nwords == table->nfreqslots || table->nwords == table->nwordslots) {
		int v;
		v = wordtable_enlarge (table);
		if (v > 0) return v;
	}
	table->words[table->nwords] = word;
	table->frequencies[table->nwords] = freq;
	table->nwords += 1;
	return 0;
}

int 
wordtable_add_word_nofreq (wordtable *table, unsigned long long word, unsigned int wordlength)
{
#if 0
	if (table->wordlength != wordlength) {
		return GT_INCOMPATIBLE_WORDLENGTH_ERROR;
	}
#endif
	if (table->nwords >= table->nwordslots) {
		int v;
		v = wordtable_enlarge_nofreq (table);
		if (v > 0) return v;
	}
	table->words[table->nwords] = word;
	table->nwords += 1;
	return 0;
}

int 
wordtable_merge (wordtable *table, wordtable *other)
{
	long long i = 0L, j = 0L, k;
	unsigned long long nnew, incr, nequals;
	int v;

	if (table->wordlength != other->wordlength) return GT_INCOMPATIBLE_WORDLENGTH_ERROR;

	nequals = 0L;
	nnew = 0L;
	while ((unsigned long long) j < other->nwords && (unsigned long long) i < table->nwords) {

		if (table->words[i] == other->words[j]) {
			table->frequencies[i] += other->frequencies[j];
			i += 1;
			j += 1;
			nequals += 1;
		} else if (table->words[i] < other->words[j]) {
			i += 1;
		} else {
			j += 1;
		}
	}
	nnew = other->nwords - nequals;
	/* if (nnew == 0) return 0; */
	incr = nnew;
	if ((table->nwords + incr) > table->nwordslots) {
		unsigned long long step = (table->nwordslots + 7) >> 3;
		if (incr < step) incr = step;
	}

	v = wordtable_ensure_size (table, table->nwords + incr, table->nwords + incr);
	if (v > 0) return v;

	i = table->nwords - 1;
	j = other->nwords - 1;

	for (k = table->nwords + nnew - 1; k >= 0; k--) {

		if (j >= 0 && i >= 0 && table->words[i] == other->words[j]) {
			table->words[k] = table->words[i];
			table->frequencies[k] = table->frequencies[i];
			i -= 1;
			j -= 1;
		} else if ((j < 0) || (i >= 0 && table->words[i] > other->words[j])) {
			table->words[k] = table->words[i];
			table->frequencies[k] = table->frequencies[i];
			i -= 1;
		} else {
			table->words[k] = other->words[j];
			table->frequencies[k] = other->frequencies[j];
			j -= 1;
		}
	}

	table->nwords += nnew;
	return 0;
}

void 
wordtable_sort (wordtable *table, int sortfreqs)
{
	unsigned int firstshift = 0;
	if (table->nwords == 0) return;

	/* calculate the number of shifted positions for making radix sort faster (no need to sort digits that are all zeros)*/
	while (firstshift + 8 < table->wordlength * 2) {
		firstshift += 8;
	}

	if (sortfreqs) {
		hybridInPlaceRadixSort256 (table->words, table->words + table->nwords, table->frequencies, firstshift);
		return;
	}
	hybridInPlaceRadixSort256 (table->words, table->words + table->nwords, NULL, firstshift);
	return;
}

int 
wordtable_find_frequencies (wordtable *table)
{
	unsigned long long ri, wi, count, nunique;
	int v;

	wi = 0;
	count = 1;
	if (table->nwords == 0) return 0;
	nunique = wordtable_count_unique (table);
	v = wordtable_ensure_size (table, nunique, nunique);
	if (v > 0) return v;

	for (ri = 1; ri < table->nwords; ri++) {
		if (table->words[ri] == table->words[ri - 1]) {
			count += 1;
		} else {
			table->words[wi] = table->words[ri - 1];
			table->frequencies[wi] = count;
			count = 1;
			wi += 1;
		}
	}

	table->words[wi] = table->words[ri - 1];
	table->frequencies[wi] = count;
	table->nwords = nunique;
	return 0;
}

void 
wordtable_merge_frequencies (wordtable *table)
{
	unsigned long long ri, wi, count, nunique;

	wi = 0;
	if (table->nwords == 0) return;
	nunique = wordtable_count_unique (table);

	count = table->frequencies[0];
	for (ri = 1; ri < table->nwords; ri++) {
		if (table->words[ri] == table->words[ri - 1]) {
			count += table->frequencies[ri];
		} else {
			table->words[wi] = table->words[ri - 1];
			table->frequencies[wi] = count;
			count = table->frequencies[ri];
			wi += 1;
		}
	}

	table->words[wi] = table->words[ri - 1];
	table->frequencies[wi] = count;
	table->nwords = nunique;
	return;
}

unsigned long long 
wordtable_count_unique (wordtable *table)
{
	unsigned long long i, count;
	count = 0;
	for (i = 0; i < table->nwords; i++) {
		if (i == 0 || table->words[i] != table->words[i - 1]) {
			count += 1;
		}
	}
	return count;
}

#define BSIZE 10000

unsigned int
wordtable_write_to_file (wordtable *table, const char *outputname, unsigned int cutoff)
{
	unsigned long long i, count, totalfreq;
	char fname[256]; /* the length of the output name is limited and checked in main(..) method */
	FILE *f;
	GT4ListHeader h;
	char b[BSIZE + 12];
	unsigned int bp;
	if (table->nwords == 0) return 0;

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
	for (i = 0; i < table->nwords; i++) {
		if (table->frequencies[i] >= cutoff) {
			memcpy (b + bp, &table->words[i], 8);
			bp += 8;
			memcpy (b + bp, &table->frequencies[i], 4);
			bp += 4;
			if (bp >= BSIZE) {
				fwrite (b, 1, bp, f);
				bp = 0;
			}
#if 0
			fwrite (&table->words[i], sizeof (table->words[i]), 1, f);
			fwrite (&table->frequencies[i], sizeof (table->frequencies[i]), 1, f);
#endif
			count += 1;
			totalfreq += table->frequencies[i];
		}
	}
	if (bp) {
		fwrite (b, 1, bp, f);
	}

	h.code = GT4_LIST_CODE;
	h.version_major = VERSION_MAJOR;
	h.version_minor = VERSION_MINOR;
	h.wordlength = table->wordlength;
	h.nwords = count;
	h.totalfreq = totalfreq;
	h.padding = sizeof (GT4ListHeader);
	fseek (f, 0, SEEK_SET);
	fwrite (&h, sizeof (GT4ListHeader), 1, f);
	fclose (f);
#if 0
	free (b);
#endif
	return 0;
}

void write_word_to_file (unsigned long long word, unsigned freq, FILE *f)
{
	fwrite (&word, sizeof (unsigned long long), 1, f);
	fwrite (&freq, sizeof (unsigned int), 1, f);
	return;
}

unsigned long long generate_mismatches (wordtable *mmtable, unsigned long long word, unsigned int wordlength,
		unsigned int givenfreq, unsigned int nmm, unsigned int startsite, int usesmallercomplement, int countonly,
		int equalmmonly)
{
	unsigned long long mask = 0L, count = 0L, mismatch = 0L;
	unsigned int i;

	/* first I put the current word into the table */
	if (!countonly && (nmm == 0 || !equalmmonly)) {
		if (usesmallercomplement) {
			word = get_canonical_word (word, wordlength);
		}
		wordtable_add_word (mmtable, word, givenfreq, wordlength);
	}
	if (nmm == 0) return 1;

	/* generating mm-s */
	for (i = startsite; i < wordlength; i++) {
		for (mismatch = 1; mismatch < 4; mismatch++) {
			if (!countonly) {
				mask = mismatch << (2 * i);
			}
			count += generate_mismatches (mmtable, word ^ mask, wordlength, givenfreq, nmm - 1, i + 1, usesmallercomplement, countonly, equalmmonly);
		}
	}
	return count;
}

/* Generates filename from wordtable and prefix */
/* Returns the length of full filename */
/* If c is NULL only length is returned */
/* If length is insufficient string is truncated but full length returned */

unsigned int
wordtable_build_filename (wordtable *table, char *c, unsigned int length, const char *prefix)
{
	unsigned int p, q;
	char suffix[32];
	sprintf (suffix, "_%d.list", table->wordlength);
	if (c) {
		p = 0;
		q = 0;
		while (prefix[p] && (q < (length - 1))) c[q++] = prefix[p++];
		p = 0;
		while (suffix[p] && (q < (length - 1))) c[q++] = suffix[p++];
		c[q] = 0;
	}
	return strlen (prefix) + strlen (suffix);
}
