#ifndef WORDTABLE_H_
#define WORDTABLE_H_

/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 * 
 * Copyright (C) 2014 University of Tartu
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

#include <stdio.h>
#include "wordmap.h"

#ifndef WORDTABLE_C
extern unsigned int debug_tables;
#endif

typedef struct _wordtable {
	unsigned int wordlength;	/* initial length of the considered words */
        char id[16];	/* Text ID for debugging */
	unsigned long long nwordslots;	/* number of slots */
	unsigned long long nfreqslots;
	unsigned long long nwords;		/* filled slots */
	unsigned long long *words;
	unsigned int *frequencies;
} wordtable;

wordtable *wordtable_new (unsigned int wordlength, unsigned long long size);

void wordtable_delete (wordtable *table);

void wordtable_empty (wordtable *table);

int wordtable_enlarge (wordtable *table);
int wordtable_enlarge_nofreq (wordtable *table);

int wordtable_ensure_size (wordtable *table, unsigned long long size, unsigned long long freqsize);

int wordtable_add_word (wordtable *table, unsigned long long word, unsigned int freq, unsigned int wordlength);
int wordtable_add_word_nofreq (wordtable *table, unsigned long long word, unsigned int wordlength);

int wordtable_merge (wordtable *table, wordtable *other);

void wordtable_sort (wordtable *table, int sortfreqs);

int wordtable_find_frequencies (wordtable *table);

void wordtable_merge_frequencies (wordtable *table);

unsigned long long wordtable_count_unique(wordtable *table);

unsigned int wordtable_write_to_file (wordtable *table, const char *outputname, unsigned int cutoff);

void write_word_to_file (unsigned long long word, unsigned freq, FILE *f);

unsigned int wordtable_build_filename (wordtable *table, char *c, unsigned int length, const char *prefix);

#endif /* WORDTABLE_H_ */
