#ifndef __GT4_WORD_TABLE_H__
#define __GT4_WORD_TABLE_H__

/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 * 
 * Copyright (C) 2014-2018 University of Tartu
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

typedef struct _GT4WordTable GT4WordTable;

#include <stdio.h>

#ifndef __GT4_WORD_TABLE_C__
extern unsigned int debug_tables;
#endif

struct _GT4WordTable {
	unsigned int wordlength;
	unsigned long long nwordslots;
	unsigned long long nfreqslots;
	unsigned long long nwords;
	unsigned long long *words;
	unsigned int *freqs;
	/* Text ID for debugging */
        char _id[16];
};

GT4WordTable *wordtable_new (unsigned int wordlength, unsigned long long size);

void wordtable_delete (GT4WordTable *table);

void wordtable_empty (GT4WordTable *table);

int wordtable_enlarge (GT4WordTable *table);
int wordtable_enlarge_nofreq (GT4WordTable *table);

int wordtable_ensure_size (GT4WordTable *table, unsigned long long size, unsigned long long freqsize);

int wordtable_add_word (GT4WordTable *table, unsigned long long word, unsigned int freq, unsigned int wordlength);
int wordtable_add_word_nofreq (GT4WordTable *table, unsigned long long word, unsigned int wordlength);

int wordtable_merge (GT4WordTable *table, GT4WordTable *other);

void wordtable_sort (GT4WordTable *table, int sortfreqs);

int wordtable_find_frequencies (GT4WordTable *table);

void wordtable_merge_frequencies (GT4WordTable *table);

unsigned long long wordtable_count_unique(GT4WordTable *table);

unsigned int wordtable_write_to_file (GT4WordTable *table, const char *outputname, unsigned int cutoff);

void write_word_to_file (unsigned long long word, unsigned freq, FILE *f);

unsigned int wordtable_build_filename (GT4WordTable *table, char *c, unsigned int length, const char *prefix);

#endif
