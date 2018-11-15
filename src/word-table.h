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
typedef struct _GT4WordTableLocation GT4WordTableLocation;

#include <stdint.h>
#include <stdio.h>

#ifndef __GT4_WORD_TABLE_C__
extern unsigned int debug_tables;
#endif

struct _GT4WordTableLocation {
  unsigned int file_idx;
  unsigned int seq_idx;
  unsigned int seq_pos;
};

struct _GT4WordTable {
  unsigned long long n_word_slots;
  unsigned long long n_words;
  unsigned long long *words;
  unsigned char *data;
  unsigned int data_size;
  unsigned int wordlength;
};

GT4WordTable *gt4_word_table_new (unsigned int wordlen, unsigned long long n_slots, unsigned int data_size);
void gt4_word_table_delete (GT4WordTable *table);

void gt4_word_table_clear (GT4WordTable *table);
int gt4_word_table_ensure_size (GT4WordTable *table, unsigned long long size);

int gt4_word_table_add_word (GT4WordTable *table, unsigned long long word, void *data);
int gt4_word_table_add_word_nofreq (GT4WordTable *table, unsigned long long word);

int wordtable_merge (GT4WordTable *table, GT4WordTable *other);

void wordtable_sort (GT4WordTable *table, int sortfreqs);

int wordtable_find_frequencies (GT4WordTable *table);

void wordtable_merge_frequencies (GT4WordTable *table);

unsigned long long wordtable_count_unique(GT4WordTable *table);

unsigned int wordtable_write_to_file (GT4WordTable *table, const char *outputname, unsigned int cutoff);

#endif
