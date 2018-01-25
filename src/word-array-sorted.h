#ifndef __GT4_WORD_ARRAY_SORTED_H__
#define __GT4_WORD_ARRAY_SORTED_H__

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

/*
 * An interface for accessing words and counts sequentially
 *
 * After initialization it should point to the first word
 */

typedef struct _GT4WordSArrayImplementation GT4WordSArrayImplementation;
typedef struct _GT4WordSArrayInstance GT4WordSArrayInstance;
typedef struct _GT4WordSArrayClass GT4WordSArrayClass;

#define GT4_TYPE_WORD_SARRAY (gt4_word_sarray_get_type ())

#include <az/interface.h>

struct _GT4WordSArrayImplementation {
  AZImplementation implementation;
  /* Invariant: (idx == 0) && (num_words > 0) */
  unsigned int (* get_first_word) (GT4WordSArrayImplementation *impl, GT4WordSArrayInstance *inst);
  /* Invariant: (idx == next) && (idx <= num_words) */
  unsigned int (* get_next_word) (GT4WordSArrayImplementation *impl, GT4WordSArrayInstance *inst);
};

struct _GT4WordSArrayInstance {
  unsigned long long num_words;
  unsigned long long idx;
  unsigned long long word;
  unsigned int count;
  unsigned int word_length;
};

struct _GT4WordSArrayClass {
  AZInterfaceClass interface_class;
};

unsigned int gt4_word_sarray_get_type (void);

/* All methods return 1 on success, 0 on error */

/* Enforces invariant */
unsigned int gt4_word_sarray_get_first_word (GT4WordSArrayImplementation *impl, GT4WordSArrayInstance *inst);
/* Enforces invariant */
unsigned int gt4_word_sarray_get_next_word (GT4WordSArrayImplementation *impl, GT4WordSArrayInstance *inst);

#endif
