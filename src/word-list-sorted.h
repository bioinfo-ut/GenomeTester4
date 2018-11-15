#ifndef __GT4_WORD_LIST_SORTED_H__
#define __GT4_WORD_LIST_SORTED_H__

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

typedef struct _GT4WordSListImplementation GT4WordSListImplementation;
typedef struct _GT4WordSListInstance GT4WordSListInstance;
typedef struct _GT4WordSListClass GT4WordSListClass;

#define GT4_TYPE_WORD_SLIST (gt4_word_slist_get_type ())

#include <stdint.h>
#include <az/interface.h>

struct _GT4WordSListImplementation {
  AZImplementation implementation;
  /* Invariant: (idx == 0) && (num_words > 0) */
  unsigned int (* get_first_word) (GT4WordSListImplementation *impl, GT4WordSListInstance *inst);
  /* Invariant: (idx == next) && (idx <= num_words) */
  unsigned int (* get_next_word) (GT4WordSListImplementation *impl, GT4WordSListInstance *inst);
};

struct _GT4WordSListInstance {
  unsigned long long num_words;
  unsigned long long sum_counts;
  unsigned long long idx;
  unsigned long long word;
  unsigned int count;
  unsigned int word_length;
};

struct _GT4WordSListClass {
  AZInterfaceClass interface_class;
};

unsigned int gt4_word_slist_get_type (void);

/* All methods return 1 on success, 0 on error */

/* Enforces invariant */
unsigned int gt4_word_slist_get_first_word (GT4WordSListImplementation *impl, GT4WordSListInstance *inst);
/* Enforces invariant */
unsigned int gt4_word_slist_get_next_word (GT4WordSListImplementation *impl, GT4WordSListInstance *inst);

#endif
