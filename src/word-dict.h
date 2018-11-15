#ifndef __GT4_WORD_DICT_H__
#define __GT4_WORD_DICT_H__

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

typedef struct _GT4WordDictImplementation GT4WordDictImplementation;
typedef struct _GT4WordDictInstance GT4WordDictInstance;
typedef struct _GT4WordDictClass GT4WordDictClass;

#define GT4_TYPE_WORD_DICT (gt4_word_dict_get_type ())

#include <stdint.h>
#include <az/interface.h>

struct _GT4WordDictImplementation {
  AZImplementation implementation;
  unsigned int (* lookup) (GT4WordDictImplementation *impl, GT4WordDictInstance *inst, unsigned long long word);
};

struct _GT4WordDictInstance {
  unsigned long long num_words;
  unsigned int value;
  unsigned int word_length;
};

struct _GT4WordDictClass {
  AZInterfaceClass interface_class;
};

unsigned int gt4_word_dict_get_type (void);

/* All methods return 1 on success, 0 on error */
unsigned int gt4_word_dict_lookup (GT4WordDictImplementation *impl, GT4WordDictInstance *inst, unsigned long long word);

#endif
