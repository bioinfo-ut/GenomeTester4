#ifndef __GT4_FILE_ARRAY_H__
#define __GT4_FILE_ARRAY_H__

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

typedef struct _GT4FileArrayImplementation GT4FileArrayImplementation;
typedef struct _GT4FileArrayInstance GT4FileArrayInstance;
typedef struct _GT4FileArrayClass GT4FileArrayClass;

#define GT4_TYPE_FILE_ARRAY (gt4_file_array_get_type ())

#include <az/interface.h>

struct _GT4FileArrayImplementation {
  AZImplementation implementation;
  
  unsigned int (* get_file) (GT4FileArrayImplementation *impl, GT4FileArrayInstance *inst, unsigned int idx);
  unsigned int (* get_sequence) (GT4FileArrayImplementation *impl, GT4FileArrayInstance *inst, unsigned long long idx);
};

struct _GT4FileArrayInstance {
  unsigned int num_files;
  /* File data */
  const unsigned char *file_name;
  unsigned long long file_size;
  unsigned long long n_sequences;
  /* Sequence_data */
  unsigned long long name_pos;
  unsigned int name_len;
  unsigned long long seq_pos;
  unsigned long long seq_len;
};

struct _GT4FileArrayClass {
  AZInterfaceClass interface_class;
};

unsigned int gt4_file_array_get_type (void);

/* All methods return 1 on success, 0 on error */
unsigned int gt4_file_array_get_file (GT4FileArrayImplementation *impl, GT4FileArrayInstance *inst, unsigned int idx);
unsigned int gt4_file_array_get_sequence (GT4FileArrayImplementation *impl, GT4FileArrayInstance *inst, unsigned long long idx);

#endif
