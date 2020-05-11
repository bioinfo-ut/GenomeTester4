#ifndef __GT4_WORD_INDEX_H__
#define __GT4_WORD_INDEX_H__

/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 * 
 * Copyright (C) 2014-2020 University of Tartu
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
 * An interface for accessing location data for words
 * Word lookup has to be implemented by separate interface
 */

typedef struct _GT4WordIndexImplementation GT4WordIndexImplementation;
typedef struct _GT4WordIndexInstance GT4WordIndexInstance;
typedef struct _GT4WordIndexClass GT4WordIndexClass;

#define GT4_TYPE_WORD_INDEX (gt4_word_index_get_type ())

#include <stdint.h>
#include <az/interface.h>

struct _GT4WordIndexImplementation {
  AZImplementation implementation;
  unsigned int (*get_location) (GT4WordIndexImplementation *impl, GT4WordIndexInstance *inst, unsigned long long idx);
};

struct _GT4WordIndexInstance {
  uint32_t n_locations;
  uint32_t file_idx;
  uint32_t seq_idx;
  uint32_t dir;
  uint64_t pos;
};

struct _GT4WordIndexClass {
  AZInterfaceClass interface_class;
};

unsigned int gt4_word_index_get_type (void);

/* All methods return 1 on success, 0 on error */
unsigned int gt4_word_index_get_location (GT4WordIndexImplementation *impl, GT4WordIndexInstance *inst, unsigned long long idx);

#endif
