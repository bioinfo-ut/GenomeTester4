#ifndef __GT4_WORD_LIST_H__
#define __GT4_WORD_LIST_H__

/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 * 
 * Copyright (C) 2014-2016 University of Tartu
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

#include <stdint.h>

/*
 * GenomeTester4 list file format
 */

#ifndef __GT4_WORD_LIST_C__
/* List tag "GT4C" encoded to big-endian 32-bit integer */
extern unsigned int GT4_LIST_CODE;
#endif

typedef struct _GT4ListHeader_4_4 GT4ListHeader;

struct _GT4ListHeader_4_0 {
  uint32_t code;
  uint32_t version_major;
  uint32_t version_minor;
  uint32_t wordlength;
  uint64_t nwords;
  uint64_t totalfreq;
  uint64_t padding;
};
                                                        
struct _GT4ListHeader_4_2 {
  uint32_t code;
  uint32_t version_major;
  uint32_t version_minor;
  uint32_t word_length;
  uint64_t n_words;
  uint64_t total_count;
  /* Offset (from the beginning of header) of word/count data */
  uint64_t list_start;
};

struct _GT4ListHeader_4_4 {
  uint32_t code;
  uint32_t version_major;
  uint32_t version_minor;
  uint32_t word_length;
  uint64_t n_words;
  uint64_t total_count;
  /* Offset (from the beginning of header) of word/count data */
  uint64_t list_start;
  uint32_t word_bytes;
  uint32_t count_bytes;
};

void gt4_list_header_init (GT4ListHeader *hdr, unsigned int word_length);

#endif
