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

/*
 * GenomeTester4 list file format
 */

#ifndef __GT4_WORD_LIST_C__
/* List tag "GT4C" encoded to big-endian 32-bit integer */
extern unsigned int GT4_LIST_CODE;
#endif

typedef struct _GT4ListHeader GT4ListHeader;

struct _GT4ListHeader {
  unsigned int code;
  unsigned int version_major;
  unsigned int version_minor;
  unsigned int word_length;
  unsigned long long n_words;
  unsigned long long total_count;
  /* Offset (from the beginning of header) of word/count data */
  unsigned long long list_start;
};

void gt4_list_header_init (GT4ListHeader *hdr, unsigned int word_length);

#endif
