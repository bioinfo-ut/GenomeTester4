#define __GT4_WORD_LIST_C__

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

#include <string.h>

#include "version.h"
#include "word-list.h"

unsigned int GT4_LIST_CODE = 'G' << 24 | 'T' << 16 | '4' << 8 | 'C';

void
gt4_list_header_init (GT4ListHeader *hdr, unsigned int word_length)
{
  memset (hdr, 0, sizeof (GT4ListHeader));
  hdr->code = GT4_LIST_CODE;
  hdr->version_major = VERSION_MAJOR;
  hdr->version_minor = VERSION_MINOR;
  hdr->word_length = word_length;
  hdr->list_start = sizeof (GT4ListHeader);
  hdr->word_bytes = 8;
  hdr->count_bytes = 4;
}
