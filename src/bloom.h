#ifndef __GT4_BLOOM_H__
#define __GT4_BLOOM_H__

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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
 
 /*
 * GT4WordMap is the most basic list container
 */

typedef struct _GT4Bloom GT4Bloom;

#include <stdint.h>

struct _GT4Bloom {
  uint64_t mask;
  uint32_t n_hashes;
  uint32_t _filler;
  uint64_t bits[1];
};

GT4Bloom *gt4_bloom_new (unsigned int substr_bits, unsigned int n_hashes);
void gt4_bloom_delete (GT4Bloom *bloom);

unsigned int gt4_bloom_insert (GT4Bloom *bloom, uint64_t word);
unsigned int gt4_bloom_query (GT4Bloom *bloom, uint64_t word);

#endif
