#define __GT4_BLOOM_C__

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

#include <malloc.h>
#include <string.h>

#include "bloom.h"

GT4Bloom *
gt4_bloom_new (unsigned int substr_bits, unsigned int n_hashes)
{
  GT4Bloom *bloom;
  uint64_t n_bits = 1, size;
  while (substr_bits) {
    n_bits = (n_bits << 1);
    substr_bits -= 1;
  }
  size = ((n_bits + 63) / 64) * 8;
  bloom = (GT4Bloom *) malloc (sizeof (GT4Bloom) + size - 8);
  memset (bloom, 0, sizeof (GT4Bloom) + size - 8);
  bloom->mask = n_bits - 1;
  bloom->n_hashes = n_hashes;
  return bloom;
}

void
gt4_bloom_delete (GT4Bloom *bloom)
{
  free (bloom);
}

unsigned int
gt4_bloom_insert (GT4Bloom *bloom, uint64_t word)
{
  unsigned int i, result = 1;
  for (i = 0; i < bloom->n_hashes; i++) {
    uint64_t p = (word & bloom->mask) / 64;
    uint64_t m = (1ULL << (word & 63));
    result &= ((bloom->bits[p] & m) != 0);
    bloom->bits[p] |= m;
    word = (word >> 1);
  }
  return result;
}

unsigned int
gt4_bloom_query (GT4Bloom *bloom, uint64_t word)
{
  unsigned int i;
  for (i = 0; i < bloom->n_hashes; i++) {
    uint64_t p = (word & bloom->mask) / 64;
    uint64_t m = (1ULL << (word & 63));
    if (!(bloom->bits[p] & m)) return 0;
    word = (word >> 1);
  }
  return 1;
}

