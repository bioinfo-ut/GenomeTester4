#define __GT4_SET_OPERATIONS_C__

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

extern int debug;

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <libarikkei/arikkei-utils.h>

#include "utils.h"
#include "version.h"

#include "set-operations.h"

#define TMP_BUF_SIZE (12 * 1024)

unsigned int
gt4_write_union (AZObject *arrays[], unsigned int n_arrays, unsigned int cutoff, int ofile, GT4ListHeader *header)
{
  GT4WordSListImplementation *impls[GT4_MAX_SETS];
  GT4WordSListInstance *insts[GT4_MAX_SETS];
  unsigned int n_sources;
  unsigned long long total = 0;
  unsigned int j;

  arikkei_return_val_if_fail (n_arrays <= GT4_MAX_SETS, 1);
  
  n_sources = 0;
  for (j = 0; j < n_arrays; j++) {
    impls[n_sources] = (GT4WordSListImplementation *) az_object_get_interface (AZ_OBJECT(arrays[j]), GT4_TYPE_WORD_SLIST, (void **) &insts[n_sources]);
    if (insts[n_sources]->num_words) {
      gt4_word_slist_get_first_word (impls[n_sources], insts[n_sources]);
      total += insts[n_sources]->num_words;
      n_sources += 1;
    }
  }

  gt4_list_header_init (header, insts[0]->word_length);

  if (ofile) write (ofile, header, sizeof (GT4ListHeader));

  if (n_sources) {
    unsigned long long word;
    unsigned char b[TMP_BUF_SIZE];
    unsigned int bp = 0;
    double t_s, t_e;
    /* Get timestamp */
    t_s = get_time ();
    /* Find first word */
    word = insts[0]->word;
    for (j = 1; j < n_sources; j++) if (insts[j]->word < word) word = insts[j]->word;
    /* Iterate until all lists are exhausted */
    while (n_sources) {
      unsigned long long next = 0xffffffffffffffff;
      unsigned int freq = 0;
      j = 0;
      while (j < n_sources) {
        if (insts[j]->word == word) {
          freq += insts[j]->count;
          if (!gt4_word_slist_get_next_word (impls[j], insts[j])) {
            n_sources -= 1;
            if (n_sources > 0) {
              impls[j] = impls[n_sources];
              insts[j] = insts[n_sources];
              continue;
            } else {
              break;
            }
          }
        }
        if (insts[j]->word < next) next = insts[j]->word;
        j += 1;
      }
      /* Now we have word and freq */
      if (freq >= cutoff) {
        if (ofile) {
          memcpy (&b[bp], &word, 8);
          memcpy (&b[bp + 8], &freq, 4);
          bp += 12;
          if (bp >= TMP_BUF_SIZE) {
            write (ofile, b, bp);
            bp = 0;
          }
        }
        header->n_words += 1;
        header->total_count += freq;
        if (debug && !(header->n_words % 100000000)) {
          fprintf (stderr, "Words written: %uM\n", (unsigned int) (header->n_words / 1000000));
        }
      }
      word = next;
    }
    if (ofile) {
      if (bp) write (ofile, b, bp);
      pwrite (ofile, header, sizeof (GT4ListHeader), 0);
    }
    t_e = get_time ();

    if (debug > 0) {
      fprintf (stderr, "Combined %u arrays: input %llu (%.3f Mwords/s) output %llu (%.3f Mwords/s)\n", n_arrays, total, total / (1000000 * (t_e - t_s)), header->n_words, header->n_words / (1000000 * (t_e - t_s)));
    }
  }
  
  return 0;
}
