#define __GDISTRIBUTION_C__

/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 * 
 * Copyright (C) 2014 University of Tartu
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

#define debug 1

#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>

#include "wordmap.h"

#define ELSIZE (sizeof (unsigned long long) + sizeof (unsigned int))
#define MAP_WORD(m,i) *((unsigned long long *) ((m)->wordlist + (i) * ELSIZE))
#define MAP_COUNT(m,i) *((unsigned int *) ((m)->wordlist + (i) * ELSIZE + sizeof (unsigned long long)))

typedef struct {
  float freq;
  unsigned int count;
} Freq;

static void get_distribution (wordmap *maps[2]);

static void
print_usage (FILE *ofs) {
  fprintf (ofs, "gdistribution LIST LIST2\n");
}

int
main (int argc, const char *argv[])
{
  const char *names[2];
  wordmap *maps[2];
  Freq *freqs;
  
  if (argc < 3) {
    print_usage (stderr);
    exit (1);
  }

  names[0] = argv[1];
  names[1] = argv[2];

  if (debug) fprintf (stderr, "%s %s\n", names[0], names[1]);
  
  maps[0] = wordmap_new (names[0], 1);
  maps[1] = wordmap_new (names[1], 1);
  
  get_distribution (maps);
  
  return 0;
}

static int compare (const void *lhs, const void *rhs) {
  if (*((float *) lhs) < *((float *) rhs)) return -1;
  if (*((float *) lhs) == *((float *) rhs)) return 0;
  return 1;
}

static void
get_distribution (wordmap *maps[2])
{
  unsigned long long size, i0, i1, fidx;
  unsigned int j, count;
  float *flist, current;
  
  size = maps[0]->header->nwords + maps[1]->header->nwords;
  if (debug) fprintf (stderr, "Total size %llu\n", size);
  flist = (float *) malloc (size * sizeof (float));

  i0 = 0;
  i1 = 0;
  fidx = 0;

  if (debug) fprintf (stderr, "Finding intersection\n");
  while ((i0 < maps[0]->header->nwords) && (i1 < maps[1]->header->nwords)) {
    if (MAP_WORD(maps[0], i0) == MAP_WORD(maps[1], i1)) {
      unsigned int c0, c1;
      float freq;
      c0 = MAP_COUNT(maps[0], i0);
      c1 = MAP_COUNT(maps[1], i1);
      /* freq = (float) c1 / c0; */
      freq = (float) c1;
      flist[fidx++] = freq;
      if (debug > 1) fprintf (stderr, "%llu %u %llu %u Freq %.2f\n", MAP_WORD(maps[0], i0), c0, MAP_WORD(maps[1], i1), c1, freq);
      i0 += 1;
      i1 += 1;
    } else if (MAP_WORD(maps[0], i0) < MAP_WORD(maps[1], i1)) {
      flist[fidx++] = 0;
      i0 += 1;
    } else {
      i1 += 1;
    }
  }
  if (debug) fprintf (stderr, "Size %llu\n", fidx);
  if (fidx == 0) {
    return;
  }

  /* Sort */
  if (debug) fprintf (stderr, "Sorting\n");
  qsort (flist, fidx, sizeof (float), compare);
  if (debug) fprintf (stderr, "Done\n");
  
  if (debug > 1) {
    for (j = 0; j < fidx; j++) {
      fprintf (stderr, "%.2f\n", flist[j]);
    }
  }

  j = 0;
  while(j < fidx) {
    count = 0;
    current = flist[j];
    while ((j < fidx) && (flist[j] == current)) {
      count += 1;
      j += 1;
    }
    fprintf (stdout, "%g\t%u\n", current, count);
  }
}

