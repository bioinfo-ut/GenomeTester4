/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 * 
 * Cpyright (C) 2014 University of Tartu
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

#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>

#include "utils.h"

const unsigned char *
gt4_mmap (const char *filename, unsigned long long *size)
{
  struct stat st;
  int status, handle;
  const unsigned char *data;

  status = stat (filename, &st);
  if (status < 0) {
    perror ("gt4_mmap (stat)");
    return NULL;
  }

  handle = open (filename, O_RDONLY);
  if (handle < 0) {
    perror ("gt4_mmap (open)");
    return NULL;
  }

  data = mmap (NULL, st.st_size, PROT_READ, MAP_PRIVATE, handle, 0);
  if (data == (const unsigned char *) -1) {
    perror ("gt4_mmap (mmap)");
    return NULL;
  } else {
    *size = st.st_size;
  }

  close (handle);
  return data;
}

void
gt4_munmap (const unsigned char *cdata, unsigned long long csize)
{
  munmap ((void *) cdata, csize);
}

static unsigned int finish = 0;

static void *
scout_map (void *arg)
{
  MapData *map = (MapData *) arg;
  unsigned long long i, val = 0;
  for (i = 0; i < map->csize; i += 1000) {
    val += map->cdata[i];
    if (finish) break;
  }
  free (map);
  return (void *) val;
}

void
scout_mmap (const unsigned char *cdata, unsigned long long csize)
{
  pthread_t thread;
  MapData *map = (MapData *) malloc (sizeof (MapData));
  map->cdata = cdata;
  map->csize = csize;
  pthread_create (&thread, NULL, scout_map, map);
}

void
delete_scouts () {
  finish = 1;
}

/* this implementation is based on:
  http://www.drdobbs.com/architecture-and-design/algorithm-improvement-through-performanc/221600153?pgno=1

  In-place MSD hybrid radix sort (with insertion sort) */

/* used for small buckets */
void insertionSort (unsigned long long *begin, unsigned long long *end, unsigned char *data, unsigned int data_size)
{
  unsigned long long *p, *q;
  unsigned long long temp;
  unsigned char tmp[32];

  for (p = begin + 1; p != end; p++) {
    for (q = p; q != begin && *q < *(q - 1); --q) {
      temp = *q;
      *q = *(q - 1);
      *(q - 1) = temp;
      if (data) {
        memcpy (tmp, data + (q - begin) * data_size, data_size);
        memcpy (data + (q - begin) * data_size, data + (q - begin - 1) * data_size, data_size);
        memcpy (data + (q - begin - 1) * data_size, tmp, data_size);
      }
    }
  }
}

void hybridInPlaceRadixSort256 (unsigned long long *begin, unsigned long long *end, unsigned char *data, unsigned int data_size, unsigned int shift)
{
  unsigned long long *p;
  unsigned long long temp, position, digit;
  unsigned char tmp[32];

  size_t bins[256];  /* for counts */
  size_t positions[256]; /* for starting positions */
  size_t binsize[256]; /* for counting the filled positions */
  unsigned long long i;

  if (end - begin <= 32) {
    insertionSort (begin, end, data, data_size);
    return;
  }

  memset (bins, 0, sizeof (bins));
  memset (binsize, 0, sizeof (binsize));

  /* calculating counts for every bin */
  for (p = begin; p != end; p++) {
    digit = (*p >> shift) & 255;
    bins[digit] += 1;
  }

  /* finding starting positions for each bin */
  positions[0] = 0;
  for (i = 1; i < 256; i++) {
    positions[i] = positions[i - 1] + bins[i - 1];
  }

  /* swapping words and locations */
  for (i = 0; i < end - begin; )  {
    p = begin + i;
    digit = (*p >> shift) & 255;

    if (bins[digit] <= binsize[digit]) {
      i += 1;
      continue;
    }
    position = positions[digit] + binsize[digit];
    binsize[digit] += 1;
    if (p == begin + position) {
      i += 1;
      continue;
    }
    temp = *p;
    *p = begin[position];
    begin[position] = temp;

    if (data) {
      memcpy (tmp, data + i * data_size, data_size);
      memcpy (data + i * data_size, data + position * data_size, data_size);
      memcpy (data + position * data_size, tmp, data_size);
    }

  }

  /* recursive step */
  if (shift > 0) {
    for (i = 0; i < 256; i++) {
      if (bins[i] > 1) {
        if (data) {
          hybridInPlaceRadixSort256 (begin + positions[i], begin + positions[i] + bins[i], data + positions[i] * data_size, data_size, shift - 8);
        } else {
          hybridInPlaceRadixSort256 (begin + positions[i], begin + positions[i] + bins[i], NULL, 0, shift - 8);
        }
      }
    }
  }
  return;
}

double
get_time (void)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  double dtval = (tv.tv_sec + tv.tv_usec / 1000000.0);
  return dtval;
}

unsigned long long
rand_long_long (unsigned long long min, unsigned long long max)
{
  unsigned long long delta = max - min + 1;
  return min + (unsigned long long) (delta * (rand () / (RAND_MAX + 1.0)));
}

unsigned int
split_line_chr (const unsigned char *cdata, unsigned long long csize, const unsigned char *tokenz[], unsigned int lengths[], unsigned int max_tokens, unsigned int chr)
{
  unsigned int i = 0;
  unsigned long long s = 0;
  while ((i < max_tokens) && (cdata[s] != '\n')) {
    unsigned long long e = s;
    tokenz[i] = cdata + s;
    while ((e < csize) && (cdata[e] != chr) && (cdata[e] != '\n')) e += 1;
    lengths[i] = e - s;
    i += 1;
    s = e;
    if (cdata[s] != '\n') s += 1;
  }
  return i;
}

unsigned int
split_line (const unsigned char *cdata, unsigned long long csize, const unsigned char *tokenz[], unsigned int lengths[], unsigned int max_tokens)
{
  unsigned int i = 0;
  unsigned long long s = 0;
  while ((i < max_tokens) && (cdata[s] != '\n')) {
    unsigned long long e = s;
    tokenz[i] = cdata + s;
    while ((e < csize) && (cdata[e] >= ' ')) e += 1;
    lengths[i] = e - s;
    i += 1;
    s = e;
    if (cdata[s] != '\n') s += 1;
  }
  return i;
}

unsigned int
number_to_binary (char buf[], unsigned long long number, unsigned int ndigits)
{
  unsigned int pos;
  if (!ndigits == 0) {
    if (number) {
      unsigned long long tmp = number;
      while (tmp) {
        ndigits += 1;
        tmp = tmp >> 1;
      }
    } else {
      ndigits = 1;
    }
  } else if (ndigits > 64) {
    ndigits = 64;
  }
  buf[ndigits + 1] = 0;
  pos = ndigits;
  while (pos) {
    buf[pos] = '0' + (char) (number & 1);
    number = number >> 1;
    pos -= 1;
  }
  return ndigits;
}
