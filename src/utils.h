#ifndef UTILS_H_
#define UTILS_H_

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

#include <stdlib.h>

/* Memory-map a given file */
const unsigned char *gt4_mmap (const char *filename, unsigned long long *csize);

/* MMap scouting */

typedef struct _MapData MapData;
struct _MapData {
  const unsigned char *cdata;
  unsigned long long csize;
};

/* Create new thread that reads mmap into resident memory */
void scout_mmap (const unsigned char *cdata, unsigned long long csize);
void delete_scouts ();

/* Memory-unmap a previously mapped file */
void gt4_munmap (const unsigned char *cdata, unsigned long long csize);

void insertionSort (unsigned long long *begin, unsigned long long *end, unsigned char *data, unsigned int data_size);
void hybridInPlaceRadixSort256 (unsigned long long *begin, unsigned long long *end, unsigned char *data, unsigned int data_size, unsigned int shift);

double get_time (void);
unsigned long long rand_long_long (unsigned long long min, unsigned long long max);

/* Split line into tokens */
/* Line ends with \n (or at csize), tokens end with \t */
/* Return number of tokens */
unsigned int split_line (const unsigned char *cdata, unsigned long long csize, const unsigned char *tokenz[], unsigned int lengths[], unsigned int max_tokens);
unsigned int split_line_chr (const unsigned char *cdata, unsigned long long csize, const unsigned char *tokenz[], unsigned int lengths[], unsigned int max_tokens, unsigned int chr);

/* Print number as binary with given number of digits (0 - start from leftmost 1) */
/* Returns the length of string (not counting terminating 0) */
unsigned int number_to_binary (char buf[], unsigned long long number, unsigned int ndigits);

#endif /* UTILS_H_ */
