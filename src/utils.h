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

/* memory-map a given file */
const char * mmap_by_filename (const char *filename, size_t *size);

/* memory-unmap a given file */
int munmap_by_file (const char *file, size_t *size);

void insertionSort (unsigned long long *begin, unsigned long long *end, unsigned int *begfreq);

void hybridInPlaceRadixSort256 (unsigned long long *begin, unsigned long long *end, unsigned int *begfreq, unsigned int shift);

double get_time (void);

#endif /* UTILS_H_ */
