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

#include <string.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>

#include "utils.h"

const char * mmap_by_filename (const char *filename, size_t *size)
{
	struct stat st;
	int status, handle;
	const char *data;

	status = stat (filename, &st);
	if (status < 0) {
		return NULL;
	}

	handle = open (filename, O_RDONLY);
	if (handle < 0) {
		return NULL;
	}

	data = (const char *) mmap (NULL, st.st_size, PROT_READ, MAP_PRIVATE, handle, 0);
	if (data == (const char *) -1) {
		return NULL;
	} else {
		*size = (size_t)st.st_size;
	}

	close (handle);
	return data;
}

int munmap_by_file (const char *file, size_t *size)
{
	munmap ((void *) file, *size);
	return 0;
}

/* this implementation is based on:
	http://www.drdobbs.com/architecture-and-design/algorithm-improvement-through-performanc/221600153?pgno=1

	In-place MSD hybrid radix sort (with insertion sort) */

/* used for small buckets */
void insertionSort (unsigned long long *begin, unsigned long long *end, unsigned int *begfreq)
{
	unsigned long long *p, *q;
	unsigned long long temp;
	unsigned int temp_freq;

	for (p = begin + 1; p != end; p++) {
		for (q = p; q != begin && *q < *(q - 1); --q) {
			temp = *q;
			*q = *(q - 1);
			*(q - 1) = temp;
			if (begfreq) {
				temp_freq = begfreq[q - begin];
				begfreq[q - begin] = begfreq[q - begin - 1];
				begfreq[q - begin - 1] = temp_freq;
			}
		}
	}
}

void hybridInPlaceRadixSort256 (unsigned long long *begin, unsigned long long *end, unsigned int *begfreq, unsigned int shift)
{
	unsigned long long *p;
	unsigned long long temp, position, digit;
	unsigned int temp_freq;

	size_t bins[256];	/* for counts */
	size_t positions[256]; /* for starting positions */
	size_t binsize[256]; /* for counting the filled positions */
	int i;

	if (end - begin <= 32) {
		insertionSort (begin, end, begfreq);
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

		if (begfreq) {
			temp_freq = begfreq[i];
			begfreq[i] = begfreq[position];
			begfreq[position] = temp_freq;
		}

	}

	/* recursive step */
	if (shift > 0) {
		for (i = 0; i < 256; i++) {
			if (bins[i] > 1) {
				if (begfreq) {
					hybridInPlaceRadixSort256 (begin + positions[i], begin + positions[i] + bins[i], begfreq + positions[i], shift - 8);
				} else {
					hybridInPlaceRadixSort256 (begin + positions[i], begin + positions[i] + bins[i], NULL, shift - 8);
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
