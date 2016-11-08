#define __WORDMAP_C__

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>

#include "wordmap.h"
#include "wordtable.h"
#include "sequence.h"
#include "common.h"
#include "utils.h"
#include "queue.h"

unsigned int debug_wordmap = 0;

unsigned int GT4_LIST_CODE = 'G' << 24 | 'T' << 16 | '4' << 8 | 'C';

GT4WordMap * 
gt4_wordmap_new (const char *listfilename, unsigned int scout)
{
	const unsigned char *cdata;
	unsigned long long csize;
	GT4WordMap *map = (GT4WordMap *) malloc (sizeof (GT4WordMap));
	if (!map) {
		fprintf (stderr, "gt4_wordmap_new: could not allocate map\n");
		return NULL;
	}
	memset (map, 0, sizeof (GT4WordMap));

	map->filename = strdup (listfilename);
	cdata = gt4_mmap (listfilename, &csize);
	if (!cdata) {
		fprintf (stderr, "gt4_wordmap_new: could not mmap file %s\n", listfilename);
		gt4_wordmap_delete (map);
		return NULL;
	}
	map->file_map = cdata;
	map->file_size = csize;
	map->header = (GT4ListHeader *) cdata;
	if (map->header->code != GT4_LIST_CODE) {
		fprintf (stderr, "gt4_wordmap_new: invalid file tag (%x, should be %x)\n", map->header->code, GT4_LIST_CODE);
		gt4_munmap (cdata, csize);
		gt4_wordmap_delete (map);
		return NULL;
	}
	if (csize != sizeof (GT4ListHeader) + map->header->nwords * 12) {
		fprintf (stderr, "gt4_wordmap_new: invalid file size (%llu, should be %llu)\n", csize, sizeof (GT4ListHeader) + map->header->nwords * 12);
		gt4_munmap (cdata, csize);
		gt4_wordmap_delete (map);
		return NULL;
	}
	map->wordlist = cdata + sizeof (GT4ListHeader);
	if (scout) {
		scout_mmap ((const unsigned char *) cdata, csize);
	}
	return map;
}

void
gt4_wordmap_release (GT4WordMap *map)
{
	if (map->filename) {
		free (map->filename);
		map->filename = NULL;
	}
	if (map->file_map) {
		munmap ((void *) map->file_map, map->file_size);
		map->file_map = NULL;
		map->file_size = 0;
	}
	map->header = NULL;
	map->wordlist = NULL;
	map->user_data = NULL;
}

void
gt4_wordmap_delete (GT4WordMap *map)
{
	gt4_wordmap_release (map);
	free (map);
}

unsigned int 
wordmap_search_query (GT4WordMap *map, unsigned long long query, parameters *p, int printall, unsigned int equalmmonly, unsigned int dosubtraction, GT4WordMap *querymap)
{
	static wordtable mm_table = {0};
	unsigned long long i, nwords = 0L;
	unsigned int count = 0L, currentcount = 0L, querycount = 0L;

	/* if no mismatches */
	if (!p->nmm) {
		return gt4_wordmap_lookup (map, query);
	}

	mm_table.wordlength = p->wordlength;

	/* find and set table size */
	if (!mm_table.nwords) {
		nwords = generate_mismatches (NULL, query, p->wordlength, 0, p->nmm, p->pm3, 0, 1, 0);
		wordtable_ensure_size (&mm_table, nwords, 0);
		if (debug_wordmap > 1) {
			fprintf (stderr, "MM Table size %llu, num mismatches %llu\n", mm_table.nwords, nwords);
		}
	}
	generate_mismatches (&mm_table, query, p->wordlength, 0, p->nmm, p->pm3, 0, 0, equalmmonly);
	if (debug_wordmap > 1) {
		fprintf (stderr, "MM Table size %llu\n", mm_table.nwords);
	}

	for (i = 0; i < mm_table.nwords; i++) {
		if (dosubtraction) {
			querycount = gt4_wordmap_lookup (querymap, mm_table.words[i]);
			currentcount = gt4_wordmap_lookup (map, mm_table.words[i]);
			if (currentcount > querycount) {
				if (debug_wordmap > 1) {
					fprintf (stderr, "%llu %llu %llu querycount %u currentcount %u\n", query, i, mm_table.words[i], querycount, currentcount);
				}
				mm_table.nwords = 0;
				return ~0L;
			}
			count += (currentcount - querycount);
		} else {
			currentcount = gt4_wordmap_lookup (map, mm_table.words[i]);
			count += currentcount;
			if (printall && currentcount > 0) {
				fprintf (stdout, "%s\t%u\n", word_to_string (mm_table.words[i], mm_table.wordlength), currentcount);
			}
		}
	}
	
	mm_table.nwords = 0;
	return count;
}

unsigned int 
gt4_wordmap_lookup_canonical (GT4WordMap *map, unsigned long long query)
{
	unsigned long long word, low, high, mid;
	low = 0;
	high = map->header->nwords - 1;
	mid = (low + high) / 2;
			
	while (low <= high) {
		word = WORDMAP_WORD (map, mid);
		if (word < query) {
			low = mid + 1;
		} else if (word > query) {
			if (mid == 0) break;
			high = mid - 1;
		} else {
			return WORDMAP_FREQ (map, mid);
		}
		mid = (low + high) / 2;
	}
	return 0;
}

unsigned int 
gt4_wordmap_lookup (GT4WordMap *map, unsigned long long query)
{
	unsigned long long rev = get_reverse_complement (query, map->header->wordlength);
	if (rev < query) query = rev;
	return gt4_wordmap_lookup_canonical (map, query);
}




