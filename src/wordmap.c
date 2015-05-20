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

#include "wordmap.h"
#include "wordtable.h"
#include "sequence.h"
#include "common.h"
#include "utils.h"

unsigned int debug_wordmap = 0;

unsigned int glistmaker_code_match = 'G' << 24 | 'T' << 16 | '4' << 8 | 'C';


wordmap * 
wordmap_new (const char *listfilename)
{
	const char *content;
	size_t size;
	wordmap *map = (wordmap *) malloc (sizeof (wordmap));
	if (!map) return NULL;
	
	memset (map, 0, sizeof (wordmap));
	map->filename = listfilename;
	
	content = mmap_by_filename (listfilename, &size);
	if (!content) {
		wordmap_delete (map);
		return NULL;
	}
	
	map->header = (header *) content;
	if (map->header->code != glistmaker_code_match || size != sizeof (header) + map->header->nwords * 12) {
		free ((void *) content);
		wordmap_delete (map);
		return NULL;
	}
	map->wordlist = content + sizeof (header);
	return map;
}

void 
wordmap_delete (wordmap *map)
{
	if (map->wordlist) free ((void *) map->wordlist);
	if (map->header) free ((void *) map->header);
	if (map->additional_data) free ((void *) map->additional_data);
	if (map) free ((void *) map);
	return;
}

void 
wordmap_set_additional_data (wordmap *map, void *additional_data)
{
	map->additional_data = additional_data;
}

unsigned int 
wordmap_search_query (wordmap *map, unsigned long long query, parameters *p, int printall, int equalmmonly,
		int dosubtraction, wordmap *querymap)
{
	static wordtable mm_table = {0};
	unsigned long long i, nwords = 0L;
	unsigned int count = 0L, currentcount = 0L, querycount = 0L;

	/* if no mismatches */
	if (!p->nmm) {
		return search_query_from_both_strands (map, query);
	}

	mm_table.wordlength = p->wordlength;

	/* find and set table size */
	if (mm_table.nwordslots == 0) {
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
			querycount = search_query_from_both_strands (querymap, mm_table.words[i]);
			currentcount = search_query_from_both_strands (map, mm_table.words[i]);
			if (currentcount > querycount) {
				if (debug_wordmap > 1) {
					fprintf (stderr, "%llu %llu %llu querycount %u currentcount %u\n", query, i, mm_table.words[i], querycount, currentcount);
				}
				mm_table.nwords = 0;
				return ~0L;
			}
			count += (currentcount - querycount);

		} else {
			currentcount = search_query_from_both_strands (map, mm_table.words[i]);
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
binary_search (wordmap *map, unsigned long long query)
{
	unsigned long long word, low, high, mid;
	unsigned int freq;
	low = 0;
	high = map->header->nwords - 1;
	mid = (low + high) / 2;
			
	while (low <= high) {
		word = *((unsigned long long *) (map->wordlist + mid * (sizeof (unsigned long long) + sizeof (unsigned int))));
				
		if (word < query) {
			low = mid + 1;
		} else if (word > query) {
			if (mid == 0) break;
			high = mid - 1;
		} else {
			freq = *((unsigned int *) (map->wordlist + mid * (sizeof (unsigned long long) + sizeof (unsigned int)) + sizeof (unsigned long long)));
			return freq;
		}
		mid = (low + high) / 2;
	}
	return 0;
}

unsigned int 
search_query_from_both_strands (wordmap *map, unsigned long long query)
{
	unsigned int freqfw, freqrv;
	freqfw = binary_search (map, query);
	if (!freqfw) {
		freqrv = binary_search (map, get_reverse_complement (query, map->header->wordlength));
		return freqrv;
	}
	return freqfw;
}




