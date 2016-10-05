#ifndef __WORDMAP_H__
#define __WORDMAP_H__

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


#ifndef __WORDMAP_C__
extern unsigned int debug_wordmap;
extern unsigned glistmaker_code_match;
#endif

#define WORDMAP_ELEMENT_SIZE (sizeof (unsigned long long) + sizeof (unsigned int))

typedef struct _header {
	unsigned int code;
	unsigned int version_major;
	unsigned int version_minor;
	unsigned int wordlength;
	unsigned long long nwords;
	unsigned long long totalfreq;
	unsigned long long padding;
} header;

typedef struct _parameters {
	unsigned int wordlength;
	unsigned int nmm;
	unsigned int pm3;
	double coef;
} parameters;

typedef struct _wordmap {
	const char *filename;
	unsigned char *file_map;
	unsigned long long file_size;
	header *header;
	const char *wordlist;
	void *additional_data;
} wordmap;

#define WORDMAP_WORD(w,i) (*((unsigned long long *) ((w)->wordlist + 12 * (i))))
#define WORDMAP_FREQ(w,i) (*((unsigned int *) ((w)->wordlist + 12 * (i) + 8)))

wordmap *wordmap_new (const char *listfilename, unsigned int scout);
void wordmap_release (wordmap *map);

void wordmap_delete (wordmap *map);

void wordmap_set_additional_data (wordmap *map, void *data);

unsigned int wordmap_search_query (wordmap *wmap, unsigned long long query, parameters *p,
		int printall, int equalmmonly, int subtract, wordmap *querywmap);

unsigned int search_query_from_both_strands (wordmap *wmap, unsigned long long query);

unsigned int binary_search (wordmap *wmap, unsigned long long query);

#endif /* WORDMAP_H_ */
