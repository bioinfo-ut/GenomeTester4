#ifndef __WORDMAP_H__
#define __WORDMAP_H__

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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * GT4WordMap is the most basic list container
 */

#ifndef __WORDMAP_C__
/* Defaults to 0, can be increased to print debug information */
extern unsigned int debug_wordmap;
/* List tag "GT4C" encoded to big-endian 32-bit integer */
extern unsigned GT4_LIST_CODE;
#endif

#define WORDMAP_ELEMENT_SIZE (sizeof (unsigned long long) + sizeof (unsigned int))

typedef struct _GT4ListHeader GT4ListHeader;
typedef struct _GT4WordMap GT4WordMap;

struct _GT4ListHeader {
	unsigned int code;
	unsigned int version_major;
	unsigned int version_minor;
	unsigned int wordlength;
	unsigned long long nwords;
	unsigned long long totalfreq;
	unsigned long long padding;
};

struct _GT4WordMap {
	char *filename;
	const unsigned char *file_map;
	unsigned long long file_size;
	GT4ListHeader *header;
	const unsigned char *wordlist;
	void *user_data;
};

typedef struct _parameters {
	unsigned int wordlength;
	unsigned int nmm;
	unsigned int pm3;
	double coef;
} parameters;

#define WORDMAP_WORD(w,i) (*((unsigned long long *) ((w)->wordlist + 12 * (i))))
#define WORDMAP_FREQ(w,i) (*((unsigned int *) ((w)->wordlist + 12 * (i) + 8)))

/* Creates new GT4WordMap by memory-mapping file, returns NULL if error */
/* If "scout" is true, a new thread is created that sequentially prefetces the map into virtual memory */
GT4WordMap *gt4_wordmap_new (const char *listfilename, unsigned int scout);
/* Releases allocated and mapped memory and cleans data fields */
void gt4_wordmap_release (GT4WordMap *map);
/* Releases wordmap and frees the structure */
void gt4_wordmap_delete (GT4WordMap *map);

unsigned int wordmap_search_query (GT4WordMap *wmap, unsigned long long query, parameters *p, int printall, unsigned int equalmmonly, unsigned int dosubtraction, GT4WordMap *querymap);

unsigned int gt4_wordmap_lookup_canonical (GT4WordMap *wmap, unsigned long long query);
unsigned int gt4_wordmap_lookup (GT4WordMap *wmap, unsigned long long query);

#endif /* WORDMAP_H_ */
