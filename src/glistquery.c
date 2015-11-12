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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "utils.h"
#include "sequence.h"
#include "wordtable.h"
#include "wordmap.h"
#include "fasta.h"
#include "common.h"

typedef struct _querystructure {
	wordmap *map;
	parameters *p;
	unsigned int minfreq;
	unsigned int maxfreq;
	int printall;
} querystructure;

void search_one_query_string (wordmap *map, const char *querystring, parameters *p, unsigned int minfreq, unsigned int maxfreq, int printall);
int search_n_query_strings (wordmap *map, const char *queryfile, parameters *p, unsigned int minfreq, unsigned int maxfreq, int printall);
int search_fasta (wordmap *map, const char *seqfilename, parameters *p, unsigned int minfreq, unsigned int maxfreq, int printall);
int search_list (wordmap *map, const char *querylistfilename, parameters *p, unsigned int minfreq, unsigned int maxfreq, int printall);
int process_word (FastaReader *reader, unsigned long long word, void *data);
int print_full_map (wordmap *map);
void get_statistics (wordmap *map);
void print_help (int exitvalue);

int debug = 0;

int main (int argc, const char *argv[])
{

	int argidx, v = 0;
	const char *listfilename = NULL;
	const char *querystring = NULL, *queryfilename = NULL, *seqfilename = NULL, *querylistfilename = NULL;
	parameters p = {0};
	wordmap *map;
	char *end;
	int printall = 0, getstat = 0;
	unsigned int minfreq = 0, maxfreq = UINT_MAX;

	/* parsing commandline */
	
	/* FIXME: null argumenti!!! */
	
	
	for (argidx = 1; argidx < argc; argidx++) {

		if (!strcmp (argv[argidx], "-v") || !strcmp (argv[argidx], "--version")) {
			fprintf (stdout, "glistquery v%d.%d\n", VERSION_MAJOR, VERSION_MINOR);
			return 0;

		} else if (!strcmp (argv[argidx], "-h") || !strcmp (argv[argidx], "--help") || !strcmp (argv[argidx], "-?")) {
			print_help (0);

		} else if (argidx == 1) {
			/* get the locations of the input files */
			if (argv[argidx][0] == '-') {
				fprintf(stderr, "Error: No list file specified!\n");
				print_help (0);
			}
			listfilename = argv[argidx];
		} else if (!strcmp(argv[argidx], "-s") || !strcmp(argv[argidx], "--seqfile")) {
			if (!argv[argidx + 1] || argv[argidx + 1][0] == '-') {
				fprintf(stderr, "Warning: No sequence file name specified!\n");
				argidx += 1;
				continue;
			}
			seqfilename = argv[argidx + 1];
			argidx += 1;
		} else if (!strcmp(argv[argidx], "-l") || !strcmp(argv[argidx], "--listfile")) {
			if (!argv[argidx + 1] || argv[argidx + 1][0] == '-') {
				fprintf(stderr, "Warning: No query list file name specified!\n");
				argidx += 1;
				continue;
			}
			querylistfilename = argv[argidx + 1];
			argidx += 1;
		} else if (!strcmp(argv[argidx], "-f") || !strcmp(argv[argidx], "--queryfile")) {
			if (!argv[argidx + 1] || argv[argidx + 1][0] == '-') {
				fprintf(stderr, "Warning: No query file name specified!\n");
				argidx += 1;
				continue;
			}
			queryfilename = argv[argidx + 1];
			argidx += 1;
		} else if (!strcmp(argv[argidx], "-q") || !strcmp(argv[argidx], "--query")) {
			if (!argv[argidx + 1] || argv[argidx + 1][0] == '-') {
				fprintf(stderr, "Warning: No query specified!\n");
				argidx += 1;
				continue;
			}
			querystring = argv[argidx + 1];
			argidx += 1;
		} else if (!strcmp(argv[argidx], "-p") || !strcmp(argv[argidx], "--perfectmatch")) {
			if (!argv[argidx + 1]) {
				fprintf(stderr, "Warning: No number of 3 prime perfect matches specified! Using the default value: %d.\n", p.pm3);
				argidx += 1;
				continue;
			}
			p.pm3 = strtol (argv[argidx + 1], &end, 10);
			if (*end != 0) {
				fprintf(stderr, "Error: Invalid number of 3 prime perfect matches: %s! Must be an integer.\n", argv[argidx + 1]);
				print_help (1);
			}
			argidx += 1;
		} else if (!strcmp(argv[argidx], "-mm") || !strcmp(argv[argidx], "--mismatch")) {
			if (!argv[argidx + 1]) {
				fprintf(stderr, "Warning: No number of mismatches specified! Using the default value: %d.\n", p.nmm);
				argidx += 1;
				continue;
			}
			p.nmm = strtol (argv[argidx + 1], &end, 10);
			if (*end != 0) {
				fprintf(stderr, "Error: Invalid number of mismatches: %s! Must be an integer.\n", argv[argidx + 1]);
				print_help (1);
			}
			argidx += 1;
		} else if (!strcmp(argv[argidx], "-min") || !strcmp(argv[argidx], "--minfreq")) { 
			if (!argv[argidx + 1]) {
				fprintf(stderr, "Warning: No minimum frequency specified! Using the default value: %d.\n", minfreq);
				argidx += 1;
				continue;
			}
			minfreq = strtol (argv[argidx + 1], &end, 10);
			if (*end != 0) {
				fprintf(stderr, "Error: Invalid minimum frequency: %s! Must be a positive integer.\n", argv[argidx + 1]);
				print_help (1);
			}
			argidx += 1;
		} else if (!strcmp(argv[argidx], "-max") || !strcmp(argv[argidx], "--maxfreq")) {
			if (!argv[argidx + 1]) {
				fprintf(stderr, "Warning: No maximum frequency specified! Using the default value: %d.\n", maxfreq);
				argidx += 1;
				continue;
			}
			maxfreq = strtol (argv[argidx + 1], &end, 10);
			if (*end != 0) {
				fprintf(stderr, "Error: Invalid maximum frequency: %s! Must be a positive integer.\n", argv[argidx + 1]);
				print_help (1);
			}
			argidx += 1;
		} else if (!strcmp(argv[argidx], "-D")) {
			debug += 1;
		} else if (!strcmp(argv[argidx], "-all")) {
			printall = 1;
		} else if (!strcmp(argv[argidx], "-stat")) {	
			getstat = 1;
		} else {
			fprintf(stderr, "Error: Unknown argument: %s!\n", argv[argidx]);
			print_help (1);
		}
	}
	
	if (!listfilename) {
		fprintf(stderr, "Error: Missing a list file!\n");
		print_help (1);			
	}

	/* checking the parameters */
	if (p.nmm < 0) {
		fprintf(stderr, "Error: Invalid number of mismatches: %d! Must be between 0 and the length of the query.\n", p.nmm);
		print_help (1);
	}
	if (p.pm3 < 0) {
		fprintf(stderr, "Error: Invalid number of 3 prime perfect matches: %d! Must be between 0 and the length of the query.\n", p.pm3);
		print_help (1);
	}
	if (minfreq < 0) {
		fprintf(stderr, "Error: Invalid number of minimum frequency: %d! Must be positive integer.\n", minfreq);
		print_help (1);	  
	}
	if (maxfreq < 0) {
		fprintf(stderr, "Error: Invalid number of maximum frequency: %d! Must be positive integer.\n", maxfreq);
		print_help (1);	  
	}

	map = wordmap_new (listfilename, !getstat);
	if (!map) {
		fprintf (stderr, "Error: Could not make wordmap from file %s!\n", listfilename);
		return 1;
	}
	
	if (map->header->version_major > VERSION_MAJOR || map->header->version_minor > VERSION_MINOR) {
		fprintf (stderr, "Error: %s is created with a newer glistmaker version.", map->filename);
		exit (1);
	}

	if (map->header->code != glistmaker_code_match) {
		fprintf (stderr, "Error: %s is not a glistmaker v.4 file.\n", map->filename);
		exit (1);
	}
	
	if (getstat) {
		get_statistics (map);
		exit (0);
	}

	p.wordlength = map->header->wordlength;

	/* glistquery options */
	if (seqfilename) { /* fasta input */

		search_fasta (map, seqfilename, &p, minfreq, maxfreq, printall);

	} else if (querylistfilename) { /* list input */	
		
		search_list (map, querylistfilename, &p, minfreq, maxfreq, printall);
	  
	} else if (queryfilename) { /* list of queries */

		v = search_n_query_strings (map, queryfilename, &p, minfreq, maxfreq, printall);
		if (v) return v;

	} else if (querystring) { /* one query */

		/* checking possible errors */
		if (p.wordlength != strlen (querystring)) {
			fprintf (stderr, "Error: Incompatible wordlengths! Wordlength in list: %u, query length: %lu\n", p.wordlength, strlen (querystring));
			return 1;
		}
		if (p.wordlength - p.pm3 < p.nmm) {
			fprintf(stderr, "Error: Number or mismatches specified is too large for %s with %d nucleotides long 3 prime perfect match.\n", querystring, p.pm3);
			return 1;
		}
		search_one_query_string (map, querystring, &p, 0, UINT_MAX, printall);

	} else { /* no query */

		print_full_map (map);
	}
	exit (0);
}

/* print the whole list */
int print_full_map (wordmap *map)
{
	unsigned long long i, word;
	unsigned int freq;
	char *word_str;

	for (i = 0; i < map->header->nwords; i++) {
		word = *((unsigned long long *) (map->wordlist + i * (sizeof (unsigned long long) + sizeof (unsigned))));
		freq = *((unsigned *) (map->wordlist + i * (sizeof (unsigned long long) + sizeof (unsigned)) + sizeof (unsigned long long)));
		word_str = word_to_string (word, map->header->wordlength);
		fprintf (stdout, "%s\t%u\n", word_str, freq);
	}
	fprintf (stdout, "NUnique\t%llu\nNTotal\t%lld\n", map->header->nwords, map->header->totalfreq);
	return 0;
}

/* use one specific query sequence */
void search_one_query_string (wordmap *map, const char *querystring, parameters *p, unsigned int minfreq, unsigned int maxfreq, int printall)
{
	unsigned int freq = 0;
	unsigned long long query;

	query = string_to_word (querystring, p->wordlength);
	freq = wordmap_search_query (map, query, p, printall, 0, 0, NULL);
	if (!printall && freq >= minfreq && freq <= maxfreq) fprintf (stdout, "%s\t%u\n", querystring, freq);
	return;
}

/* use file which consist of many queries */
int search_n_query_strings (wordmap *map, const char *queryfile, parameters *p, unsigned int minfreq, unsigned int maxfreq, int printall)
{
	FILE *f;
	char querystring[256];
	int beg = 1;

	f = fopen (queryfile, "r");
	if (f == NULL) {
		fprintf (stderr, "Error: Cannot open file %s.\n", queryfile);
		return 1;
	}

	while (fscanf (f, "%s\n", querystring) != EOF) {
		if (beg) {
			/* checking possible errors */
			if (p->wordlength != strlen (querystring)) {
				fprintf (stderr, "Error: Incompatible wordlengths! Wordlength in list: %u, query length: %lu\n", p->wordlength, strlen (querystring));
				return 1;
			}
			if (p->wordlength - p->pm3 < p->nmm) {
				fprintf (stderr, "Error: Number or mismatches specified is too large for %s with %d nucleotides long 3 prime perfect match.\n", querystring, p->pm3);
				return 1;
			}
			beg = 0;
		}
		search_one_query_string (map, querystring, p, minfreq, maxfreq, printall);
	}
	return 0;
}

int search_fasta (wordmap *map, const char *seqfilename, parameters *p, unsigned int minfreq, unsigned int maxfreq, int printall)
{
	const char *ff;
	size_t size;
	querystructure qs = {0};
	int v;
	FastaReader r;

	qs.map = map;
	qs.p = p;
	qs.minfreq = minfreq;
	qs.maxfreq = maxfreq;
	qs.printall = printall;
	ff = mmap_by_filename (seqfilename, &size);

	fasta_reader_init_from_data (&r, p->wordlength, 0, (const unsigned char *) ff, size);
	v = fasta_reader_read_nwords (&r, size, NULL, NULL, NULL, NULL, process_word, (void *) &qs);
	/* v = fasta_reader_read_nwords (ff, size, p->wordlength, (void *) &qs, 0, 0, NULL, NULL, process_word); */
	fasta_reader_release (&r);
	if (v) return v;
	return 0;
}

int search_list (wordmap *map, const char *querylistfilename, parameters *p, unsigned int minfreq, unsigned int maxfreq, int printall)
{
	wordmap *qmap;
	unsigned long long i, word = 0L;
	unsigned int freq;
	
	qmap = wordmap_new (querylistfilename, 1);
	
	if (map->header->wordlength != qmap->header->wordlength) return GT_INCOMPATIBLE_WORDLENGTH_ERROR;
	
	for (i = 0; i < qmap->header->nwords; i++) {
		word = *((unsigned long long *) (qmap->wordlist + i * (sizeof (unsigned long long) + sizeof (unsigned))));
		freq = wordmap_search_query (map, word, p, printall, 0, 0, NULL);
		if (!printall && freq >= minfreq && freq <= maxfreq) fprintf (stdout, "%s\t%u\n", word_to_string (word, map->header->wordlength), freq);
	}
	return 0;
}

int process_word (FastaReader *reader, unsigned long long word, void *data)
{
	unsigned int freq = 0;
	querystructure *qs = (querystructure *) data;
	freq = wordmap_search_query (qs->map, word, qs->p, qs->printall, 0, 0, NULL);
	if (!qs->printall && freq >= qs->minfreq && freq <= qs->maxfreq) fprintf (stdout, "%s\t%u\n", word_to_string (word, reader->wordlength), freq);
	return 0;
}

void get_statistics (wordmap *map)
{
	fprintf (stdout, "Statistics of %s <<Built with glistmaker version %d.%d>>\n", map->filename, map->header->version_major, map->header->version_minor);
	fprintf (stdout, "Wordlength\t%u\n", map->header->wordlength);
	fprintf (stdout, "NUnique\t%llu\n", map->header->nwords);
	fprintf (stdout, "NTotal\t%llu\n", map->header->totalfreq);
	return;
}

void print_help (int exit_value)
{
	fprintf (stderr, "Usage: glistquery <INPUTLIST> [OPTIONS]\n");
	fprintf (stderr, "Options:\n");
	fprintf (stderr, "    -v, --version             - print version information and exit\n");
	fprintf (stderr, "    -h, --help                - print this usage screen and exit\n");
	fprintf (stderr, "    -stat                     - print statistics of the list file and exit\n");
	fprintf (stderr, "    -q, --query               - single query word\n");
	fprintf (stderr, "    -f, --queryfile           - list of query words in a file\n");
	fprintf (stderr, "    -s, --seqfile             - FastA/FastQ file\n");
	fprintf (stderr, "    -l, --listfile            - list file made by glistmaker\n");
	fprintf (stderr, "    -mm, --mismatch NUMBER    - specify number of mismatches (default 0, can be used with -d and -dd)\n");
	fprintf (stderr, "    -p, --perfectmatch NUMBER - specify number of 3' perfect matches (default 0)\n");
	fprintf (stderr, "    -min, --minfreq NUMBER    - minimum frequency of the printed words (default 0)\n");
	fprintf (stderr, "    -max, --maxfreq NUMBER    - maximum frequency of the printed words (default MAX_UINT)\n");
	fprintf (stderr, "    -all                      - in case of mismatches prints all found words\n");
	fprintf (stderr, "    -D                        - increase debug level\n");
	exit (exit_value);
}


